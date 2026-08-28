#!/usr/bin/env bash
# make_results_bundle.sh -- pack a finished run into a small, self-describing
# archive suitable for a GitHub release asset or a Zenodo upload.
#
# What goes in: the derived catalogue, the recurrence intervals, the complete
# model inputs (mesh + loading + FE files), the summary figures, and a manifest
# recording binary version, thread count and cycle range.
#
# What stays out: totalop.txt*, which is ~99% of the bytes (0.9-1.5 GB per run,
# ~4.6x gzip). Ship that separately with a DOI if it is needed; the catalogue
# is what analysis actually consumes.
#
# Usage: bash scripts/make_results_bundle.sh <case_dir> [out_dir]

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CASE="$(cd "$1" && pwd)"; NAME="$(basename "$CASE")"
OUT="${2:-$ROOT/work/bundles}"; mkdir -p "$OUT"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
B="$TMP/$NAME"; mkdir -p "$B"

# Collect inputs and small outputs. Patterns are matched with find inside each
# candidate directory -- a bare glob in a `for` list would expand against the
# caller's cwd, not the case dir, and silently copy nothing.
for d in "$CASE" "$CASE/aRawSimuData" "$CASE/fem_mesh_output"; do
    [ -d "$d" ] || continue
    find "$d" -maxdepth 1 -type f \( \
        -name 'interval.txt*' -o -name 'cyclelog.txt*' -o -name 'meshGeneralInfo.txt' -o \
        -name 'nsmpGeoPhys.txt' -o -name 'nsmp.txt' -o -name 'vert.txt' -o \
        -name 'fac.txt' -o -name 'nsmpTanLen.txt' -o -name 'FE_*.txt' -o \
        -name 'Rate_direction.txt' -o -name 'x[0-9]_1.txt' -o \
        -name 'user_defined_params.py' -o -name 'run_*.log' \) \
        -exec cp -n {} "$B/" \; 2>/dev/null || true
done

# Figures and the catalogue may live in the case dir or in its snapshot
# (snapshot_case.py writes to <case>_snap, and for a restarted run that is the
# only place the merged, whole-sequence catalogue exists).
mkdir -p "$B/aPlots"
for a in "$CASE/aPlots" "${CASE}_snap/aPlots"; do
    [ -d "$a" ] || continue
    find "$a" -maxdepth 1 -type f \( -name '*.csv' -o -name 'catalog_analysis.png' \
        -o -name 'Figure*.png' -o -name 'faults.png' -o -name 'loading_inputs.png' \) \
        -exec cp {} "$B/aPlots/" \; 2>/dev/null || true
done
# Prefer the snapshot's merged interval/cyclelog when a restart split them.
if [ -f "${CASE}_snap/interval.txt1" ]; then
    cp "${CASE}_snap/interval.txt1" "$B/interval.merged.txt" 2>/dev/null || true
fi
rmdir "$B/aPlots" 2>/dev/null || true

IVS=$(ls "$CASE"/interval.txt* "$CASE"/aRawSimuData/interval.txt* 2>/dev/null || true)
cyc=$( { [ -n "$IVS" ] && cat $IVS; } 2>/dev/null | wc -l | tr -d ' ')
yrs=$( { [ -n "$IVS" ] && cat $IVS; } 2>/dev/null | awk '{s+=$1}END{printf "%.0f", s}')
{
  echo "# $NAME"
  echo
  echo "| | |"
  echo "|---|---|"
  echo "| Cycles | $cyc |"
  echo "| Simulated time | ${yrs} yr |"
  echo "| Code version | $(cat "$ROOT/VERSION") |"
  echo "| Binary | $(ls "$CASE" | grep -m1 run_eqdyna2d || echo unknown) |"
  echo "| Threads | ${OMP_NUM_THREADS:-recorded in run log; results are thread-count dependent, see PROJECT_RULES R30} |"
  echo "| Packed | $(date -u '+%Y-%m-%dT%H:%M:%SZ') |"
  echo
  echo "\`totalop.txt*\` is deliberately excluded (see make_results_bundle.sh)."
  echo "Everything needed to re-run is here: restore the inputs into a case"
  echo "directory, build the matching version, and run at the same thread count."
} > "$B/MANIFEST.md"

tar czf "$OUT/$NAME.tar.gz" -C "$TMP" "$NAME"
echo "$OUT/$NAME.tar.gz  $(du -h "$OUT/$NAME.tar.gz" | cut -f1)"
