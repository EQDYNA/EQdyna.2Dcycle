#!/usr/bin/env bash
# Fetch the GSRM v2.1 strain-rate grid and cut it to the Xianshuihe region.
#
# The data is public and citable, so it is fetched rather than vendored:
#   Kreemer, C., G. Blewitt, E.C. Klein (2014), A geodetic plate motion and
#   Global Strain Rate Model, G-cubed 15, 3849-3889, doi:10.1002/2014GC005407
#   Data: https://geodesy.unr.edu/GSRM/    Licence: CC-BY-NC-SA 3.0, GEM
#
# The global file is 88 MB compressed / 529 MB expanded; only the regional
# subset (~240 KB) is kept.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$HERE/strain_rate_input"
OUT="$OUT_DIR/GSRM_v2.1_xianshuihe_region.txt"

LAT_MIN=27.5; LAT_MAX=33.0
LON_MIN=99.0; LON_MAX=104.0
URL="https://geodesy.unr.edu/GSRM/GSRM_strain.txt.Z"
HEADER_LINES=26      # GSRM licence header, ending with the column-name line

mkdir -p "$OUT_DIR"
if [ -f "$OUT" ]; then
  echo "$OUT already present ($(wc -l < "$OUT") cells); delete it to refetch."
  exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Fetching $URL (88 MB) ..."
if command -v curl >/dev/null; then curl -fSL --retry 3 -o "$TMP/g.txt.Z" "$URL"
elif command -v wget >/dev/null; then wget -q --show-progress -O "$TMP/g.txt.Z" "$URL"
else echo "error: need curl or wget" >&2; exit 1; fi

echo "Decompressing ..."
gunzip -c "$TMP/g.txt.Z" > "$TMP/g.txt" 2>/dev/null || uncompress -c "$TMP/g.txt.Z" > "$TMP/g.txt"

echo "Cutting to ${LON_MIN}-${LON_MAX} E / ${LAT_MIN}-${LAT_MAX} N ..."
awk -v hl="$HEADER_LINES" -v a="$LAT_MIN" -v b="$LAT_MAX" -v c="$LON_MIN" -v d="$LON_MAX" \
    'NR>hl && $1>a && $1<b && $2>c && $2<d' "$TMP/g.txt" > "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") cells)"
echo "Columns (units 1e-9/yr): lat lon exx eyy exy vorticity RL-NLC LL-NLC e1 e2 azi_e1"
echo
echo "Next: python3 strain_rate_loading.py"
