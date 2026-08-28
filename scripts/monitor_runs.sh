#!/usr/bin/env bash
# monitor_runs.sh -- periodically snapshot running cases and re-plot them.
#
# Each pass calls snapshot_case.py, which truncates every segment to whole
# cycles and merges restart segments (totalop.txt1 + totalop.txt73 + ...) into
# one continuous sequence, then runs the plot suite on the snapshot. Plotting
# the live directory instead would catch half-written cycles and would silently
# report only one restart segment.
#
# Usage:
#   bash scripts/monitor_runs.sh                          # defaults below, 20 min
#   bash scripts/monitor_runs.sh -i 600 work/a work/b     # custom interval + cases
#   bash scripts/monitor_runs.sh -1 work/a                # single pass, then exit
#   nohup bash scripts/monitor_runs.sh > monitor.log 2>&1 &

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

INTERVAL=1200          # 20 minutes
ONCE=0
MINMAG="${MIN_PLOT_MAGNITUDE:-6.0}"
STAGES="catalog figure4 analysis rupture"

while [ $# -gt 0 ]; do
    case "$1" in
        -i) INTERVAL="$2"; shift 2 ;;
        -1) ONCE=1; shift ;;
        -m) MINMAG="$2"; shift 2 ;;
        -s) STAGES="$2"; shift 2 ;;
        *)  break ;;
    esac
done

CASES=("$@")
if [ ${#CASES[@]} -eq 0 ]; then
    for d in work/*/; do
        d="${d%/}"
        case "$d" in *_snap) continue ;; esac
        [ -f "$d/totalop.txt1" ] && CASES+=("$d")
    done
fi
[ ${#CASES[@]} -eq 0 ] && { echo "no cases with totalop.txt1 under work/"; exit 1; }

echo "monitoring: ${CASES[*]}"
echo "interval ${INTERVAL}s, stages: $STAGES, M>=$MINMAG"

while true; do
    echo "================ $(date '+%Y-%m-%d %H:%M:%S') ================"
    for d in "${CASES[@]}"; do
        [ -d "$d" ] || continue
        snap="${d}_snap"
        running=$(pgrep -f "$(basename "$d")" >/dev/null 2>&1 && echo yes || echo "no")
        if ! python3 scripts/snapshot_case.py "$d" "$snap" --quiet; then
            echo "[$d] snapshot failed"; continue
        fi
        if python3 scripts/make_paper_figures.py "$snap" \
               --only $STAGES --min-magnitude "$MINMAG" >/dev/null 2>&1; then
            n=$(( $(wc -l < "$snap/aPlots/catalog.csv") - 1 ))
            big=$(awk -F, 'NR>1 && $2>=6.5' "$snap/aPlots/catalog.csv" | wc -l)
            mx=$(awk -F, 'NR>1{if($2>m)m=$2}END{printf "%.2f", m}' "$snap/aPlots/catalog.csv")
            yr=$(awk '{s+=$1}END{printf "%.0f", s}' "$snap/interval.txt1")
            echo "[$d] $n events over ${yr} yr, Mmax $mx, M>=6.5: $big  -> $snap/aPlots/"
        else
            echo "[$d] plotting failed; see $snap/aPlots/logs/"
        fi
    done
    [ "$ONCE" -eq 1 ] && break
    echo "sleeping ${INTERVAL}s..."
    sleep "$INTERVAL"
done
