#!/usr/bin/env bash
# monitor_runs.sh — periodically print run status and re-plot figure 4 for all
# active paper.saf.A.* and saf.gmsh.lite-derived cases under work/.
#
# Usage:
#   bash scripts/monitor_runs.sh                 # poll every 600s (10 min)
#   bash scripts/monitor_runs.sh 300             # custom interval (seconds)
#   nohup bash scripts/monitor_runs.sh > monitor.log 2>&1 &   # background

set -u
INTERVAL="${1:-600}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

while true; do
    ts="$(date '+%Y-%m-%d %H:%M:%S')"
    echo "================ $ts ================"

    for d in work/paper.saf.A.omp* work/paper.saf.A.demo work/paper.saf.A.test \
             work/saflite work/saf.gmsh.lite*; do
        [ -d "$d" ] || continue
        nint="$(wc -l < "$d/interval.txt1" 2>/dev/null || echo 0)"
        nint="${nint// /}"
        bytes="$(wc -c < "$d/totalop.txt1" 2>/dev/null || echo 0)"
        bytes="${bytes// /}"
        last_step="$(grep '^\[' "$d"/run_*.log 2>/dev/null | tail -1 | awk '{print $3,$4,$5}')"
        echo "[$d]  cycles=$nint  totalop=${bytes}B  last: ${last_step:-no log yet}"

        if [ -s "$d/totalop.txt1" ]; then
            python3 scripts/plot_event_slips_overtime_fig4.py "$d" --threshold 0 > /dev/null 2>&1 \
                && echo "    -> $d/aPlots/Figure4_7_8_slipdist_*.png" \
                || echo "    plot_saf_figure4 failed (likely missing nsmp.txt or partial output)"
            ( cd "$d" && MIN_PLOT_MAGNITUDE=6.0 python3 "$ROOT/scripts/plotRuptureDynamics" > /dev/null 2>&1 ) \
                && echo "    -> $d/aPlots/cRuptureDynamics*.png (M>6.0)" \
                || echo "    plotRuptureDynamics failed"
        fi
    done

    echo
    echo "Sleeping ${INTERVAL}s..."
    sleep "$INTERVAL"
done
