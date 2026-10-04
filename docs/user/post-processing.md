# Post-processing

The scripts live in `scripts/` and are copied into every case. Run them from
inside the case directory, or pass the case path.

## The full figure set in one command

```bash
python3 make_paper_figures.py .                 # everything
python3 make_paper_figures.py . --skip rupture  # skip the per-event figures
python3 make_paper_figures.py . --only catalog figure4 analysis
python3 make_paper_figures.py --list            # stages
```

Each stage logs to `aPlots/logs/<stage>.log`. A stage that cannot apply to a
case (the San Andreas paleoseismic sites, say) fails on its own without
stopping the others.

## Individual tools

| tool | what it makes |
|---|---|
| `plotRuptureDynamics` | per-event figure: shear, normal stress, slip, rupture time along strike. Saves events with Mw ≥ `MIN_PLOT_MAGNITUDE` (default 6.5); `FORCE_REPLOT=1` redraws |
| `CATALOG=1 plotRuptureDynamics` | `aPlots/catalog.csv` only, about 10× faster |
| `analyze_catalog.py . --mmax 7.0` | b-value, magnitude-frequency, magnitude and nucleation through time |
| `plot_event_slips_overtime_fig4.py .` | slip-distribution stacks through time (Figure 4 style); `--tstart`, `--duration` in kyr, `--threshold` minimum slip |
| `plot_saf_figure3.py` | long-term slip rate against geological rates |
| `plot_saf_figure6.py` | cumulative moment and magnitude-frequency |
| `plot_saf_figure9.py` | characteristic-event slip distributions, chosen by rupture footprint |
| `paleo_site_stats.py` | recurrence and slip statistics at paleoseismic sites |
| `compare_cycle_over_strike.py` | one cycle from several cases overlaid |
| `plot_stf.py . <eqid>` | rupture time histories of one event (needs `outputSTF`) |

The SAF-specific tools apply to `paper.saf.A` and `saf.gmsh.lite`.

`plot_event_slips_overtime_fig4.py` exits with an error when a window holds
no event above the slip threshold, rather than saving an empty figure. Event
times start at the first event, so the last event sits earlier than the total
simulated time.

## Comparing with the published results

```bash
bash fetch_published_reference.sh
```

downloads the published Liu et al. (2022) software (Zenodo
10.5281/zenodo.5823021) and results (PANGAEA 10.1594/PANGAEA.940262) by DOI,
with checksums. The figure tools can then read the published output for
side-by-side comparison.
