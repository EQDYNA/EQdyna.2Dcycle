# EQdyna.2Dcycle User Guide

EQdyna.2Dcycle is a 2D finite-element code for physics-based multicycle
earthquake simulation on geometrically complex fault systems. Each earthquake
cycle couples a quasi-static interseismic loading phase with a fully dynamic
rupture, so a single run produces thousands of earthquakes whose sizes,
recurrence and rupture extents emerge from fault geometry, friction and
loading rather than being prescribed.

The solver core is Fortran 90 with OpenMP; case set-up, meshing and
post-processing are Python.

## Where to go

* [Getting started](getting-started.md) — install, build, and run a first case.
* [Running a case](running-a-case.md) — the case workflow, restarts, threads, monitoring.
* [Compsets](compsets.md) — the ready-made fault systems and the two meshing modes.
* [Model inputs](inputs.md) — fault geometry and tectonic loading files.
* [Parameters](parameters.md) — every configurable default, generated from the code.
* [Output files](outputs.md) — what each output contains, column by column.
* [Post-processing](post-processing.md) — catalogues, figures, and the paper figure set.
* [Source time functions](source-time-functions.md) — full rupture time histories for ground-motion work.
* [Testing and reproducibility](testing.md) — the test suite and why thread count matters.
* [Troubleshooting](troubleshooting.md) — symptoms and what they mean.
* [Citing](citing.md) — how to cite EQdyna.2Dcycle.

## Getting the code

EQdyna.2Dcycle is released under the MIT License at
[github.com/EQDYNA/EQdyna.2Dcycle](https://github.com/EQDYNA/EQdyna.2Dcycle).
