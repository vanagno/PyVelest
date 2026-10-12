# Ioannina 2026: minimum 1-D velocity models

Final products of four full PyVelest runs, one per starting model. Method, settings and interpretation
are in [`docs/ioannina_workflow.md`](../../docs/ioannina_workflow.md).

| folder | starting model | final RMS |
|---|---|---|
| `kiratzi87/` (recommended) | `kiratzi87.nd`, linear (TauP) reading | 0.125 s |
| `ionio/` | `ionio.nd`, linear reading | 0.130 s |
| `kiratzi87_step/` | `kiratzi87.nd`, rows read as layer tops | 0.136 s |
| `ionio_step/` | `ionio.nd`, rows read as layer tops | 0.128 s |

Each folder contains:

* `model.mod` / `model.csv`: the final model (VELEST format / CSV with the Vs used, Vp/1.82).
* `stations.sta` / `stations.csv`: P delays and the applied S delays (1.82 × P); JAN is the reference.
* `events.csv`: the final relocation of all 10,362 events, with depth-quality flags. Filter on
  `depth_reliable`.
* `starting_model.dat`: the starting model the run began from.
* `summary.json`: per-run statistics.
* `figures/`: the pipeline's own figures, downscaled.

`comparison/` holds the cross-run figures. The phase-level relocations (`relocated.csv` / `.out`) are not
committed; rerun stage 5 to regenerate them.
