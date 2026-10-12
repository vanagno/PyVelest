# PyVelest for the Ioannina 2026 sequence (PhaseNet → REAL → hypoSA input)

This branch adapts PyVelest to your data and runs it. It contains:

* **Input converters:** your REAL/hypoSA files → `inputs/catalog.pha` + `inputs/stations.csv`; your TauP models
  `kiratzi87.nd` / `ionio.nd` → VELEST starting models.
* **Pipeline fixes:** several of them change results, including bugs that were already in the original code
  (see [§7](#7-what-changed-in-the-code)).
* **Results:** a minimum 1-D model, station delays and relocations for the target area 39.5–39.8°N /
  20.4–20.8°E, from **four starting models**, in `results/ioannina/`.

Every number below comes from the full runs in `results/ioannina/*/summary.json`.

---

## 1. Quick start

```sh
# VELEST: the bundled src/exc/velest is an Intel build (needs libimf.so). If it does not run:
sh src/exc/build_velest.sh                      # -> src/exc/velest_gfortran
#   then set GLOBAL_SETTINGS.velest_exec: src/exc/velest_gfortran in the config

# Inputs (already committed; rebuild only if the REAL/hypoSA files change)
python tools/prepare_real_inputs.py \
    --stations "inputs/My inputs/ioannina/ioan_stations.dat" \
    --phases "inputs/My inputs/ioannina_final.txt" \
    --sa-catalog "inputs/My inputs/ioannina/catalogSA_allday.txt" \
    --sa-phases "inputs/My inputs/ioannina/phaseSA_allday_kiratzi87.txt"

# Starting models (already committed)
python tools/nd_to_velest.py "inputs/My inputs/kiratzi87.nd" --tops -3 0 1 2 3 4 6 8 10 13 16 20 25 \
    --moho 30 --vtop 3.5 --min-increment 0.05 --out inputs/vmodel_kiratzi87.dat           # add --interp step
python tools/nd_to_velest.py "inputs/My inputs/ionio.nd" --tops -3 0 1 2 3 4 6 8 10 13 16 20 25 30 35 \
    --moho 40 --vtop 3.5 --min-increment 0.05 --out inputs/vmodel_ionio.dat               # for *_step.dat

# Run: one config per copy of the repository (all configs write to outputs/ and figures/)
python pyvelest.py                                    # configs/user.yaml  (kiratzi87 start)
python pyvelest.py --config configs/user_ionio.yaml   # ionio start
python pyvelest.py --stages 5                         # redo single stages; earlier outputs are kept
```

A full run is 5 × 50 + 5 × 25 VELEST inversions of 500 events. One inversion takes about 15 s, so a run
takes about 60 min with `n_jobs: 1` and about 20 min with `n_jobs: 4`. The final products land in
`outputs/stage_fn/`:

* `model.mod` / `model.csv`: the minimum 1-D model (`model.csv` includes the Vs that was used).
* `stations.sta` / `stations.csv`: station delays; P delays plus the S delay that was applied.
* `relocated.csv`: phase-level relocations.
* `events.csv`: one row per event, with quality flags.

---

## 2. Inputs and how they were built

| File | Built from | Notes |
|---|---|---|
| `inputs/stations.csv` | `ioan_stations.dat` | 16 stations, elevation km → m. Only 13 have picks; KBN, LSK and SRN have none. |
| `inputs/catalog.pha` | hypoSA hypocentres and origin times (`catalogSA_allday.txt`) + picks and PhaseNet probabilities (`ioannina_final.txt`) | Events are matched by order and checked against `phaseSA_allday_kiratzi87.txt`. Travel times are re-referenced to the hypoSA origin (they match hypoSA's to 1 ms). The hypoSA residual goes into the header's RMS field. Magnitudes come from `ioannina_final.txt`; its picks are identical to `phase_allday.txt`. |
| `inputs/vmodel_*.dat` | `kiratzi87.nd`, `ionio.nd` | See §3. |

Facts about the data that drive the settings:

* **Catalogue size:** 10,362 events and 130,843 picks (69,153 P / 61,690 S). 9,540 events (92 %) lie in
  the target box.
* **Network geometry:** no station lies inside the target box. The nearest are JAN (~21 km) and IGT
  (~26 km), the farthest with picks is KZN (~124 km). Depth resolution is the weak point (§6).
* **Pinned depths:** **hypoSA pins 44 % of the events at exactly 0 km.** These events are left out of the
  stage-1 event selection. In stages 2 and 5 VELEST starts them at 10 km
  (`INPUTS.unconstrained_depth_below_km` / `start_depth_km`), because a relocation started at the free
  surface stays there.
* **Pick weights:** PhaseNet probabilities become VELEST weight classes (`PICK_WEIGHTS`):

  | probability | ≥ 0.85 | 0.70–0.85 | 0.55–0.70 | 0.40–0.55 | < 0.40 |
  |---|---|---|---|---|---|
  | VELEST class | 0 | 1 | 2 | 3 | 4 (not used in inversions) |
  | share of P picks | 43 % | 22 % | 15 % | 13 % | 7 % |
  | share of S picks | 7 % | 24 % | 30 % | 24 % | 15 % |

* **Vp/Vs = 1.82**, from a Wadati fit of 27k P–S pairs in the box; the median per-event slope is 1.83.
  VELEST uses the S picks with this fixed ratio (`phase_mode: 3`), so the S delay of each station is
  1.82 × its P delay.

---

## 3. Starting models

TauP `.nd` files are interpolated **linearly** between the listed depths; a jump needs two rows at the
same depth. `kiratzi87.nd` and `ionio.nd` list their crust as single rows. TauP (and REAL's travel-time
tables) therefore see them as **gradient** models, even if they were meant as layer tops. I built both
readings:

* `linear` = TauP semantics (the default of `tools/nd_to_velest.py`);
* `step` = each row is the top of a constant-velocity layer.

That gives four starting models.

Common choices:

* **Top layer:** −3 to 0 km at 3.5 km/s in all four models. It only covers station elevations (up to
  1.13 km), it is heavily damped and it is never merged. When it differed between models it acted as a
  different velocity floor in each, because low-velocity layers are not allowed.
* **Layering:** 1–2 km layers in the upper crust, 3–5 km deeper.
* **Velocity steps:** every step is at least 0.05 km/s, so PyVelest's random-model generator stays
  within ~0.06 km/s of the starting model. As a consequence the flat 6.0 km/s section of kiratzi87 starts
  at 6.05–6.20.
* **Moho:** at 30 km (kiratzi87) / 40 km (ionio), with 8.0 km/s below. Your files put the first
  "mantle" row at 70 km, which a TauP reader turns into a 6→8 km/s gradient; that cannot be expressed in
  VELEST, and this geometry does not resolve it anyway (§6).

![starting and final models](../results/ioannina/comparison/models_comparison.png)

---

## 4. Pipeline settings for this dataset

| Stage | What runs | Key settings |
|---|---|---|
| 1 | 5 runs × 50 simultaneous inversions from perturbed starting models; each run's average model is merged and seeds the next run | Events: box, depth > 0.1 km, ≥ 8 stations with usable picks, gap ≤ 180° (usable stations), hypoSA RMS ≤ 0.5 s → 739 candidates → 500 spread over a 6×6×6 grid. Damping: 1.0 (velocity), 0.1 (station delays). Reference station: JAN (highest score) |
| 2 | Single-event relocation of all 10,362 events with the stage-1 model and delays, then reselection | `RESELECTION`: erh ≤ 4 km, 0.5 ≤ erz ≤ 4.5 km, RMS ≤ 0.3 s → 1,375–1,555 candidates → 500 |
| 3 | 5 runs × 25 inversions on the reselected events | as stage 1 |
| 4 | Stability test: unperturbed reference inversion (`run_00`), then 5 inversions with hypocentres shifted N(0, 5 km) in x, y and z | Final positions are compared with `run_00` |
| 5 | Final single-event relocation of all events with the final model and delays → `outputs/stage_fn` | Station delays are used (`use_station_corrections: true`) |

Other settings, in all stages:

* VELEST coordinate origin: 39.65°N 20.60°E; VELEST counts East longitudes as negative, so `olon` = −20.60.
* Maximum epicentral distance: 150 km.
* `phase_mode: 3`, Vp/Vs 1.82.

---

## 5. Results

### 5.1 Velocity model

All 375 inversions of every run converged (5 × 50 in stage 1, 5 × 25 in stage 3). Final models (Vp in
km/s at each depth):

| depth (km) | kiratzi87 linear | ionio linear | kiratzi87 layer-top | ionio layer-top | spread |
|---|---|---|---|---|---|
| 0.5 | 3.51 | 3.50 | 3.52 | 3.51 | 0.02 |
| 1.5 | 3.51 | 3.50 | 4.23 | 4.30 | 0.80 |
| 2.5 | 3.51 | 3.92 | 4.23 | 4.30 | 0.79 |
| 3.5 | 5.51 | 5.46 | 4.23 | 5.66 | 1.43 |
| 5–7 | 5.51 | 5.46 | 5.72 | 5.66 | 0.26 |
| 9 | 6.12 | 5.92 | 5.72 | 5.66 | 0.46 |
| 11 | 6.12 | 6.11 | 6.25 | 5.96 | 0.29 |
| 13–15 | 6.33 | 6.11 | 6.25 | 5.96 | 0.37 |
| 18 | 6.33 | 6.52 | 6.25 | 6.41 | 0.27 |
| 22 | 6.33 | 6.95 | 6.91 | 6.41 | 0.62 |
| 27+ | 7.96 | 6.95 | 7.93 | 6.41 | 1.55 |
| **final RMS (stage 3)** | **0.125 s** | 0.130 s | 0.136 s | 0.128 s | |

What the data constrain:

* **The average velocity of the upper crust.** The models place their layer boundaries differently, but
  the **vertical travel time from 0 to 10 km is 2.03–2.08 s in all four, and from 0 to 15 km it is
  2.875–2.879 s**. In other words, 4.8–4.9 km/s averaged over the top 10 km and 5.2 km/s over the top
  15 km.
* **A slow near-surface layer, about 2–3 km thick.** Vertical travel time through the top 3 km is
  0.76–0.83 s in all four models. Its velocity sits at the 3.5 km/s floor of the top layer in every run,
  so the data would accept it even slower. That layer trades off with the station delays, so do not read
  its exact velocity or thickness. If the shallow structure matters, test a lower `--vtop` (e.g. 3.0).
* **About 5.5–5.7 km/s from ~3 to ~8 km, rising to 5.9–6.3 km/s by 9–13 km.**
* **Not resolved below ~20 km.** Both kiratzi87 runs moved their 25–30 km layer to ~7.9 km/s, which merged
  with the half-space at 24–25 km. The ionio runs keep 6.4–6.95 km/s down to their fixed 40 km Moho.
  First arrivals within 125 km are crustal phases, so no Moho depth can be inferred.

**Recommended model: kiratzi87, linear reading** (`results/ioannina/kiratzi87/model.mod`). It has the
lowest RMS and the simplest layering, and it is stable in the stability test. Use the spread of the four
models as the uncertainty.

```
 3.50  -3.00 km   (dummy layer, station elevations)
 3.51   0.00 km
 5.51   2.78 km
 6.12   8.35 km
 6.33  12.91 km
 7.96  24.32 km   (half-space; not resolved)
```

### 5.2 Station delays

![station delays](../results/ioannina/comparison/station_delays.png)

The reference station is JAN, whose delay is 0 by definition. The pattern is consistent across the
starting models:

* **Negative delays** (faster than the model) at IGT, KSTR and TETR: −0.2 to −0.6 s.
* **Positive delays** at LFKM (+0.3 to +0.6 s), NEST and the most distant stations, OTHN and KZN
  (+0.3 to +1.5 s).
* **Weakest at the distant stations:** their delays trade off with the unresolved deep structure, which
  is why OTHN and KZN vary most between runs.
* **Shift of the ionio layer-top run:** it carries a common shift of about −0.2 to −0.3 s at most stations.
  That is the same trade-off, with its slower 3.5–11 km layer.

The S delays VELEST applied are 1.82 × the P delays; see `stations.csv`.

### 5.3 Relocations

![seismicity before/after](../results/ioannina/comparison/seismicity_comparison.png)

For the recommended model, inside the target box:

* **Fit:** relocated RMS median 0.13 s, against a hypoSA median of 0.40 s.
* **Free depths:** 6,977 of 9,526 events have a free depth (`depth_reliable`), with depths
  3.2 / 7.6 / 14.4 km at the 10th / 50th / 90th percentile.
* **Shape of the sequence:** it collapses into a narrow N–S zone at 20.62–20.65°E (39.55–39.72°N), mostly
  2–13 km deep, plus a second cluster around 20.55°E / 39.66–39.72°N.

The other three starting models give the same picture (median depth 7.1–9.2 km).

`events.csv` flags depths that are not a free solution:

* `at_layer_boundary`: within 0.1 km of a model layer top. VELEST hypocentres pile up there, about three
  times more often than by chance.
* `depth_error_near_zero`: erz < 0.5 km, meaning the depth was held.
* `depth_not_adjusted`: one iteration only, or still at the 10 km start depth.
* `at_surface`: depth below 0.05 km.

`depth_reliable` is true only when none of these applies. VELEST's formal errors are large (median
erh ≈ 4.0 km, erz ≈ 4.1 km), so use them to rank events, not as uncertainties.

### 5.4 Stability test (stage 4)

The test shifts the 500 events by N(0, 5 km) in each direction and re-inverts. The relocated events
return to the reference positions:

| start model | mean offset lat / lon / depth (km) | σ lat / lon / depth (km) |
|---|---|---|
| kiratzi87 linear | 0.06 / 0.04 / −0.06 | 0.27 / 0.14 / 1.67 |
| ionio linear | 0.13 / 0.22 / −0.33 | 0.24 / 0.13 / 1.48 |
| kiratzi87 layer-top | −0.13 / −0.06 / 0.09 | 0.21 / 0.12 / 1.25 |
| ionio layer-top | 0.22 / 0.23 / −0.21 | 0.26 / 0.11 / 2.58 |

Epicentres are stable to a few hundred metres. Depths are reproducible to about ±1.5–2.5 km, which is
the realistic depth uncertainty for this network.

---

## 6. Caveats

* **Depth resolution is limited.** The nearest station is ~21 km from the sequence. Take ±2 km from the
  stability test as the depth uncertainty, and expect trade-offs between depth, origin time and the
  shallow velocities.
* **Lateral structure is averaged.** The rays cross different tectonic units (Ionian zone and the Pindus
  front) on their way to stations up to 125 km away. The station delays absorb this lateral structure, so
  they belong with the model and should not be used without it.
* **The deep crust and the Moho are not constrained** (see §5.1).
* **VELEST relocations are a by-product.** For the final catalogue, use the minimum 1-D model with the
  station delays in your locator, as in the LOC-FLOW chain: HypoInverse or NonLinLoc, then hypoDD or
  GrowClust. At minimum, filter `events.csv` on `depth_reliable`.
* **Runs need separate directories.** Both configs write to `outputs/` and `figures/`, so use a separate
  copy of the repository for each config.

---

## 7. What changed in the code

**Input and configuration**

* `tools/prepare_real_inputs.py` and `tools/nd_to_velest.py` (`--interp linear|step`, `--vtop`,
  `--min-increment`, `--smooth`).
* Event box (`EVENT_SELECTION`) separate from the station box (`STUDY_AREA`). Depth limits and a
  minimum depth error for selection.
* `RESELECTION` overrides for stages 2 and 5, and `PICK_WEIGHTS`.
* Start depth for unconstrained catalogue depths.
* A fixed reference station (`REFERENCE_STATION.code`), `n_jobs`, `velest_exec`, and
  `fix_first_layer_bottom`.
* The YAML keys now match the VELEST config fields. Five keys used to be ignored silently; unknown keys
  are now errors.

**Correctness fixes.** Three review rounds, each finding checked independently, found these. Those
marked † were already in the original code and would affect other datasets too.

* † **Ensemble averaging:** the run average used one global iteration or state type, which silently
  dropped most inversions. Some stage-1 "averages" came from 2 of 50 models. Each inversion now
  contributes its own final state.
* † **Layer merging:** merges used chained 50/50 means; they now use the thickness-weighted slowness,
  which preserves travel time.
* † **Stage-1 damping:** 0.1 / 0.01 made up to 86 % of inversions stop with "error in (damped) matrix G",
  and the average came from the survivors. It is now 1.0 / 0.1, and every inversion converges.
* † **Silent failures:** a failed VELEST run could pass silently and leave results from an earlier run in
  the average. Old outputs are now removed before each run, failures are reported, and only this run's
  models are read.
* † **Time parsing:** VELEST times with minute or second = 60 crashed the parser.
* † **Station delays in relocation:** stages 2 and 5 relocated without the station delays.
* † **P+S models:** reading and writing a model with an S part crashed.
* † **Horizontal error filter:** `erh` was never read from Nordic input, so the filter did nothing.
* **Stability-test reference:** stage 4 compared against stage-2 positions computed with a different
  model; it now has an unperturbed reference run. Perturbed depths are no longer clipped to 0 km, which
  had inflated the depth scatter.
* **Pinned depths:** the start-depth fix for hypoSA's 0 km depths.
* **Usable picks only:** station counts and gaps now count only stations with usable picks.
* **S cards:** S continuation cards in VELEST output were parsed with distance and azimuth 0.
* **Reselection crash:** reselection crashed when VELEST could not locate an event.
* **`--stages`:** a partial run no longer deletes earlier stages.
* **Product honesty:** the final products report the S delays and Vs that were actually used, and
  `events.csv` carries the depth-quality flags.

**Practical changes**

* VELEST runs inside each model directory, with relative paths (it reads at most 80 characters per path),
  so runs can go in parallel.
* `pyvelest.py --config/--stages`.
* `src/exc/build_velest.sh` builds VELEST with gfortran.

---

## 8. Open questions for you

1. **`.nd` intent:** were `kiratzi87.nd` / `ionio.nd` meant as layer tops or as TauP gradients? If
   layer tops, the travel-time tables REAL built from them were gradient models.
2. **Mantle row:** what do you intend for the Moho and the "mantle" row at 70 km?
3. **Magnitudes:** I used the magnitudes in `ioannina_final.txt` rather than REAL's. Is that right?
4. **Unused stations:** do KBN, LSK and SRN have data? Adding picks from them would help depth
   resolution: KBN and LSK are north of the sequence, SRN to its north-west.
