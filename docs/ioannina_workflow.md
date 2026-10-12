# Using PyVelest with the Ioannina catalog: review and proposed workflow

This note covers three things:

1. How the inputs in `inputs/My inputs/` compare with what PyVelest expects.
2. Problems in the code that affect a run on this dataset.
3. A step-by-step workflow, and the information still needed.

To check the claims below, I ran parts of the pipeline on `ioannina_final.txt` in a scratch copy of the repo:
`load_catalog`, a reduced `stage_01` (1 run × 2–3 models) and `stage_02` with all 10,362 events.
**The station coordinates in that test were placeholders.** It shows that the file formats go through
VELEST without errors. The velocity models and locations it produced mean nothing.

---

## 1. Your inputs at a glance

| Input | What the code expects | What you have | Status |
|---|---|---|---|
| Catalog | ObsPy-readable file (the example uses Nordic `sample.out`) with picks, location quality (erh/erz/rms/gap/nst) and Nordic pick weights 0–4 | `ioannina_final.txt`: HypoDD phase format, 10,362 events, 2026-03-01 → 2026-08-31, 130,843 picks (69,153 P / 61,690 S), 13 stations | **Readable as is** (ObsPy reads the HypoDD phase format). Quality info and weights are lost; see §2 |
| Stations | `inputs/stations.csv` with `code,lat,lon,elv` (elv in m, ≥ 0) | **missing** | **Blocker** |
| Velocity model | VELEST `.mod` text file (`inputs/vmodel.dat`) | `kiratzi87.nd` (TauP `.nd` format, whole-Earth) | Needs conversion; see §3 |

Facts about the catalog that shape the workflow:

* **The header columns `eh`, `ez` and `rms` are all 0.0**, so the catalog has no location-quality information.
* **4,999 events (48 %) are at exactly 0.0 km depth.** All depths lie on a 2 km grid (0, 2, 4 … 30 km).
  That looks like an associator's grid search (e.g. REAL) rather than a final location, and depth 0 probably
  means "unconstrained / top of grid". These hypocentres are poor starting points for a simultaneous inversion.
* **Pick weights range from 0.30 to 0.987.** That looks like ML picker probabilities (0.3 is PhaseNet's
  default threshold), not HypoDD weights.
* **Phases per event:** 8–24 (median 12, about 6 stations).
  Events with ≥ 8 P picks: 3,003. Events with ≥ 10 P picks: 1,167. Events with ≥ 12 stations: only 182.
* **Picks per station:** KSTR 19,822 · IGT 19,292 · PRMD 19,246 · TETR 18,936 · JAN 18,228 · PENT 10,858 ·
  KPRO 7,350 · KEK 7,296 · AMPL 3,518 · NEST 3,324 · LFKM 2,314 · KZN 540 · OTHN 119.
  All station codes are ≤ 4 characters, which VELEST requires.
* **Vp/Vs from a Wadati fit** on 29k P–S pairs (P prob ≥ 0.7, S prob ≥ 0.55) is **≈ 1.82**.
  That is consistent with the 1.80 in `kiratzi87.nd`.

---

## 2. Code review: problems that affect your run

File references are relative to the repo root.

### Blocking or high impact

1. **No station coordinates.** Everything depends on them: distances, azimuths, the reference station,
   the station file and VELEST itself.

2. **Crash on VELEST origin times with minute = 60.** `add_origin` (`src/io/catalog.py:410`) parses
   VELEST's single-event output with `strptime`. VELEST sometimes prints times such as `15:60:00.262`.
   Stage 2 on your catalog **crashed** on this. The ~120 example events never triggered it.
   `cnv2csv` (`catalog.py:522`) has a similar second = 60 workaround that clips the time to 59.99 s.
   *Fix:* parse hour, minute and second separately and add them as a `timedelta`.
   I tested this in the scratch copy.

3. **S picks are thrown away.** `phase_mode: 1` (`configs/stage_0X/velest.yaml:108`) tells VELEST
   "P only", and VELEST then skips every S pick (`velest.f:1966`). That is 47 % of your data.
   Worse, the stage-2 and stage-5 single-event outputs list only the phases VELEST used.
   The reselected `select.cnv` therefore contains only P, and stages 3–4 cannot recover S even if you
   change the setting later.
   *Fix:* use `phase_mode: 3` (S from a fixed Vp/Vs; S constrains depth) with `vp_vs_ratio: 1.80` in
   **all five** stage configs. I tested that VELEST runs this way with your data.

4. **Pick probabilities are discarded.** The weight is read only from Nordic `nordic_pick_weight`
   (`catalog.py:70-74` and `catalog.py:355-362`). ObsPy stores your probabilities in `arrival.time_weight`,
   so every pick becomes class 0 (full weight). A proposed mapping:

   | prob. | VELEST class | P picks | S picks |
   |---|---|---|---|
   | ≥ 0.85 | 0 | 43 % | 7 % |
   | 0.70–0.85 | 1 | 22 % | 24 % |
   | 0.55–0.70 | 2 | 15 % | 30 % |
   | 0.40–0.55 | 3 | 13 % | 24 % |
   | < 0.40 | 4 (weight 0) | 7 % | 15 % |

5. **The event-selection filters have no effect on this catalog.** `select` (`catalog.py:151-164`)
   keeps every event where `erh/erz/rms/gap/nst` is missing (`| df[column].isna()`). For your catalog
   all five columns are empty. The within-cell ranking (`select3D.py:60-63`) sorts by those same empty
   columns, which leaves magnitude as the only criterion.
   In the test, 23 % of the 500 events picked for stage 1 had depth 0.
   Also, with only 13 stations, `min_used_stations: 12` would leave 182 events.

6. **Results from the example run can leak into yours.**
   * The repo has 5,185 output and figure files from the example dataset committed under `outputs/` and
     `figures/`, and `clean_workspace: false` is set (`configs/user.yaml:74`).
   * `finalized_inverted_models` averages **every** `velocity_summary.csv` it finds under a run directory
     (`src/io/model.py:614`).
   * `run_velest` silently skips a model whose VELEST run wrote to stderr (`src/core/velest.py:35`).

   So if you use fewer than 50 models, or if a run fails, the old example results are averaged into your
   model without any warning. **Empty `outputs/` and `figures/` before the first run**, or set
   `clean_workspace: true`.

7. **The origin longitude is still the example's, and VELEST's sign convention is easy to miss.**
   VELEST's `olon` is **positive West, negative East** (`velest.f:1320-1324`). All five stage configs have
   `origin_longitude: -46.75` (line 23), which is 46.75°E from the example. For your area use about
   `origin_latitude: 39.65`, `origin_longitude: -20.75`. Adjust these to the network centre once the
   station coordinates are known.

8. **Input paths don't match the repo layout and are partly hard-coded.**
   * `configs/user.yaml:2` still points to `inputs/sample.out`, which no longer exists since the inputs
     moved into subfolders. The repo as committed doesn't run.
   * `prepare_station` hard-codes `inputs/stations.csv` (`src/io/station.py:560`), and `prepare_model`
     hard-codes `inputs/vmodel.dat` (`src/io/model.py:508`). Both ignore `INPUTS`.
   * `catalog.csv` is always written to `inputs/` (`catalog.py:96`).

   For now, put the three input files directly in `inputs/` under those names.

### Medium

9. **Stages 2 and 5 relocate without station corrections.** Those configs set `use_station_corrections: false`
   (`configs/stage_02|05/velest.yaml:159`), and VELEST then zeroes the input delays (`velest.f:1636-1639`).
   The final catalog is therefore relocated with the minimum 1-D model but **without** the station terms
   inverted together with it. I recommend `true`.

10. **Some YAML settings are silently ignored.** Several keys don't match the pydantic field names,
    and pydantic ignores unknown keys:

    | YAML key | Model field |
    |---|---|
    | `resolution_matrix` | `compute_resolution` |
    | `topography` | `use_topography` |
    | `shot_corrections` | `use_shot_corrections` |
    | `write_summary_cards` | `write_summary_file` |
    | `ale` | `averaging_kernel` |

    For example, stages 2 and 5 ask for a resolution matrix and summary cards, but the generated
    `velest.cmn` has 0 for both.

11. **The reference-station ranking expects PG/SG phase names.** `calculate_reference_station_scores`
    (`station.py:286-294`) only counts `PG`/`SG`. Your phases are `P`/`S`, so the residual and phase-count
    parts score zero. Your catalog has no residuals anyway. The `REFERENCE_STATION` weights in `user.yaml`
    are never passed in (`station.py:484`).
    In effect, the reference station (the highest ICC, `velest.f:1652-1654`) is chosen by event coverage
    and distance only. It is safer to set it by hand: a central station with many picks. JAN, KSTR, IGT,
    PRMD and TETR are the best-recorded; which is most central depends on the coordinates.

12. **One box filters both events and stations.** `STUDY_AREA` selects events and also drops stations
    outside it (`station.py:21`). If you shrink the box to focus on the seismicity, peripheral stations
    such as OTHN, KZN or LFKM could disappear and take their picks with them. Keep the box wide enough for
    all stations and limit events with `max_gap_deg` instead.

13. **The random-model generator smooths sharp velocity jumps.** `generate_random_models`
    (`src/model/Models.py:363-453`) works as follows:
    * It perturbs the velocity *increments* and rejects any increment below 0.01 km/s.
    * It rescales the increments so the top and bottom velocities stay fixed.

    With a step-like model (big jumps next to near-zero increments), the ensemble mean drifts away from the
    starting model. With a Kiratzi-faithful step model I measured up to **0.28 km/s** of bias; with a
    smoothed version, ≤ 0.07 km/s. Two consequences:
    * Start from a smooth model (see §3).
    * The top of the half-space (the Moho depth) and its velocity are **never** perturbed. The Moho depth
      has to be tested with separate starting models.

14. **Unclear crash when nothing passes the selection.** If no event passes, `select3D.py:97-99` fails with
    `KeyError: ['ilat','ilon','idep']`. If you see that, the thresholds are too strict.

### Minor / for awareness

* **The bundled binary won't run here.** `src/exc/velest` was built with the Intel compiler and needs
  `libimf.so`. Unless you have the Intel runtime, rebuild it with gfortran. This compiled cleanly and ran
  your data:
  `cd src/exc && gfortran -O2 -std=legacy -ffixed-line-length-none -fno-automatic velest.f -o velest`
* **VELEST array limits** (`src/exc/vel_com.f:13-19`): 4,500 events in simultaneous mode (500 is fine),
  200 stations, 100 layers. Single-event mode (stages 2 and 5) has no event limit.
* **Model limits:** the model validators allow Vp ≤ 8.5 km/s and layer tops ≤ 60 km (`Models.py:157-158`).
  That is fine for a crustal model.
* **Workspace directories:** `prepare_workspace` sizes the `stage_03` folders with `STAGE_01.n_runs` and
  `SYNTHETIC_MODELS.n_models` (`src/io/workspace.py:21-29`). Keep `STAGE_03` values ≤ those.
* **Selection semantics change between stages:** after stage 2 the `nst` column is VELEST's *number of
  observations* (P + S), not a station count. The same `min_used_stations` threshold means different
  things in stage 1 and stage 2.
* **Run time:** `load_catalog` takes about 1 min on your file. Each VELEST simultaneous inversion
  (500 events) took about 15 s here, so the default stage 1 (5 runs × 50 models) is roughly an hour.

---

## 3. Velocity model: from `kiratzi87.nd` to VELEST

`kiratzi87.nd` is a TauP model of the whole Earth. VELEST only needs the crust and the uppermost mantle.
Read as layer tops, the crust is:

| top (km) | Vp | Vs | Vp/Vs |
|---|---|---|---|
| 0 | 4.20 | 2.33 | 1.80 |
| 2 | 5.40 | 3.00 | 1.80 |
| 10 | 6.00 | 3.33 | 1.80 |
| 30 | 6.00 | 3.33 | 1.80 |
| *mantle* 70 | 8.00 | 4.49 | — |

**The mantle part needs clarifying.** In TauP `.nd` semantics, the lines `30 km → 6.00` then `70 km → 8.00`
mean a *linear gradient* from 6.0 to 8.0 km/s between 30 and 70 km, with no Moho discontinuity. Most likely
the intent was a Moho at about 30 km with Pn ≈ 8.0 km/s. This matters: Pn arrives first beyond roughly
100–150 km, which is about your network's aperture.

**Proposed starting model** (P only, to be used with `phase_mode: 3` and `vp_vs_ratio: 1.80`).
It keeps the Kiratzi values at 2, 10 and 30 km but smooths between them, so the random-model generator
doesn't bias it (item 13). The top layer starts above the highest station, as VELEST requires.

```
KIRATZI87_GRADIENT
 14        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)
 4.20     -3.00    50.00           P-VELOCITY MODEL
 4.80      1.00    1.000
 5.20      2.00    1.000
 5.40      4.00    1.000
 5.60      6.00    1.000
 5.80      8.00    1.000
 6.00     10.00    1.000
 6.10     13.00    1.000
 6.20     16.00    1.000
 6.30     19.00    1.000
 6.40     22.00    1.000
 6.50     26.00    1.000
 6.60     28.00    1.000
 8.00     30.00    10.00
```

The minimum-1-D-model practice (Kissling et al., 1994; Kissling, 1995) is to **run several starting models**
and check that they converge. Good candidates:

* the model above;
* the same model with the Moho at 35 km and at 40 km (the generator never moves it; Moho depth varies
  across Epirus and the Pindus);
* a faster and a slower variant (± 0.2 km/s);
* if you have it, a published minimum 1-D model for the area. Haslinger et al. (1999, *Tectonophysics*)
  derived one for the Gulf of Arta / Ionian region with VELEST.

---

## 4. Proposed workflow

### Step 0: Environment
* Python ≥ 3.10 with `obspy pandas pydantic>=2 pyyaml tqdm matplotlib ultraplot`.
* Rebuild VELEST with gfortran (command above).
* Empty `outputs/` and `figures/`, or set `clean_workspace: true` for the first run.

### Step 1: Inputs (flat in `inputs/`)
* `inputs/stations.csv`: `code,lat,lon,elv` for the 13 stations. **You provide this.**
* `inputs/ioannina_final.txt`: the catalog unchanged; the code changes in Step 2 handle weights and quality.
* `inputs/vmodel.dat`: the starting model from §3.

### Step 2: Small code changes (I can make these once you confirm)
1. Map `arrival.time_weight` (pick probability) to VELEST classes 0–4 (table in §2, item 4).
2. Fix the minute/second = 60 parsing (§2, item 2).
3. Add a **stage 0**: relocate *all* events once in single-event mode with the starting model.
   This reuses the stage-2 machinery (`reselect_best`). Stage 1 then selects from VELEST's real
   erh/erz/rms/gap/nst instead of the empty catalog values.
   It also replaces the depth-0 / 2 km-grid hypocentres with real locations before the simultaneous inversion.
   *(A simpler fallback: compute nst and the azimuthal gap from the picks and the station coordinates,
   and leave depth-0 events out of the stage-1 selection.)*
4. Accept `P`/`S` in the reference-station scoring and pass the `REFERENCE_STATION` weights,
   or allow a fixed reference station.
5. Fix the YAML key names (§2, item 10); set `use_station_corrections: true` in stages 2 and 5.
6. Use the `INPUTS` paths everywhere; raise clear errors on an empty selection or a failed VELEST run.

### Step 3: Configuration
`configs/user.yaml` (starting values, to tune after stage 0):
```yaml
INPUTS:
  catalog: inputs/ioannina_final.txt
  stations: inputs/stations.csv
  vmodel: inputs/vmodel.dat
STUDY_AREA:              # must contain ALL stations (the box filters stations too)
  min_lat: 38.60
  max_lat: 40.70
  min_lon: 19.40
  max_lon: 22.20
  lat_hemisphere: N
  lon_hemisphere: E
EVENT_SELECTION:
  max_events: 500
  cell_count_lat: 6
  cell_count_lon: 6
  cell_count_depth: 5
  max_horizontal_error_km: 2.0
  max_depth_error_km: 3.0
  max_rms_sec: 0.3
  max_gap_deg: 180
  min_used_stations: 8   # after a VELEST relocation this is the number of P+S observations
FIGURE_SETTINGS:
  min_sei_lat: 38.60
  max_sei_lat: 40.70
  min_sei_lon: 19.40
  max_sei_lon: 22.20
  max_sei_depth: 30.0
  max_vel_depth: 40.0
GLOBAL_SETTINGS:
  clean_workspace: true
```
All five `configs/stage_0X/velest.yaml`:
```yaml
COORDINATE_SYSTEM:
  origin_latitude: 39.65
  origin_longitude: -20.75   # VELEST: East longitudes are NEGATIVE
VELOCITY:
  phase_mode: 3
  vp_vs_ratio: 1.80
RELOCATION:
  maximum_distance: 150.0    # check against the network geometry; limits Pn and 3-D effects
```
In stages 2 and 5 only: `use_station_corrections: true`.

### Step 4: Stage 0, initial relocation (new)
Locate all 10,362 events with the starting model.
* Look at the distributions of RMS, GAP, erh, erz and the number of observations, and set the
  `EVENT_SELECTION` thresholds from them.
* Check how many of the former depth-0 events get a constrained depth.

### Step 5: Stage 1, minimum 1-D model (ensemble)
`n_runs` × `n_models` VELEST inversions on the about 500 best, spatially balanced events, starting from
perturbed versions of the model.

Check:
* `figures/models_st01_run_XX.png` and `statistics_st01_run_XX.png`: the ensemble should converge and the
  RMS should drop.
* The station delays in `run_XX/summary/stations.sta`: they should be small and geologically sensible.
* The layers merged by `merging_threshold_km_s`.

### Step 6: Stage 2, relocate and reselect
Relocate all events with the stage-1 model and station delays, then reselect the best events using
VELEST's errors.

### Step 7: Stage 3, refinement
Repeat the ensemble inversion on the reselected events, starting from the stage-2 model.

### Step 8: Stage 4, stability test
Shift the hypocentres randomly (`dx/dy/dz_km`, 10 km by default) and invert again. The events should move
back and the model shouldn't change; see `figures/dislocation_XX.png`. Also run Steps 5–8 for the other
starting models (§3): a robust minimum 1-D model should not depend on which one you start from.

### Step 9: Stage 5, final products
Relocate everything with the final model and station delays. The outputs go to `outputs/stage_fn/`:
`model.mod`, `model.csv`, `stations.sta`, `stations.csv`, `relocated.csv` and `relocated.out`.

---

## 5. Information needed from you

1. **The station list** (required): network, code, latitude, longitude and elevation in metres for
   KSTR, IGT, PRMD, TETR, JAN, PENT, KPRO, KEK, AMPL, NEST, LFKM, KZN and OTHN.
   StationXML or a plain table is fine.
2. **How the catalog was made:** which picker (are the weights pick probabilities?), which
   associator/locator (REAL, GaMMA, HypoInverse, NonLinLoc …), and what the 0.0 km depths mean.
   If you have locator output with errors, RMS and GAP, it can replace stage 0.
3. **The velocity model:**
   * Is the Kiratzi model meant to have a Moho at 30 km with 8.0 km/s, or the 30→70 km gradient?
   * Its reference.
   * Any other regional models, or Moho depths (e.g. from receiver functions), you'd like to test.
4. **The target area:** the whole 38.6–40.7°N / 19.4–22.1°E cloud, or a specific sequence?
   This sets the study box and the coordinate origin.
5. **Known station issues** during March–August 2026: timing problems, station moves or replacements.
