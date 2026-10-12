from __future__ import annotations

import re
from pathlib import Path

from obspy import UTCDateTime as utc
from obspy.core.event import (
    Catalog,
    Event,
    Origin,
    Magnitude,
    Pick,
    Arrival,
    WaveformStreamID,
)

from obspy import read_events
from pandas import DataFrame, read_csv, to_numeric
from src.io.select3D import select_events_3d
from datetime import timezone
from obspy.geodetics.base import gps2dist_azimuth as gps
from obspy.geodetics.base import kilometer2degrees as k2d
from shutil import copy


from datetime import datetime as dt, timedelta
from math import sqrt, cos, radians
from numpy import nan, random
from src.io.station import read_station

import warnings

warnings.filterwarnings("ignore")


CATALOG_SUMMARY = Path("inputs/catalog.csv")


def pick_weight(pick, arrival, pick_weights=None, use_pick_time_errors=False):
    """VELEST weight class (0 = best ... 4 = unused) of one pick.

    Priority: a Nordic operator weight; then a pick probability in
    ``arrival.time_weight`` (0-1, e.g. PhaseNet via a HypoDD phase file),
    mapped with ``pick_weights``; then, optionally, the pick time
    uncertainty; otherwise 0.
    """
    extra = getattr(pick, "extra", None) or {}
    if "nordic_pick_weight" in extra:
        try:
            return int(extra["nordic_pick_weight"]["value"])
        except Exception:
            pass

    probability = getattr(arrival, "time_weight", None)
    if (
        pick_weights is not None
        and pick_weights.enabled
        and probability is not None
        and 0.0 <= probability <= 1.0
    ):
        return pick_weights.weight_class(probability)

    if use_pick_time_errors and pick.time_errors:
        err = pick.time_errors.uncertainty
        if err is not None:
            if err <= 0.05:
                return 0
            elif err <= 0.10:
                return 1
            elif err <= 0.20:
                return 2
            elif err <= 0.50:
                return 3
            return 4
    return 0


def azimuthal_gap(azimuths):
    """Largest azimuthal gap (degrees) between the given azimuths."""
    azimuths = sorted({round(a % 360.0, 3) for a in azimuths})
    if len(azimuths) < 2:
        return 360.0
    gaps = [b - a for a, b in zip(azimuths, azimuths[1:])]
    gaps.append(azimuths[0] + 360.0 - azimuths[-1])
    return max(gaps)


def station_coordinates(config):
    """{code: (lat, lon)} of the stations inside STUDY_AREA."""
    stations_df = read_station(config, config.INPUTS.stations)
    return {
        str(row.code): (row.lat, row.lon)
        for row in stations_df.itertuples(index=False)
    }


def get_picks_info(event, eid, station_info, pick_weights=None):
    """Extract VELEST-ready pick rows from one ObsPy event.

    ``station_info`` maps station code -> (lat, lon). Distances are in
    kilometers, azimuths are from the event to the station (as in VELEST
    output) and travel times are in seconds. Picks whose phase or station
    cannot be resolved are ignored.
    """
    origin = event.preferred_origin()
    if origin is None:
        if not event.origins:
            raise ValueError(f"Event {eid} has no origin")
        origin = event.origins[0]

    phases = []
    for arrival in origin.arrivals:
        phase = (arrival.phase or "").upper()
        if not phase or phase[0] not in "PS" or arrival.pick_id is None:
            continue

        pick = arrival.pick_id.get_referred_object()
        if pick is None or pick.waveform_id is None:
            continue
        sta = pick.waveform_id.station_code
        if sta not in station_info:
            continue

        distance_m, azimuth, _ = gps(
            origin.latitude, origin.longitude, *station_info[sta]
        )
        dist = distance_m * 1e-3
        weight = pick_weight(pick, arrival, pick_weights)
        phases.append(
            {
                "eid": eid,
                "sta": sta,
                "dis": dist,
                "azm": azimuth,
                "pha": phase,
                "wet": weight,
                "ptt": pick.time - origin.time,
                "res": arrival.time_residual,
            }
        )
    return phases


def _horizontal_error_km(origin):
    """Horizontal location error in km, from whatever the reader filled."""
    ou = origin.origin_uncertainty
    if ou is not None:
        for value in (
            ou.horizontal_uncertainty,
            ou.max_horizontal_uncertainty,
        ):
            if value:
                return value * 1e-3  # QuakeML stores metres
    lat_err = getattr(origin.latitude_errors, "uncertainty", None)
    lon_err = getattr(origin.longitude_errors, "uncertainty", None)
    if lat_err and lon_err:
        dy = lat_err * 111.19
        dx = lon_err * 111.19 * cos(radians(origin.latitude))
        return sqrt(dx**2 + dy**2)
    return None


def summarize_catalog(config):

    # inputs
    cat_inpfile = config.INPUTS.catalog

    # outputs
    cat_outfile = CATALOG_SUMMARY

    # processing
    cat = read_events(cat_inpfile)
    station_info = station_coordinates(config)
    data = []
    for eid, event in enumerate(cat, 1):
        po = event.preferred_origin() or event.origins[0]
        pm = event.preferred_magnitude() or None
        quality = po.quality
        phase_info = get_picks_info(
            event, eid, station_info, config.PICK_WEIGHTS
        )
        # Stations with at least one usable pick (class 4 is dropped from
        # simultaneous inversions) define nst and the azimuthal gap.
        picked = {}
        for phase in phase_info:
            if int(phase["wet"]) < 4:
                picked[phase["sta"]] = phase["azm"]

        ort = po.time
        lat = po.latitude
        lon = po.longitude
        dep = po.depth * 1e-3
        mag = pm.mag if pm else None
        # Use the catalog's own values when present, otherwise compute
        # them from the picks at stations inside STUDY_AREA.
        gap = getattr(quality, "azimuthal_gap", None) or (
            azimuthal_gap(picked.values()) if picked else None
        )
        rms = getattr(quality, "standard_error", None) or None
        nst = getattr(quality, "used_station_count", None) or (
            len(picked) or None
        )
        erh = _horizontal_error_km(po)
        erz = getattr(po.depth_errors, "uncertainty", None)
        erz = erz * 1e-3 if erz else None
        event_info = {
            "eid": eid,
            "ort": ort,
            "lat": lat,
            "lon": lon,
            "dep": dep,
            "mag": mag,
            "gap": gap,
            "rms": rms,
            "nst": nst,
            "erh": erh,
            "erz": erz,
        }
        phase_info.append(event_info)
        data.extend(phase_info)
    cat_df = DataFrame(data)
    cat_df.to_csv(cat_outfile, index=False, float_format="%.3f")
    return cat


def select(
    config, catalog, catalog_sum_file, select_cat_outfile, select_csv_outfile
):

    # process
    df = read_csv(catalog_sum_file)
    # Events in ``catalog`` are in ascending eid order, including events
    # VELEST could not locate (no origin time; make_catalog keeps a
    # placeholder for them), so take the eids before dropping those.
    all_eids = sorted(df["eid"].dropna().astype(int).unique())
    df.dropna(subset=["ort"], inplace=True)
    if len(all_eids) != len(catalog):
        raise ValueError(
            f"{catalog_sum_file} lists {len(all_eids)} events but the "
            f"catalog has {len(catalog)}"
        )
    event_by_eid = dict(zip(all_eids, catalog))

    es, sa = config.EVENT_SELECTION, config.STUDY_AREA

    def bound(value, fallback):
        return fallback if value is None else value

    mask = df["lat"].between(
        bound(es.min_lat, sa.min_lat), bound(es.max_lat, sa.max_lat)
    ) & df["lon"].between(
        bound(es.min_lon, sa.min_lon), bound(es.max_lon, sa.max_lon)
    )
    if es.min_depth_km is not None:
        mask &= df["dep"] >= es.min_depth_km
    if es.max_depth_km is not None:
        mask &= df["dep"] <= es.max_depth_km

    limits = (
        ("erh", config.EVENT_SELECTION.max_horizontal_error_km, "le"),
        ("erz", config.EVENT_SELECTION.max_depth_error_km, "le"),
        ("erz", config.EVENT_SELECTION.min_depth_error_km, "ge"),
        ("rms", config.EVENT_SELECTION.max_rms_sec, "le"),
        ("gap", config.EVENT_SELECTION.max_gap_deg, "le"),
        ("nst", config.EVENT_SELECTION.min_used_stations, "ge"),
    )
    for column, limit, operator in limits:
        if limit is None:
            continue
        comparison = (
            df[column] <= limit if operator == "le" else df[column] >= limit
        )
        mask &= comparison | df[column].isna()

    filtered = df[mask].copy()
    if filtered.empty:
        raise ValueError(
            f"No event out of {len(df)} passes the EVENT_SELECTION "
            f"criteria ({catalog_sum_file}); relax the thresholds."
        )
    selected = select_events_3d(filtered, config)
    print(
        f"+++ {len(filtered)} of {len(df)} events pass the selection "
        f"criteria; {len(selected)} selected"
    )
    selected.to_csv(select_csv_outfile, index=False, float_format="%.3f")
    cat_sel = Catalog()
    for eid in selected.eid:
        cat_sel.events.append(event_by_eid[int(eid)])
    cat_sel.write(
        select_cat_outfile,
        format="NORDIC",
        high_accuracy=False,
        nordic_format="OLD",
    )
    return cat_sel


def make_origin(origin_data):
    origin = Origin(
        time=origin_data["time"],
        latitude=origin_data["latitude"],
        longitude=origin_data["longitude"],
        depth=(
            origin_data["depth_km"] * 1000.0
            if origin_data["depth_km"]
            else None
        ),
    )
    return origin


def make_magnitude(origin_data):
    magnitude = Magnitude(
        mag=origin_data["magnitude"],
        magnitude_type="ML",
    )
    return magnitude


def make_catalog(events_df):
    catalog = Catalog()
    events_df.replace({nan: None}, inplace=True)
    for eid, event_df in events_df.groupby(["eid"]):

        event_info = event_df.iloc[-1]
        event = Event()

        if not event_info.ort:
            origin_data = {
                "time": utc(1900, 1, 1),
                "latitude": 0,
                "longitude": 0,
                "depth_km": 0,
                "magnitude": 0,
            }
        else:
            origin_data = {
                "time": utc(event_info.ort),
                "latitude": event_info.lat,
                "longitude": event_info.lon,
                "depth_km": event_info.dep,
                "magnitude": event_info.mag,
            }

        origin = make_origin(origin_data)
        event.origins.append(origin)
        event.preferred_origin_id = origin.resource_id

        magnitude = make_magnitude(origin_data)
        event.magnitudes.append(magnitude)
        event.preferred_magnitude_id = magnitude.resource_id
        catalog.events.append(event)

        for r, row in event_df.iterrows():
            if row.ptt == None:
                continue
            pick = Pick(
                time=origin_data["time"] + row.ptt,
                waveform_id=WaveformStreamID(station_code=row.sta),
                phase_hint=row.pha,
                evaluation_mode="automatic",
            )
            pick.extra = {
                "nordic_pick_weight": {
                    "value": row.wet or 0,
                    "namespace": "velest",
                }
            }
            event.picks.append(pick)

            arrival = Arrival(
                pick_id=pick.resource_id,
                phase=row.pha,
                time_residual=row.res,
            )
            origin.arrivals.append(arrival)
    return catalog


def catalog_to_cnv(
    catalog,
    stations_inpfile,
    root_out,
    cnv_outfile,
    pick_weights=None,
    use_pick_time_errors=False,
    unconstrained_depth_below_km=None,
    start_depth_km=10.0,
):
    """
    Convert an ObsPy Catalog to VELEST CNV format.

    Parameters
    ----------
    catalog : obspy.core.event.Catalog
    csv_outfile : str or Path
    default_weight : int
        Weight assigned when no uncertainty exists.
    use_pick_time_errors : bool
        Convert pick.time_errors.uncertainty to VELEST weights.
    """

    # outputs
    root_out = Path(root_out)
    unused_st_file = root_out / "unused_st.csv"

    # processing
    stations_df = read_csv(stations_inpfile)
    station_codes = set(stations_df["code"].astype(str))
    unused_st = []

    with open(cnv_outfile, "w") as f:

        for event in catalog:

            if event.preferred_origin():
                origin = event.preferred_origin()
            else:
                origin = event.origins[0]

            if event.preferred_magnitude():
                magnitude = event.preferred_magnitude().mag
            elif event.magnitudes:
                magnitude = event.magnitudes[0].mag
            else:
                magnitude = 0.0
            magnitude = magnitude or 0.0

            ot = origin.time.datetime.replace(tzinfo=timezone.utc)

            lat = abs(origin.latitude)
            lon = abs(origin.longitude)

            ns = "N" if origin.latitude >= 0 else "S"
            ew = "E" if origin.longitude >= 0 else "W"

            depth = origin.depth / 1000.0 if origin.depth else 0.0
            # Unconstrained (e.g. pinned 0 km) depths: start VELEST deeper.
            if (
                unconstrained_depth_below_km is not None
                and depth < unconstrained_depth_below_km
            ):
                depth = start_depth_km

            # -------------------------
            # Event header
            # -------------------------

            header = (
                f"{ot:%y%m%d} "
                f"{ot:%H%M} "
                f"{ot.second + ot.microsecond/1e6:05.2f} "
                f"{lat:7.4f}{ns} "
                f"{lon:8.4f}{ew} "
                f"{depth:7.2f}  "
                f"{magnitude:5.2f}"
            )

            f.write(header + "\n")

            observations = []
            repetitions = []

            for arrival in origin.arrivals:
                if arrival.pick_id is None:
                    continue
                pick = arrival.pick_id.get_referred_object()
                if pick is None or pick.waveform_id is None:
                    continue
                sta = pick.waveform_id.station_code[:4].ljust(4)
                if sta.strip() not in station_codes:
                    unused_st.append(sta.strip())
                    continue
                if not arrival.phase or not arrival.phase.strip():
                    continue
                phase = arrival.phase.upper()[0]
                tt = pick.time - origin.time
                weight = pick_weight(
                    pick, arrival, pick_weights, use_pick_time_errors
                )
                obs = f"{sta}{phase}{weight:d}{tt:6.2f}"
                rep = f"{sta}{phase}"
                observations.append(obs) if rep not in repetitions else None
                repetitions.append(rep)

            # six observations per line
            for i in range(0, len(observations), 6):
                f.write("".join(observations[i : i + 6]) + "\n")

            f.write("\n")

        f.write("9999\n")

    unused_st = sorted(set(unused_st))
    unused_st_df = DataFrame({"code": unused_st})
    with open(unused_st_file, "w") as f:
        f.write(
            "# Stations in the catalog but outside the specified region are listed below:\n"
        )
        unused_st_df.to_csv(f, index=False)


def get_eid(line):
    args = line.split()
    return int(args[-3])


def _two_digit_year(yy):
    """Same pivot as strptime's %y: 69-99 -> 19xx, 00-68 -> 20xx."""
    return 1900 + yy if yy >= 69 else 2000 + yy


def add_origin(line, header):
    # VELEST may print minute or second values of 60 (e.g. 15:60:00.262),
    # so build the time from its parts instead of using strptime.
    date = re.sub(r"\s", "0", line[1:7])
    hh, mm, ss = line[8:20].split(":")
    ort = dt(
        _two_digit_year(int(date[:2])), int(date[2:4]), int(date[4:6])
    ) + timedelta(hours=int(hh), minutes=int(mm), seconds=float(ss))
    lat = float(line[21:28])
    lon = float(line[30:37])
    dep = float(line[39:46])
    mag = float(line[47:51])
    nobs = float(line[52:55])
    dmn = float(line[56:59])
    gap = float(line[60:63])
    rms = float(line[64:68])
    header.update(
        {
            "ort": ort,
            "lat": lat,
            "lon": lon,
            "dep": dep,
            "mag": mag,
            "nobs": nobs,
            "dmn": dmn,
            "gap": gap,
            "rms": rms,
        }
    )
    return header


def add_statistics(line, header):
    args = line.split()
    erx = float(args[0])
    ery = float(args[1])
    erz = float(args[2])
    erh = sqrt(erx**2 + ery**2)

    header.update({"erh": erh, "erz": erz})
    return header


def add_phase_info(line):
    sta = line[2:6].strip()
    dis = float(line[7:11])
    azm = float(line[12:15])
    pha = line[21:23].strip()
    wet = line[23]
    ptt = float(line[37:43])
    res = float(line[58:65])
    data = {
        "sta": sta,
        "dis": dis,
        "azm": azm,
        "pha": pha,
        "wet": wet,
        "ptt": ptt,
        "res": res,
    }
    return data


def summarize_single_events(single_events_file, sum_outfile):
    events_df = []
    header = {}
    phases = []
    eid = None
    with open(single_events_file) as f:
        for line in f:
            if "E V E N T" in line:
                header = {}
                phases = []
                eid = get_eid(line)
            elif "ERROR" in line:
                header.update({"eid": eid})
                info = {"eid": eid}
                phases.append(info)
                phases.append(header)
                events_df.extend(phases)
            elif "DATE  ORIGIN" in line:
                line = next(f)
                header = add_origin(line, header)
                header.update({"eid": eid})
            elif "ERX  ERY  ERZ" in line:
                line = next(f)
                header = add_statistics(line, header)
                header.update({"eid": eid})
            elif "STN  DIST" in line:
                line = next(f)
                line = next(f)
                previous = None
                while line.strip() and "$$$" not in line:
                    info = add_phase_info(line)
                    info.update({"eid": eid})
                    # An S printed below the P of the same station is a
                    # continuation card whose DIST/AZM columns hold AMX/PRX.
                    if (
                        previous is not None
                        and info["pha"].startswith("S")
                        and previous["sta"] == info["sta"]
                        and previous["pha"].startswith("P")
                    ):
                        info["dis"] = previous["dis"]
                        info["azm"] = previous["azm"]
                    previous = info
                    phases.append(info)
                    line = next(f)
                    if "$$$" in line:
                        phases.append(header)
                        events_df.extend(phases)
    events_df = DataFrame(events_df)
    # nst = number of stations with a usable (class < 4) phase, as in
    # summarize_catalog (VELEST's own NO column, kept as nobs, counts phases).
    if "sta" in events_df and "ort" in events_df:
        usable = events_df.dropna(subset=["sta"])
        usable = usable[
            to_numeric(usable["wet"], errors="coerce").fillna(4) < 4
        ]
        n_sta = usable.groupby("eid")["sta"].nunique()
        is_header = events_df["ort"].notna()
        events_df.loc[is_header, "nst"] = (
            events_df.loc[is_header, "eid"].map(n_sta).fillna(0)
        )
    events_df.to_csv(sum_outfile, index=False, float_format="%.3f")
    return events_df


def convert_single_events(single_events_file, cat_outfile, sum_outfile):
    events_df = summarize_single_events(single_events_file, sum_outfile)
    catalog = make_catalog(events_df)
    catalog.write(
        cat_outfile, format="NORDIC", high_accuracy=False, nordic_format="OLD"
    )
    return catalog


def cnv2csv(cnv_inpfile, csv_outfile):
    data = []
    with open(cnv_inpfile) as fo:
        for line in fo:
            if len(line) > 35 and line[25] in "SN" and line[35] in "EW":
                # Fields may be blank-padded and seconds may reach 60.
                date = line[:6].replace(" ", "0")
                hhmm = line[7:11].replace(" ", "0")
                ort = dt(
                    _two_digit_year(int(date[:2])),
                    int(date[2:4]),
                    int(date[4:6]),
                ) + timedelta(
                    hours=int(hhmm[:2]),
                    minutes=int(hhmm[2:]),
                    seconds=float(line[12:17]),
                )
                lat = float(line[18:25]) * (-1 if line[25] == "S" else 1)
                lon = float(line[27:35]) * (-1 if line[35] == "W" else 1)
                dep = float(line[37:43])
                mag = float(line[46:50])
                info = {
                    "ort": ort,
                    "lat": lat,
                    "lon": lon,
                    "dep": dep,
                    "mag": mag,
                }
                data.append(info)
    cnv_df = DataFrame(data)
    cnv_df.to_csv(
        csv_outfile,
        index=False,
    )


def reselect_best(config, stage_n):

    rootpath = Path(f"outputs/stage_{stage_n:02d}")

    # inputs
    single_events_file = rootpath / "single_events.out"

    # outputs
    cat_rel_outfile = rootpath / "relocated.out"
    cat_sum_outfile = rootpath / "relocated.csv"

    select_cat_outfile = rootpath / "select.dat"
    select_cnv_outfile = rootpath / "select.cnv"
    select_csv_outfile = rootpath / "select.csv"

    cat_rel = convert_single_events(
        single_events_file, cat_rel_outfile, cat_sum_outfile
    )

    # Location errors and residuals now come from VELEST: apply the
    # RESELECTION overrides on top of EVENT_SELECTION.
    reselection = config.model_copy(
        update={
            "EVENT_SELECTION": config.RESELECTION.apply(config.EVENT_SELECTION)
        }
    )
    cat_sel = select(
        reselection,
        cat_rel,
        cat_sum_outfile,
        select_cat_outfile,
        select_csv_outfile,
    )

    stations_inpfile = rootpath / "stations_sel.csv"
    catalog_to_cnv(
        cat_sel,
        stations_inpfile,
        rootpath,
        select_cnv_outfile,
        config.PICK_WEIGHTS,
    )


def perturb_catalog(
    config, run_n, cat_cnv_file, cat_csv_file, cnv_outfile, csv_outfile
):

    cat_csv = read_csv(cat_csv_file)

    # x = east, y = north; a degree of longitude shrinks with cos(latitude)
    lat_std = k2d(config.STAGE_04.dy_km)
    lon_std = k2d(config.STAGE_04.dx_km) / cos(radians(cat_csv["lat"].mean()))
    dep_std = config.STAGE_04.dz_km
    seed = config.STAGE_04.seed

    rng = random.default_rng(seed + run_n)

    # run 0 is the unperturbed reference run
    scale = 0.0 if run_n == 0 else 1.0
    nEvt = len(cat_csv)
    cat_csv["lat"] += scale * rng.normal(0, lat_std, size=nEvt)
    cat_csv["lon"] += scale * rng.normal(0, lon_std, size=nEvt)
    cat_csv["dep"] += scale * rng.normal(0, dep_std, size=nEvt)

    # Reflect depths perturbed above the surface: clipping them to exactly
    # 0 km would start VELEST at the free surface, where they tend to stay.
    cat_csv["dep"] = cat_csv["dep"].abs()

    with open(cat_cnv_file) as fo, open(cnv_outfile, "w") as go:
        eid = 0
        for line in fo:
            if len(line) > 35 and line[25] in "NS" and line[35] in "EW":
                lat = cat_csv.iloc[eid].lat
                lon = cat_csv.iloc[eid].lon
                dep = cat_csv.iloc[eid].dep
                line = (
                    line[:18]
                    + f"{lat:7.4f}"
                    + line[25:27]
                    + f"{lon:8.4f}"
                    + line[35:38]
                    + f"{dep:6.2f}"
                    + line[44:]
                )
                eid += 1
            go.write(line)
    cat_csv.to_csv(csv_outfile, index=False, float_format="%.3f")


def prepare_catalog(config, cat_inp=None, stage_n=1, run_n=None):

    if stage_n == 1:

        root = Path("outputs/stage_01")

        # inputs
        catalog_csv_inpfile = CATALOG_SUMMARY
        stations_inpfile = root / "stations_sel.csv"

        # output
        select_cat_outfile = root / "select.dat"
        select_csv_outfile = root / "select.csv"
        select_cnv_outfile = root / "select.cnv"

        cat_sel = select(
            config,
            cat_inp,
            catalog_csv_inpfile,
            select_cat_outfile,
            select_csv_outfile,
        )
        catalog_to_cnv(
            cat_sel,
            stations_inpfile,
            root,
            select_cnv_outfile,
            config.PICK_WEIGHTS,
            unconstrained_depth_below_km=(
                config.INPUTS.unconstrained_depth_below_km
            ),
            start_depth_km=config.INPUTS.start_depth_km,
        )

    elif stage_n == 2:
        # output
        root_out = Path("outputs/stage_02")
        stations_inpfile = root_out / "stations_sel.csv"
        original_cnv_outfile = root_out / "original.cnv"
        csv_outfile = root_out / "original.csv"
        catalog_to_cnv(
            cat_inp,
            stations_inpfile,
            root_out,
            original_cnv_outfile,
            config.PICK_WEIGHTS,
            unconstrained_depth_below_km=(
                config.INPUTS.unconstrained_depth_below_km
            ),
            start_depth_km=config.INPUTS.start_depth_km,
        )
        cnv2csv(original_cnv_outfile, csv_outfile)

    elif stage_n == 3:

        root = Path("outputs/stage_03")

        cat_sel = Path("outputs/stage_02") / "select.cnv"

        # output
        run_dir = root / f"run_{run_n:02d}"

        for n_model in range(1, config.STAGE_03.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            copy(cat_sel, model_dir)

    elif stage_n == 4:

        root = Path("outputs/stage_04")

        # inputs
        cat_cnv_file = Path("outputs/stage_02") / "select.cnv"
        cat_csv_file = Path("outputs/stage_02") / "select.csv"

        # output
        run_dir = root / f"run_{run_n:02d}"
        cnv_outfile = root / f"run_{run_n:02d}" / "perturb.cnv"
        csv_outfile = root / f"run_{run_n:02d}" / "perturb.csv"

        perturb_catalog(
            config, run_n, cat_cnv_file, cat_csv_file, cnv_outfile, csv_outfile
        )

    elif stage_n == 5:
        # output
        root_out = Path("outputs/stage_05")
        stations_inpfile = root_out / "stations_sel.csv"
        original_cnv_outfile = root_out / "original.cnv"
        csv_outfile = root_out / "original.csv"
        catalog_to_cnv(
            cat_inp,
            stations_inpfile,
            root_out,
            original_cnv_outfile,
            config.PICK_WEIGHTS,
            unconstrained_depth_below_km=(
                config.INPUTS.unconstrained_depth_below_km
            ),
            start_depth_km=config.INPUTS.start_depth_km,
        )
        cnv2csv(original_cnv_outfile, csv_outfile)
