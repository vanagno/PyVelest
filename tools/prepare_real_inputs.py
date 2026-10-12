"""Convert REAL / hypoSA outputs into PyVelest inputs.

Builds two files that PyVelest reads through ``configs/user.yaml``:

* a station table (``code,lat,lon,elv,net``; elevation in metres) from a
  REAL station file (``lon lat net sta chan elev_km`` per line);
* a HypoDD phase file (ObsPy format ``HYPODDPHA``) that combines

  - the hypoSA hypocentres and origin times (``catalogSA_allday.txt``),
  - the PhaseNet pick probabilities from the REAL phase file
    (``phase_allday.txt`` or ``ioannina_final.txt``), stored in the weight
    column so PyVelest can turn them into VELEST weight classes,
  - the hypoSA travel-time residual in the ``rms`` header field.

Travel times are re-referenced to the hypoSA origin time, so absolute pick
times are unchanged.

Example::

    python tools/prepare_real_inputs.py \\
        --stations "inputs/My inputs/ioannina/ioan_stations.dat" \\
        --phases "inputs/My inputs/ioannina_final.txt" \\
        --sa-catalog "inputs/My inputs/ioannina/catalogSA_allday.txt" \\
        --sa-phases "inputs/My inputs/ioannina/phaseSA_allday_kiratzi87.txt" \\
        --out-stations inputs/stations.csv \\
        --out-catalog inputs/catalog.pha
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path


def read_real_stations(path):
    stations = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            lon, lat, net, sta, _chan, elev_km = line.split()[:6]
            stations.append(
                {
                    "code": sta,
                    "lat": float(lat),
                    "lon": float(lon),
                    "elv": int(round(float(elev_km) * 1000.0)),
                    "net": net,
                }
            )
    return stations


def read_pha(path):
    """Return a list of (header_fields, picks) in file order."""
    events = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                events.append((line.split()[1:], []))
            elif line.strip():
                sta, tt, weight, phase = line.split()[:4]
                events[-1][1].append((sta, float(tt), float(weight), phase))
    return events


def header_time(fields):
    yr, mo, dy, hr, mn = (int(v) for v in fields[:5])
    return datetime(yr, mo, dy, hr, mn) + timedelta(seconds=float(fields[5]))


def read_sa_catalog(path):
    """hypoSA catalog: yr mo dy hr mn sec lat lon dep mag nofps gap res id."""
    events = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.split()
            events.append(
                {
                    "time": header_time(p[:6]),
                    "lat": float(p[6]),
                    "lon": float(p[7]),
                    "dep": float(p[8]),
                    "mag": float(p[9]),
                    "nphase": int(p[10]),
                    "gap": float(p[11]),
                    "rms": float(p[12]),
                }
            )
    return events


def check_sa_phases(sa_events, sa_phase_file):
    """Cross-check hypoSA catalog order against the hypoSA phase file."""
    sa_phases = read_pha(sa_phase_file)
    if len(sa_phases) != len(sa_events):
        raise ValueError(
            f"{sa_phase_file}: {len(sa_phases)} events, "
            f"catalog has {len(sa_events)}"
        )
    for i, ((fields, _), event) in enumerate(zip(sa_phases, sa_events)):
        dt = abs((header_time(fields[:6]) - event["time"]).total_seconds())
        if (
            dt > 0.01
            or abs(float(fields[6]) - event["lat"]) > 1e-3
            or abs(float(fields[7]) - event["lon"]) > 1e-3
        ):
            raise ValueError(f"hypoSA files disagree at event #{i + 1}")


def write_stations(stations, outfile):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with outfile.open("w") as f:
        f.write("code,lat,lon,elv,net\n")
        for s in stations:
            f.write(
                f"{s['code']},{s['lat']:.4f},{s['lon']:.4f},"
                f"{s['elv']:d},{s['net']}\n"
            )


def write_catalog(real_events, sa_events, outfile, max_shift_s=10.0):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    n_sa = 0
    with outfile.open("w") as f:
        for i, (fields, picks) in enumerate(real_events):
            evid = int(fields[-1])
            mag = float(fields[9])
            ot_real = header_time(fields[:6])
            if sa_events is not None:
                sa = sa_events[i]
                shift = (ot_real - sa["time"]).total_seconds()
                if abs(shift) > max_shift_s:
                    raise ValueError(
                        f"Event {evid}: REAL and hypoSA origin times differ "
                        f"by {shift:.2f} s; files are not in the same order"
                    )
                ot, lat, lon, dep, rms = (
                    sa["time"],
                    sa["lat"],
                    sa["lon"],
                    sa["dep"],
                    sa["rms"],
                )
                n_sa += 1
            else:
                shift = 0.0
                ot = ot_real
                lat, lon, dep = (float(v) for v in fields[6:9])
                rms = 0.0
            sec = ot.second + ot.microsecond / 1e6
            f.write(
                f"# {ot:%Y %m %d %H %M} {sec:06.3f} {lat:.4f} {lon:.4f} "
                f"{dep:.3f} {mag:.2f} 0.0 0.0 {rms:.3f} {evid}\n"
            )
            for sta, tt, prob, phase in picks:
                f.write(f"{sta} {tt + shift:.4f} {prob:.4f} {phase}\n")
    return n_sa


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--stations", required=True)
    parser.add_argument("--phases", required=True)
    parser.add_argument("--sa-catalog")
    parser.add_argument("--sa-phases")
    parser.add_argument("--out-stations", default="inputs/stations.csv")
    parser.add_argument("--out-catalog", default="inputs/catalog.pha")
    args = parser.parse_args()

    stations = read_real_stations(args.stations)
    write_stations(stations, args.out_stations)

    real_events = read_pha(args.phases)
    sa_events = None
    if args.sa_catalog:
        sa_events = read_sa_catalog(args.sa_catalog)
        if len(sa_events) != len(real_events):
            raise ValueError(
                f"{args.sa_catalog}: {len(sa_events)} events, "
                f"{args.phases} has {len(real_events)}"
            )
        if args.sa_phases:
            check_sa_phases(sa_events, args.sa_phases)

    n_sa = write_catalog(real_events, sa_events, args.out_catalog)
    n_picks = sum(len(p) for _, p in real_events)
    print(f"{len(stations)} stations -> {args.out_stations}")
    print(
        f"{len(real_events)} events ({n_sa} with hypoSA hypocentres), "
        f"{n_picks} picks -> {args.out_catalog}"
    )


if __name__ == "__main__":
    main()
