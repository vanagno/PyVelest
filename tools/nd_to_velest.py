"""Convert the crust of a TauP ``.nd`` model into a VELEST starting model.

The crustal rows (everything before the ``mantle`` line) define a velocity
profile. With ``--interp linear`` (default, TauP's own semantics) velocity
varies linearly between consecutive rows and a discontinuity is two rows at
the same depth; with ``--interp step`` each row is the top of a
constant-velocity layer. Each thin VELEST layer gets the thickness-weighted
average of that profile, so the averages of the source model are preserved. Optionally the
profile is smoothed (running mean over ``--smooth`` layers), because
PyVelest's random-model generator perturbs velocity *increments* and is
biased for step-like models. Velocities are then forced to increase by at
least ``--min-increment`` km/s per layer, as the generator requires.

The last layer (half-space) starts at ``--moho`` with velocity ``--vmoho``.
The first layer starts at ``--top`` (it must lie above the highest station).

Example::

    python tools/nd_to_velest.py "inputs/My inputs/kiratzi87.nd" \\
        --tops -3 0 1 2 3 4 6 8 10 13 16 20 25 --moho 30 \\
        --vtop 3.5 --min-increment 0.05 --out inputs/vmodel_kiratzi87.dat

The first layer (-3 to 0 km here) only covers station elevations; giving
every starting model the same, low --vtop keeps it from acting as a
different velocity floor in each model (low-velocity layers are not
allowed, and the first layer is heavily damped).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read_nd_crust(path):
    """[(depth_top_km, vp, vs), ...] for the rows above the 'mantle' line."""
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.strip().lower().startswith("mantle"):
                break
            depth, vp, vs = (float(v) for v in line.split()[:3])
            rows.append((depth, vp, vs))
    if not rows:
        raise ValueError(f"{path}: no crustal rows found")
    return rows


def layer_average(rows, top, bottom, column=1, interp="linear"):
    """Thickness-weighted average of the profile in [top, bottom].

    Outside the listed depths the first/last value holds.
    """
    if interp == "linear":
        z = np.linspace(top, bottom, 2001)
        depths = [r[0] for r in rows]
        values = [r[column] for r in rows]
        return float(np.interp(z, depths, values).mean())
    depths = [r[0] for r in rows] + [np.inf]
    total = weighted = 0.0
    for i, row in enumerate(rows):
        lo, hi = max(top, depths[i]), min(bottom, depths[i + 1])
        if hi > lo:
            weighted += row[column] * (hi - lo)
            total += hi - lo
    return weighted / total


def build_model(
    rows, tops, moho, smooth=0, min_increment=0.02, interp="linear", vtop=None
):
    velocities = []
    for i, top in enumerate(tops):
        bottom = tops[i + 1] if i + 1 < len(tops) else moho
        if i == 0 and vtop is not None:
            velocities.append(vtop)
            continue
        velocities.append(
            layer_average(rows, max(top, rows[0][0]), bottom, interp=interp)
        )
    v = np.asarray(velocities, dtype=float)
    if smooth > 1:
        # smooth below the first layer only (it may be a fixed dummy layer)
        body = v[1:]
        padded = np.pad(body, (smooth // 2, smooth // 2), mode="edge")
        v[1:] = np.convolve(padded, np.ones(smooth) / smooth, mode="valid")
    for i in range(1, len(v)):
        v[i] = max(v[i], v[i - 1] + min_increment)
    return [float(x) for x in np.round(v, 2)]


def write_velest_model(path, name, tops, velocities, moho, vmoho, dmp):
    first_dmp, mid_dmp, last_dmp = dmp
    depths = list(tops) + [moho]
    vels = list(velocities) + [vmoho]
    n = len(depths)
    lines = [
        f"{name[:40]:40}",
        f"{n:3d}        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)",
    ]
    for i, (v, d) in enumerate(zip(vels, depths)):
        damp = first_dmp if i == 0 else last_dmp if i == n - 1 else mid_dmp
        line = f"{v:5.2f}     {d:7.2f}  {damp:7.3f}"
        if i == 0:
            line += "           P-VELOCITY MODEL"
        lines.append(line)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("nd")
    parser.add_argument("--tops", nargs="+", type=float, required=True)
    parser.add_argument("--moho", type=float, required=True)
    parser.add_argument("--vmoho", type=float, default=8.0)
    parser.add_argument(
        "--vtop",
        type=float,
        help="velocity of the first layer (e.g. a common, damped layer above "
        "sea level); default: the profile average over that layer",
    )
    parser.add_argument("--smooth", type=int, default=0)
    parser.add_argument(
        "--interp",
        choices=["linear", "step"],
        default="linear",
        help="linear: TauP .nd semantics (default); step: rows are layer tops",
    )
    parser.add_argument("--min-increment", type=float, default=0.02)
    parser.add_argument("--damping", nargs=3, type=float, default=[50, 1, 10])
    parser.add_argument("--name")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tops = sorted(args.tops)
    if tops[-1] >= args.moho:
        raise ValueError("--moho must be deeper than the last layer top")
    rows = read_nd_crust(args.nd)
    vels = build_model(
        rows,
        tops,
        args.moho,
        args.smooth,
        args.min_increment,
        args.interp,
        args.vtop,
    )
    if vels[-1] >= args.vmoho:
        raise ValueError("--vmoho must exceed the deepest crustal velocity")
    name = args.name or Path(args.nd).stem.upper()
    write_velest_model(
        args.out, name, tops, vels, args.moho, args.vmoho, args.damping
    )
    vpvs = [r[1] / r[2] for r in rows if r[2] > 0]
    print(f"{args.out}: {len(vels) + 1} layers; source Vp/Vs "
          f"{min(vpvs):.3f}-{max(vpvs):.3f}")
    for d, v in zip(tops + [args.moho], vels + [args.vmoho]):
        print(f"  {d:6.1f} km  {v:5.2f} km/s")


if __name__ == "__main__":
    main()
