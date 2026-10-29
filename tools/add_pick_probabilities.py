"""Use a located HypoDD phase file as PyVelest input, with PhaseNet pick
probabilities.

A HypoDD phase file written after a location run (for example LOC-FLOW's
HypoInverse step, NonLinLoc or hypoDD) already holds the hypocentres, origin
times, magnitudes, errors (eh, ez), RMS and picks that PyVelest reads, but
its pick weights are usually all 1.0. PyVelest turns the weights into VELEST
weight classes (PICK_WEIGHTS) and counts high-quality P picks from them
(EVENT_SELECTION.hq_p_max_class), so with weights of 1.0 every pick is a
best-class pick.

This tool keeps the located file as it is and puts back the PhaseNet
probability of each pick: that of the REAL pick with the same station and
phase whose absolute time is within --tolerance. The REAL picks come from
REAL's phase file (phase_allday.txt) or from a catalog.pha built from it with
prepare_real_inputs.py. Picks without such a REAL pick are dropped, unless
--unmatched-probability gives them one; events left without picks are
dropped. Station codes can be renamed in both files (VELEST reads codes of
at most 4 characters).

Example::

    python tools/add_pick_probabilities.py \\
        --located hypoDDevia.pha \\
        --phases phase_allday.txt \\
        --rename SAGMA=SAGM \\
        --out inputs/Evia/catalog_located.pha
"""

from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter
from datetime import datetime
from pathlib import Path

from prepare_real_inputs import header_time, parse_renames, read_pha

EPOCH = datetime(1970, 1, 1)


def pick_seconds(fields, tt):
    """Absolute pick time (seconds since 1970) from an event header and a
    travel time."""
    return (header_time(fields[:6]) - EPOCH).total_seconds() + tt


def pick_index(events, rename):
    """{(station, phase): (sorted pick times, probabilities)}."""
    picks = {}
    for fields, event_picks in events:
        for sta, tt, prob, phase in event_picks:
            key = (rename.get(sta, sta), phase.upper()[0])
            picks.setdefault(key, []).append((pick_seconds(fields, tt), prob))
    index = {}
    for key, values in picks.items():
        values.sort()
        index[key] = ([t for t, _ in values], [p for _, p in values])
    return index


def probability(index, sta, phase, t, tolerance):
    """Probability of the nearest REAL pick within tolerance, or None."""
    if (sta, phase) not in index:
        return None
    times, probs = index[(sta, phase)]
    i = bisect_left(times, t)
    nearest = min(
        (j for j in (i - 1, i) if 0 <= j < len(times)),
        key=lambda j: abs(times[j] - t),
        default=None,
    )
    if nearest is None or abs(times[nearest] - t) > tolerance:
        return None
    return probs[nearest]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--located", required=True, help="located HypoDD phase file"
    )
    parser.add_argument(
        "--phases",
        required=True,
        help="REAL phase file with PhaseNet probabilities in the weight "
        "column (phase_allday.txt or a catalog.pha built from it)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.02,
        help="largest difference of absolute pick times, s (default 0.02)",
    )
    parser.add_argument(
        "--unmatched-probability",
        type=float,
        help="keep picks without a REAL pick, with this probability "
        "(default: drop them)",
    )
    parser.add_argument(
        "--rename",
        action="append",
        default=[],
        metavar="OLD=NEW",
        help="rename a station in both files (VELEST codes have at most 4 "
        "characters); repeat for several stations",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rename = parse_renames(args.rename)
    index = pick_index(read_pha(args.phases), rename)
    located = read_pha(args.located)

    n_matched, unmatched, dropped_events = 0, Counter(), 0
    probabilities = {"P": [], "S": []}
    codes = set()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for fields, picks in located:
            lines = []
            for sta, tt, _, phase in picks:
                sta = rename.get(sta, sta)
                prob = probability(
                    index,
                    sta,
                    phase.upper()[0],
                    pick_seconds(fields, tt),
                    args.tolerance,
                )
                if prob is None:
                    unmatched[sta] += 1
                    prob = args.unmatched_probability
                    if prob is None:
                        continue
                else:
                    n_matched += 1
                    probabilities.setdefault(phase.upper()[0], []).append(prob)
                codes.add(sta)
                lines.append(f"{sta} {tt:.4f} {prob:.4f} {phase}\n")
            if not lines:
                dropped_events += 1
                continue
            f.write("# " + " ".join(fields) + "\n")
            f.writelines(lines)

    too_long = sorted(c for c in codes if len(c) > 4)
    if too_long:
        out.unlink()
        hint = " ".join(f"--rename {c}=NEW" for c in too_long)
        raise ValueError(
            f"VELEST station codes have at most 4 characters: {too_long}. "
            f"Give them unique codes of up to 4 characters with {hint}"
        )

    n_picks = n_matched + sum(unmatched.values())
    print(
        f"{len(located) - dropped_events} of {len(located)} events, "
        f"{n_matched} of {n_picks} picks with a REAL probability -> {out}"
    )
    if unmatched:
        action = (
            "dropped"
            if args.unmatched_probability is None
            else f"kept with probability {args.unmatched_probability}"
        )
        print(
            f"{sum(unmatched.values())} picks without a REAL pick ({action}):"
            + "".join(f" {s} {n}" for s, n in unmatched.most_common())
        )
    if dropped_events:
        print(f"{dropped_events} events had no picks left and were dropped")
    for phase, values in probabilities.items():
        if values:
            values.sort()
            print(
                f"{phase} probabilities: median {values[len(values) // 2]:.2f}, "
                f"{sum(v >= 0.70 for v in values)} of {len(values)} >= 0.70"
            )


if __name__ == "__main__":
    main()
