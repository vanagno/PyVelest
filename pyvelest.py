from pathlib import Path
from src.io.config import read_config
from src.io.workspace import prepare_workspace
from src.io.station import prepare_station, finalized_inverted_stations
from src.io.model import (
    prepare_model,
    finalized_inverted_models,
)
from src.io.catalog import summarize_catalog, prepare_catalog, reselect_best
from src.io.cmn import prepare_cmn
from src.core.velest import run_velest
from src.core.visualize import (
    plot_events,
    plot_models,
    plot_statistics,
    plot_dislocations,
)
from src.core.final import (
    prepare_final_model,
    prepare_final_catalog,
    prepare_final_stations,
)


class Main:
    def __init__(self, config_path: str, stages=(1, 2, 3, 4, 5)):
        self.config = read_config(config_path)
        # Later stages need the outputs of earlier ones: only wipe the
        # workspace when the run starts from stage 1.
        clean = self.config.GLOBAL_SETTINGS.clean_workspace
        if clean and min(stages) > 1:
            print(
                "+++ clean_workspace ignored: --stages does not include "
                "stage 1, earlier outputs are kept"
            )
            clean = False
        prepare_workspace(self.config, clean=clean)

    def load_catalog(self):
        self.catalog = summarize_catalog(self.config)

    def stage_01(self):
        n_runs = self.config.STAGE_01.n_runs
        for n in range(1, n_runs + 1):
            print(f"+++ Execution number {n} from stage_01 ...")
            prepare_station(self.config, stage_n=1, run_n=n)
            prepare_model(self.config, stage_n=1, run_n=n)
            prepare_catalog(self.config, self.catalog, stage_n=1, run_n=n)
            prepare_cmn(self.config, stage_n=1, run_n=n)
            run_velest(self.config, stage_n=1, run_n=n)
            finalized_inverted_models(self.config, stage_n=1, run_n=n)
            finalized_inverted_stations(self.config, stage_n=1, run_n=n)
            plot_models(self.config, stage_n=1, run_n=n)
            plot_statistics(self.config, stage_n=1, run_n=n)
        plot_events(self.config)

    def stage_02(self):
        print("+++ Executing stage_02 ...")
        prepare_station(self.config, stage_n=2)
        prepare_model(self.config, stage_n=2)
        prepare_catalog(self.config, self.catalog, stage_n=2)
        prepare_cmn(self.config, stage_n=2)
        run_velest(self.config, stage_n=2)
        reselect_best(self.config, stage_n=2)

    def stage_03(self):
        n_runs = self.config.STAGE_03.n_runs
        for n in range(1, n_runs + 1):
            print(f"+++ Execution number {n} from stage_03 ...")
            prepare_station(self.config, stage_n=3, run_n=n)
            prepare_model(self.config, stage_n=3, run_n=n)
            prepare_catalog(self.config, stage_n=3, run_n=n)
            prepare_cmn(self.config, stage_n=3, run_n=n)
            run_velest(self.config, stage_n=3, run_n=n)
            finalized_inverted_models(self.config, stage_n=3, run_n=n)
            finalized_inverted_stations(self.config, stage_n=3, run_n=n)
            plot_models(self.config, stage_n=3, run_n=n)
            plot_statistics(self.config, stage_n=3, run_n=n)

    def stage_04(self):
        # Stability test (Kissling et al., 1994). run_00 inverts the
        # unperturbed events with the same model and station delays; the
        # perturbed runs must return to its positions. Stage 4 is a
        # diagnostic: a failed run is reported, not fatal.
        n_runs = self.config.STAGE_04.n_runs
        skipped = []
        for n in range(0, n_runs + 1):
            label = "reference run" if n == 0 else f"execution number {n}"
            print(f"+++ Stage_04 {label} ...")
            figure = Path("figures") / f"dislocation_{n:02d}.png"
            figure.unlink(missing_ok=True)
            prepare_station(self.config, stage_n=4, run_n=n)
            prepare_model(self.config, stage_n=4, run_n=n)
            prepare_catalog(self.config, stage_n=4, run_n=n)
            prepare_cmn(self.config, stage_n=4, run_n=n)
            try:
                run_velest(self.config, stage_n=4, run_n=n)
                if n > 0:
                    plot_dislocations(self.config, stage_n=4, run_n=n)
            except (RuntimeError, ValueError) as exc:
                if n == 0:
                    print(
                        f"!!! stage_04 reference run failed ({exc}); "
                        "stability test skipped"
                    )
                    return
                print(f"!!! stage_04 run {n:02d} skipped: {exc}")
                skipped.append(n)
        if skipped:
            print(f"!!! stage_04 runs skipped: {skipped}")

    def stage_05(self):
        prepare_station(self.config, stage_n=5)
        prepare_model(self.config, stage_n=5)
        prepare_catalog(self.config, self.catalog, stage_n=5)
        prepare_cmn(self.config, stage_n=5)
        run_velest(self.config, stage_n=5)   
        reselect_best(self.config, stage_n=5)
        prepare_final_model(self.config)
        prepare_final_catalog(self.config)
        prepare_final_stations(self.config)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run the PyVELEST stages.")
    parser.add_argument(
        "--config",
        default=(Path("configs") / "user.yaml").as_posix(),
        help="user configuration file (default: configs/user.yaml)",
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        type=int,
        choices=[1, 2, 3, 4, 5],
        default=[1, 2, 3, 4, 5],
        help="stages to run, in order (default: all)",
    )
    args = parser.parse_args()

    app = Main(args.config, stages=args.stages)
    if {1, 2, 5} & set(args.stages):
        app.load_catalog()
    for stage_n in sorted(args.stages):
        getattr(app, f"stage_{stage_n:02d}")()
