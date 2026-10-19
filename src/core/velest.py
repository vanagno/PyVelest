from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
from subprocess import run

from src.io.velest import read_velest_inversion
from src.io.catalog import cnv2csv
from tqdm import tqdm

# Files VELEST writes; removed before each run so that a failed run can
# never leave results from an earlier run behind.
VELEST_OUTPUTS = [
    "velout.inv",
    "velout.cnv",
    "velout.sta",
    "velout.res",
    "velout.smp",
    "velout.ray",
    "velout.sum",
    "single_events.out",
]
PARSED_OUTPUTS = {
    "attempts": "attempts.csv",
    "hypocenter_changes": "hypocenter_changes.csv",
    "hypocenter_summary": "hypocenter_summary.csv",
    "velocity_changes": "velocity_changes.csv",
    # read_velest_inversion's "velocity_states" table (one row per layer,
    # iteration and state) is what the model averaging reads as
    # velocity_summary.csv; keep these names for compatibility.
    "velocity_summary": "velocity_states.csv",
    "velocity_states": "velocity_summary.csv",
}


def _velest_exec(config) -> Path:
    exe = Path(config.GLOBAL_SETTINGS.velest_exec).resolve()
    if not exe.is_file():
        raise FileNotFoundError(
            f"VELEST executable not found: {exe} "
            "(build one with src/exc/build_velest.sh)"
        )
    if not os.access(exe, os.X_OK):
        raise PermissionError(
            f"VELEST executable is not executable: {exe}\n"
            f"    run: chmod +x {exe}\n"
            "    (or build one with src/exc/build_velest.sh and set "
            "GLOBAL_SETTINGS.velest_exec)"
        )
    return exe


def _execute(exe: Path, workdir: Path, single_event: bool = False) -> str | None:
    """Run VELEST inside ``workdir`` (which holds ``velest.cmn``).

    Returns None on success, otherwise a short error description.
    """
    for name in VELEST_OUTPUTS + list(PARSED_OUTPUTS.values()):
        (workdir / name).unlink(missing_ok=True)

    result = run(
        [str(exe)],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if single_event:
        # A STOP part-way through the catalogue also exits with code 0 and
        # leaves a truncated single_events.out: require the normal end.
        product = workdir / "single_events.out"
        ok = (
            product.exists()
            and product.stat().st_size > 0
            and "End of program VELEST" in (result.stdout or "")
        )
    else:
        product = workdir / "velout.inv"
        ok = product.exists() and "End of program VELEST" in product.read_text(
            errors="replace"
        )
    if result.returncode != 0 or not ok:
        tail = (result.stderr or result.stdout or "").strip().splitlines()[-3:]
        return f"exit code {result.returncode}; " + " | ".join(tail)
    return None


def _execute_and_parse(
    exe: Path, workdir: Path, need_velocity: bool = False
) -> str | None:
    error = _execute(exe, workdir)
    if error:
        return error
    try:
        results = read_velest_inversion(workdir / "velout.inv")
        for key, name in PARSED_OUTPUTS.items():
            results[key].to_csv(workdir / name, index=False)
    except Exception as exc:  # parsing problems must not pass silently
        return f"could not parse velout.inv: {exc}"
    if need_velocity and results["velocity_states"].empty:
        # VELEST rejected every step of an early iteration and stopped
        # before adjusting the velocities: the model is still the random
        # start. Without velocity_summary.csv the inversion is left out of
        # the averages of models and station delays.
        (workdir / PARSED_OUTPUTS["velocity_states"]).unlink()
        return "stopped before its first velocity adjustment (model unchanged)"
    return None


def _run_many(config, workdirs):
    exe = _velest_exec(config)
    n_jobs = config.GLOBAL_SETTINGS.n_jobs
    desc = "+++ Running VELEST ..."
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        errors = list(
            tqdm(
                pool.map(
                    lambda d: _execute_and_parse(exe, d, need_velocity=True),
                    workdirs,
                ),
                total=len(workdirs),
                desc=desc,
            )
        )
    failed = [(d, e) for d, e in zip(workdirs, errors) if e]
    for workdir, error in failed:
        print(f"!!! VELEST failed in {workdir}: {error}")
    if len(failed) == len(workdirs):
        raise RuntimeError("All VELEST runs failed; see messages above.")
    return failed


def _check(workdir, error):
    if error:
        raise RuntimeError(f"VELEST failed in {workdir}: {error}")


def run_velest(config, stage_n, run_n=None):

    if stage_n == 1:
        run_dir = Path("outputs/stage_01") / f"run_{run_n:02d}"
        n_models = config.SYNTHETIC_MODELS.n_models
        workdirs = [run_dir / f"model_{n:02d}" for n in range(1, n_models + 1)]
        _run_many(config, workdirs)

    elif stage_n == 2:
        root = Path("outputs/stage_02")
        _check(root, _execute(_velest_exec(config), root, single_event=True))

    elif stage_n == 3:
        run_dir = Path("outputs/stage_03") / f"run_{run_n:02d}"
        n_models = config.STAGE_03.n_models
        workdirs = [run_dir / f"model_{n:02d}" for n in range(1, n_models + 1)]
        _run_many(config, workdirs)

    elif stage_n == 4:
        run_dir = Path("outputs/stage_04") / f"run_{run_n:02d}"
        (run_dir / "relocated.csv").unlink(missing_ok=True)
        _check(run_dir, _execute_and_parse(_velest_exec(config), run_dir))
        cnv2csv(run_dir / "velout.cnv", run_dir / "relocated.csv")

    elif stage_n == 5:
        root = Path("outputs/stage_05")
        _check(root, _execute(_velest_exec(config), root, single_event=True))
