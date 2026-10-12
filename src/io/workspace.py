from pathlib import Path


def prepare_workspace(config, clean=None):
    """Create the output folders; delete old outputs if ``clean`` (default:
    GLOBAL_SETTINGS.clean_workspace)."""
    if clean is None:
        clean = config.GLOBAL_SETTINGS.clean_workspace

    outputs = Path("outputs")
    outputs.mkdir(parents=True, exist_ok=True)

    for directory in [
        "stage_01",
        "stage_02",
        "stage_03",
        "stage_04",
        "stage_05",
        "stage_fn",
    ]:

        target = Path("outputs") / directory
        target.mkdir(parents=True, exist_ok=True)

        if directory in ["stage_01", "stage_03"]:
            if directory == "stage_01":
                n_runs = config.STAGE_01.n_runs
                n_models = config.SYNTHETIC_MODELS.n_models
            else:
                n_runs = config.STAGE_03.n_runs
                n_models = config.STAGE_03.n_models
            for n_run in range(1, n_runs + 1):
                run_dir = target / f"run_{n_run:02d}"
                run_dir.mkdir(parents=True, exist_ok=True)
                sum_dir = run_dir / "summary"
                sum_dir.mkdir(parents=True, exist_ok=True)
                for n_model in range(1, n_models + 1):
                    model_dir = run_dir / f"model_{n_model:02d}"
                    model_dir.mkdir(parents=True, exist_ok=True)

        if directory in ["stage_04"]:
            # run_00 is the unperturbed reference run
            for n_run in range(0, config.STAGE_04.n_runs + 1):
                run_dir = target / f"run_{n_run:02d}"
                run_dir.mkdir(parents=True, exist_ok=True)

                    
    figures = Path("figures")
    figures.mkdir(parents=True, exist_ok=True)

    if clean:
        for folder in [outputs, figures]:
            for item in folder.rglob("*"):
                if item.is_file() and item.exists():
                    item.unlink()
