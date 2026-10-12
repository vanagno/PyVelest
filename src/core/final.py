from src.io.model import read_vmodel, write_vmodel
from src.io.station import read_vstation
from src.model.Models import VELEST_CONFIG
from pathlib import Path
from pandas import DataFrame, read_csv
from shutil import copy
import yaml


def _stage05_velocity():
    """VELOCITY section of the stage-5 VELEST config (phase mode, Vp/Vs)."""
    with open(Path("configs") / "stage_05" / "velest.yaml") as f:
        return VELEST_CONFIG.model_validate(yaml.safe_load(f)).VELOCITY


def prepare_final_model(config):

    # inputs
    vmodel_inpfile = (
        Path("outputs/stage_05") / "model.mod"
    )

    # outputs
    vmodel_outfile = Path("outputs/stage_fn") / "model.mod"
    csv_outfile = Path("outputs/stage_fn") / "model.csv"

    # process
    vmodel = read_vmodel(vmodel_inpfile)
    write_vmodel(vmodel, vmodel_outfile)

    velocity = _stage05_velocity()
    vp = [layer.vel for layer in vmodel.layers_vp]
    if vmodel.layers_vs:
        vs = [layer.vel for layer in vmodel.layers_vs]
    elif velocity.phase_mode == 3:
        # S model VELEST used: Vs = Vp / (Vp/Vs)
        vs = [round(v / velocity.vp_vs_ratio, 3) for v in vp]
    else:
        vs = None
    vmodel_df = DataFrame(
        {
            "vp": vp,
            "depth": [layer.depth for layer in vmodel.layers_vp],
            "vs": vs,
        }
    )
    vmodel_df.to_csv(csv_outfile, index=False)


def prepare_final_catalog(config):

    # inputs
    reloc_csv_file = Path("outputs/stage_05") / "relocated.csv"
    reloc_cat_file = Path("outputs/stage_05") / "relocated.out"
    input_summary = Path("inputs") / "catalog.csv"
    model_file = Path("outputs/stage_05") / "model.mod"

    # outputs
    reloc_csv_outfile = Path("outputs/stage_fn") / "relocated.csv"
    reloc_cat_outfile = Path("outputs/stage_fn") / "relocated.out"
    events_outfile = Path("outputs/stage_fn") / "events.csv"

    copy(reloc_csv_file, reloc_csv_outfile)
    copy(reloc_cat_file, reloc_cat_outfile)

    # One row per event, with flags for depths that are not free solutions.
    reloc = read_csv(reloc_csv_file)
    events = reloc.dropna(subset=["ort"]).drop(
        columns=["sta", "dis", "azm", "pha", "wet", "ptt", "res"],
        errors="ignore",
    )
    # VELEST hypocentres pile up within ~0.1 km of layer boundaries (several
    # times the share expected by chance); flag them. Some events there are
    # genuine, so the flag is conservative.
    tops = [layer.depth for layer in read_vmodel(model_file).layers_vp]
    events["at_layer_boundary"] = events["dep"].apply(
        lambda d: any(abs(d - top) < 0.1 for top in tops)
    )
    min_erz = (
        config.RESELECTION.min_depth_error_km
        if config.RESELECTION.min_depth_error_km is not None
        else config.EVENT_SELECTION.min_depth_error_km
    )
    events["depth_error_near_zero"] = (
        events["erz"] < min_erz if min_erz is not None else False
    )
    events["at_surface"] = events["dep"] < 0.05
    threshold = config.INPUTS.unconstrained_depth_below_km
    if threshold is not None and input_summary.exists():
        inp = read_csv(input_summary).dropna(subset=["ort"])
        shallow = set(inp.loc[inp["dep"] < threshold, "eid"].astype(int))
        events["started_at_start_depth"] = events["eid"].astype(int).isin(
            shallow
        )
    events["depth_reliable"] = ~(
        events["at_layer_boundary"]
        | events["depth_error_near_zero"]
        | events["at_surface"]
    )
    events.to_csv(events_outfile, index=False, float_format="%.4f")
    print(
        f"+++ Final catalog: {len(events)} events, "
        f"{int(events['depth_reliable'].sum())} with a free depth "
        f"(see {events_outfile})"
    )


def prepare_final_stations(config):

    # inputs
    vstations_file = Path("outputs/stage_05") / "stations.sta"

    # outputs
    vstations_outfile = Path("outputs/stage_fn") / "stations.sta"
    stations_csvfile = Path("outputs/stage_fn") / "stations.csv"

    vstations = read_vstation(config, vstations_file)

    copy(vstations_file, vstations_outfile)

    stations_df = DataFrame([i.dict() for i in vstations.stations])

    velocity = _stage05_velocity()
    if velocity.phase_mode == 3:
        # VELEST inverts only P delays here and applies S delay =
        # P delay * Vp/Vs (velest.f, nsp=3); report the S delay it applied.
        stations_df["stcorr"] = stations_df["ptcorr"] * velocity.vp_vs_ratio

    stations_df.drop(["mn", "icc"], axis=1, inplace=True)
    stations_df.to_csv(stations_csvfile, index=False, float_format="%.4f")
