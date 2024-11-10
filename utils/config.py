import os
from utils import ROOT_DIR
from utils.general import load_json_type_safe


def get_config() -> dict:
    """Loads the config file."""
    config_obj = load_json_type_safe(
        filepath=os.path.join(ROOT_DIR, "pipeline", "config.json"), return_type="dict"
    )
    return config_obj