import os
from . import ROOT_DIR
from .general import load_json


def get_config() -> dict:
    """Loads the config file."""
    config_obj = load_json(
        filepath=os.path.join(ROOT_DIR, "pipeline", "config.json")
    )
    return config_obj
