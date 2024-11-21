import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
print(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config
from utils.general import resolve_symlink

config_obj = get_config()
print("Yay")
dl_dir = config_obj["relevant_paths"]["downloads"]
print(dl_dir)