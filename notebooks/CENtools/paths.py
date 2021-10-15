import os
from pathlib import Path

# Specify the paths to import them from the modules
path = str(Path(os.path.abspath(__file__)).parents[0]) + "/"
raw_data_path = path + "data/raw_data/"
curated_data_path = path + "data/curated_data/"
object_path = path + "data/objects/"
