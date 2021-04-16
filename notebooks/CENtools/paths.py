import os, pathlib
 
# Specify the paths to import them from the modules
path = str(pathlib.Path(__file__).absolute())
path = path.replace("paths.py","")
raw_data_path = path + "data/raw_data/"
curated_data_path = path + "data/curated_data/"
object_path = path + "data/objects/"
