import sys, os, importlib

for package in ["scipy", "pandas", "numpy", "pickle","sklearn", "networkx", "argparse",
                "matplotlib", "seaborn", "warnings", "jsonpickle"]:
    try:
        globals()[package] = importlib.import_module(package)
    except ImportError:
        print("No %s package was installed" %package)
        sys.exit()
