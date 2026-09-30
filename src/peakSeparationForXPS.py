import os
import sys
import subprocess as sp

# import time
# from glob import glob

# import numpy as np
# import pandas as pd
import time
import datetime
import csv

PYTHON = "python3"
toolname = "xps-ps-cv-20240712"
bindir = os.path.abspath(os.path.dirname(__file__))
PY_FLOW = f"{bindir}/flow.py"


def set_recipe_data():
    recipe_data = {"bg": "", "noise": "-g", "pptx": "-x", "optmode": "simple", "dev": ""}

    recipe_db = {
        "bg": {
            "auto": "",
            "shirley": "-s",
            "linear": "-l"
        },
        "noise": {
            "gauss": "-g",
            "poisson": "-p"
        },
        "pptx": {
            "1": "-x",
            "0": ""
        },
        "optmode": {
            "simple": "",
            "full": "--full"
        },
        "dev": {
            "1": "-d",
            "0": ""
        }
    }  # fmt: skip

    recipe_file = "recipe.txt"

    if os.path.exists(recipe_file):
        # Output recipe_file as it is
        with open(recipe_file, "r") as f:
            print("\nrecipe file = \n--------------------\n" + f.read() + "\n--------------------\n")

        # Read recipe file to get background option and noise option as dict format
        with open(recipe_file, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                row[0] = row[0].lower()
                row[1] = row[1].lower()

                if row[0] in recipe_db.keys():
                    if row[1] in recipe_db[row[0]].keys():
                        recipe_data[row[0]] = recipe_db[row[0]][row[1]]
                    else:
                        raise ValueError(f"Invalid key, value in {recipe_file}: {row[0]}, {row[1]}")
                else:
                    raise ValueError(f"Invalid key in {recipe_file}: {row[0]}")

    else:
        print(f"{recipe_file} not found. Set default recipe_data.")

    return recipe_data


if __name__ == "__main__":

    print(f"# Tool name = {toolname}")
    print(f"# Start time = {datetime.datetime.now()}")
    startTimeTmp = time.time()  # start time

    print(sys.argv)

    if len(sys.argv) != 2:
        print(f"Error! Usage: {sys.argv[0]} <spectrum data file>")
        sys.exit(1)

    spectrumfile = sys.argv[1]

    recipe_data = set_recipe_data()

    print("recipe_data = ", recipe_data)

    cmd = [PYTHON, PY_FLOW, spectrumfile, recipe_data["bg"], recipe_data["noise"], recipe_data["pptx"], recipe_data["optmode"], recipe_data["dev"]]
    print(cmd, flush=True)

    status = sp.run(cmd)
    print(status)

    print(f"# End time = {datetime.datetime.now()}")
    td = time.time() - startTimeTmp
    print(toolname + " run time is %.2f (sec) =  %s" % (td, datetime.timedelta(seconds=td)))
