# import os, sys
# import subprocess as sp
# import time
# import shutil

from glob import glob

read_files = glob("./num_peaks_*/opt_parameters_all_island.csv")

with open("opt_parameters_all_island.csv", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
