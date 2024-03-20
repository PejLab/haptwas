import os
import subprocess

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__),
                   os.path.pardir)



def run(prediction_fname):
    script_name = os.path.join(SCRIPTS_DIR,
                               "hello_world.R")
    subprocess.run(["Rscript", script_name])
