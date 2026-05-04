from pathlib import Path
import sys
from os.path import abspath, dirname, relpath
from inspect import getfile, getframeinfo, getsource
import os
import requests
import numpy as np
from math import trunc

from datetime import datetime


def prepare_directories(task="default"):

    """Creates plots and data folders for the specific task."""
    Path(f"../Plots/{task}").mkdir(parents=True, exist_ok=True)
    Path(f"../Data/{task}").mkdir(parents=True, exist_ok=True)


def plot_filename(base_name, task="default"):

    """Creates saving path for the task."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    return f'../Plots/{task}/{base_name}_{timestamp}.pdf'



# def plot_filename(base_name):
#     timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
#     return f'plots/{base_name}_{timestamp}.png'

