import utils
import yaml

# YAML CONFIG 
try:
    with open("CONFIG.yaml", "r") as yaml_file:
        CONFIG = yaml.safe_load(yaml_file)
except FileNotFoundError:
    CONFIG = {}

# Base task
current_task = CONFIG.get("task", "question1")

# Defaults for every task:
TASKS_DEFAULTS = {
    "question1": {
        'SpaceCraft':'JUICE',
        'CentralBody':'Ganymede',
        'GravityModel':'Point-Mass gravitational model'
    },
    "question2": {
        'SpaceCraft':'JUICE',
        'CentralBody':'Ganymede',
        'GravityModel':'Point-Mass gravitational model'
    },
    "question3": {
        'SpaceCraft':'JUICE',
        'CentralBody':'Ganymede',
        'GravityModel':'Spherical Harmonic gravitational model'
    },
    "question4": {
        'SpaceCraft':'JUICE',
        'CentralBody':'Ganymede',
        'GravityModel':'Spherical Harmonic gravitational model'
    }
}

# Extracting only the defaults of the current task
DEFAULTS = TASKS_DEFAULTS.get(current_task, TASKS_DEFAULTS["question1"])

# SETUP folders and files:
if __name__ == "__PostProcessing__":
    utils.prepare_directories(current_task)