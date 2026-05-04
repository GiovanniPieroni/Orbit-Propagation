import argparse
from config import CONFIG, TASKS_DEFAULTS

parser = argparse.ArgumentParser()

# Choose mode
parser.add_argument("--task", type=str) 
parser.add_argument("--only_plotting", type=bool) 



# Specific of tasks
parser.add_argument("--SpaceCraft", type=str) 
parser.add_argument("--CentralBody", type=str) 
parser.add_argument("--GravityModel", type=str) 


CLI = parser.parse_args()

def pick(param, override_task = None):

    active_task = override_task or CLI.task or CONFIG.get("task", "question1")

    # Special case to request the current task
    if param == "task":
        # if CLI.task: return CLI.task
        # return CONFIG.get("task", "question1")
        return active_task

    # 1. Terminal inputs have priority
    if hasattr(CLI, param) and getattr(CLI, param) not in [None, False]:
        return getattr(CLI, param)
    
    # 2. YAML config file defaults
    if param in CONFIG:
        return CONFIG[param]
    
    # 3. DEFAUTLS in config.py
    if active_task in TASKS_DEFAULTS and param in TASKS_DEFAULTS[active_task]:
        return TASKS_DEFAULTS[active_task][param]
        
    return None # O Feedback value if parameter doesn't exist