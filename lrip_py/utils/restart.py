import json

def save_state(state, restart_out_file):
    """Saves state to a file."""
    with open(restart_out_file, "w") as f:
        json.dump(state, f)