import subprocess
import time
import sys
from pathlib import Path

# --- User variables ---
N = 5  # Increased to 5 to generate a meaningful 95% Confidence Interval band
CC3D_PROJECT = str(Path("../Angiogenesis/Angiogenesis.cc3d").resolve())
CONFIG_PATH = Path("../Angiogenesis/Simulation/AngiogenesisConfig.py")
RESULTS_DIR = Path("../results")

# Ensure the directory exists
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILENAME = "angiogenesis_metrics.csv"

# Configs to run
CONFIGS = [
    "no_endo",  
    "no_hif1a",
    "branching_tuning",
]

def clear_old_results():
    """Delete old CSVs and plots in the target config folders to avoid mixing old and new runs."""
    for config in CONFIGS:
        config_dir = RESULTS_DIR / config
        if config_dir.exists():
            for file in config_dir.glob("*.csv"):
                file.unlink()
            for file in config_dir.glob("*.png"):
                file.unlink()
            print(f"Cleared old data for config: {config}")

def set_selected_preset(config_name):
    """Edit SELECTED_PRESET in AngiogenesisConfig.py to config_name."""
    lines = CONFIG_PATH.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("SELECTED_PRESET"):
            lines[i] = f'SELECTED_PRESET = "{config_name}"'
            break
    CONFIG_PATH.write_text("\n".join(lines) + "\n")

CC3D_EXECUTABLE = str(Path("~/CompuCell3D/runScript.command").expanduser())

def run_simulation():
    """Run the CC3D simulation headlessly."""
    command = [
        "bash",
        CC3D_EXECUTABLE,
        "--input", CC3D_PROJECT
    ]
    result = subprocess.run(command, shell=False)  
    if result.returncode != 0:
        raise RuntimeError("CC3D run failed!")

def move_and_tag_output(config_name, run_id):
    """Move the output CSV to a new file with config and run info."""
    config_dir = RESULTS_DIR / config_name
    config_dir.mkdir(parents=True, exist_ok=True) 

    src = RESULTS_DIR / OUTPUT_FILENAME
    dst = config_dir / f"angiogenesis_metrics_run{run_id:02d}.csv"
    if not src.exists():
        raise FileNotFoundError(f"Expected output file not found: {src}")
    
    add_run_column(src, dst, run_id)
    src.unlink()  

def add_run_column(src_path, dst_path, run_id):
    """Copy CSV and add a run column."""
    import pandas as pd
    df = pd.read_csv(src_path)
    df.insert(0, "run", run_id)
    df.to_csv(dst_path, index=False)

def main():
    clear_old_results()  # Clean the slate before we start!
    
    for config_name in CONFIGS:
        print(f"\n--- Running config: {config_name} ---")
        for run_id in range(1, N+1):
            print(f"  Run {run_id}/{N}")
            set_selected_preset(config_name)
            run_simulation()
            move_and_tag_output(config_name, run_id)
            time.sleep(1) 
    print("\nAll runs complete.")

    # Safely gather only the CSVs from our targeted config directories
    csv_files = []
    for config in CONFIGS:
        config_dir = RESULTS_DIR / config
        if config_dir.exists():
            csv_files.extend([str(f) for f in config_dir.glob("*.csv")])

    if not csv_files:
        print("No CSV files found to plot.")
        return

    # Auto-run plot_multi_run_angiogenesis.py
    # Pass the base RESULTS_DIR; the plotting script will handle placing files in the subfolders.
    plot_command = [
        sys.executable, "plot_multi_run_angiogenesis.py",
        "--inputs", *csv_files,
        "--output-dir", str(RESULTS_DIR)
    ]
    print("\nRunning plot command...")
    subprocess.run(plot_command, check=True)

if __name__ == "__main__":
    main()