# Toggleable CompuCell3D angiogenesis MVP

This model lives in `assignment3/Angiogenesis` and is designed to be built up in stages from the equations in `assignment3/equations.md`.

## What the MVP includes
- GGH / CPM mechanics through the XML `Volume` and `Contact` terms.
- An `Oxygen` diffusion field.
- A `VEGF` diffusion field.
- A frozen `BloodVessel` oxygen source.
- A tumor mass that can become `Hypoxic` when oxygen drops below a threshold.
- Hypoxic tumor cells secreting VEGF.
- Optional HIF-1α accumulation.
- Optional endothelial growth, chemotaxis, apoptosis, and mitosis.
- Optional Python-side monitoring of model metrics over time.

## How to toggle features
Open `Simulation/AngiogenesisConfig.py` and change:
- `SELECTED_PRESET = "mvp"` for the simplest validation model.
- `SELECTED_PRESET = "hypoxia_hif"` to add HIF-1α dynamics.
- `SELECTED_PRESET = "sprouting"` to add endothelial growth + chemotaxis.
- `SELECTED_PRESET = "full"` to also enable mitosis and apoptosis.

You can also edit any individual boolean in `ModelConfig` if you want a custom combination.

## Python monitoring toggles
The `MonitoringSteppable` computes time-series metrics directly in Python and can
store them in memory, print them, and/or write them to CSV.

Key toggles in `Simulation/AngiogenesisConfig.py`:
- `enable_python_monitoring`: master switch for Python monitoring.
- `monitor_frequency`: sample metrics every N MCS.
- `monitor_to_console`: print a compact monitoring line to the console.
- `monitor_to_csv`: write sampled metrics to CSV.
- `monitor_include_growth_rates`: include volume/target-volume growth rates.
- `monitor_include_field_means`: include mean oxygen and VEGF near cells.
- `monitor_include_endothelial_metrics`: include endothelial counts/volumes.
- `monitor_output_dir`: output folder, relative to `assignment3/Angiogenesis`.
- `monitor_output_filename`: CSV filename.

Current sampled metrics include:
- tumor / hypoxic / endothelial counts
- hypoxic fraction
- average and total tumor volume
- average tumor target volume
- mean HIF
- average tumor volume growth rate (`cell.volume` change per MCS)
- average tumor target-volume growth rate (`cell.targetVolume` change per MCS)
- total tumor volume growth rate
- optional endothelial volume metrics
- optional mean oxygen / VEGF values around tumor and endothelial cells

When CSV output is enabled, the file is written to:
- `assignment3/Angiogenesis/results/angiogenesis_metrics.csv`

## Plotting the monitoring CSV
A helper script is available at `assignment3/scripts/plot_angiogenesis_metrics.py`.
It reads a monitoring CSV such as `results/angiogenesis_metrics.csv` and creates
summary plots automatically.

Example:
- `python assignment3/scripts/plot_angiogenesis_metrics.py --input assignment3/Angiogenesis/results/angiogenesis_metrics.csv`

By default, figures are written to:
- `assignment3/Angiogenesis/results/plots/`

Generated plot groups include:
- cell counts and hypoxic fraction
- tumor volume / target volume / HIF
- growth rates
- tumor and endothelial field means
- endothelial metrics

Useful options:
- `--output-dir <dir>` to choose a different output folder
- `--prefix run1` to change output filenames
- `--format png|pdf|svg`
- `--show` to display figures interactively

## Suggested validation order
1. `mvp`: check oxygen stays highest near the vessel and hypoxia appears in the tumor core.
2. `hypoxia_hif`: verify `mean_HIF` increases in the log when oxygen is low.
3. `sprouting`: check the endothelial cells respond to VEGF gradients.
4. `full`: check division and low-oxygen shrinkage are stable.

## Notes
- Diffusion/decay coefficients are defined in `Simulation/Angiogenesis.xml` because CompuCell3D reads them from XML.
- Most reaction, switching, growth, and secretion logic is in `Simulation/AngiogenesisSteppables.py` so it can be toggled easily.
- The chemotaxis toggle is implemented in Python, but if the local CompuCell3D build exposes a different chemotaxis API, the rest of the model will still run and chemotaxis will stay off.
- The `MonitoringSteppable` also keeps sampled rows in Python as `self.metric_history` during the run.
