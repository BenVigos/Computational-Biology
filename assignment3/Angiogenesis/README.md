# Paper-aligned 2D angiogenesis model

This model lives in `assignment3/Angiogenesis`.
It adapts the reference implementation in `assignment3/Paper files/Simulation` so the **cell mechanics, thresholds, growth laws, and field names follow the paper model**, while the lattice remains **2D** (`z = 1`).

## Staged presets: simple to complex
The model is modular again.
You can now choose staged presets in `Simulation/AngiogenesisConfig.py` by changing:
- `SELECTED_PRESET = "paper_mvp"`
- `SELECTED_PRESET = "paper_hypoxia"`
- `SELECTED_PRESET = "paper_sprouting"`
- `SELECTED_PRESET = "tumor_oxygen_only"`
- `SELECTED_PRESET = "paper_full"`

What each preset enables:
- `paper_mvp`
  - boundary vessels + tumor mass
  - no sprouts
  - no phenotype switching
  - no growth
  - no mitosis
  - lighter monitoring
- `paper_hypoxia`
  - boundary vessels + tumor mass
  - hypoxic / necrotic switching
  - tumor growth
  - no vascular sprouting growth
  - no mitosis
- `paper_sprouting`
  - adds initial sprouts
  - enables vascular VEGF-driven growth
  - keeps mitosis off
- `tumor_oxygen_only`
  - keeps the boundary vessels as oxygen sources
  - no sprouts
  - no vascular growth or vascular switching
  - no mitosis
  - tumor growth + oxygen diffusion focus
- `paper_full`
  - full paper-aligned 2D run
  - boundary vessels, sprouts, switching, tumor growth, vascular growth, and mitosis all on

## Main toggles
You can also build your own custom stage by editing booleans in `Simulation/AngiogenesisConfig.py`:
- `enable_initial_vascular_strip`
- `enable_initial_sprouts`
- `enable_initial_tumor_mass`
- `enable_type_switching`
- `enable_tumor_growth`
- `enable_vascular_growth`
- `enable_mitosis`
- `enable_python_monitoring`
- `enable_timescale_separation`

## Timescale separation (diffusion vs growth)
The model supports decoupled timescales so diffusion can relax before cell growth updates.

Parameters in `Simulation/AngiogenesisConfig.py`:
- `enable_timescale_separation = True`
- `diffusion_relaxation_mcs = 100`
- `growth_window_mcs = 10`
- `timescale_cycle_offset_mcs = 0`

Behavior:
- for `diffusion_relaxation_mcs`, only diffusion/fields evolve,
- for `growth_window_mcs`, growth and mitosis rules execute,
- then the cycle repeats.

Monitoring now includes:
- `phase` (`diffusion` or `growth`)
- `is_growth_window` (`0` or `1`)

## Grid size
The model is currently configured for a smaller development domain:
- Python config: `100 x 100`
- CC3D XML domain: `100 x 100 x 1`
- steps: `2000`

This makes iteration faster while preserving the paper-style behavior.

## What the simulation does
At startup the model can create:
- `Vascular` boundary vessels on **all four sides** of the domain,
- several `InactiveNeovascular` sprouts attached near the left boundary vessel,
- a central tumor mass seeded as `Normal` cells.

During the run, depending on the preset:
- oxygen diffuses from `Vascular` and `ActiveNeovascular` cells,
- `Normal` tumor cells can become `Hypoxic` when oxygen drops below the paper nutrient threshold,
- tumor cells can become `Necrotic` when oxygen drops below the paper necrotic threshold,
- `Hypoxic` cells secrete `VEGF2`,
- vascular-like cells secrete `VEGF1`,
- `InactiveNeovascular` and `ActiveNeovascular` cells grow using the same VEGF-driven rules as the paper steppables,
- tumor and vascular cells can divide at the same paper doubling volumes,
- Python monitoring records growth and field metrics to CSV.

## Files and where parameters are set
### `Simulation/AngiogenesisConfig.py`
This is the main Python-side parameter file.
It contains:
- staged presets via `PRESETS`
- the active selection via `SELECTED_PRESET`
- initialization toggles
- relative placement controls (fractions of lattice dimensions)
- paper-derived thresholds and mechanics
- monitoring toggles

Important paper-derived parameters here include:
- `nutrient_thresh = 5.0`
- `necrotic_thresh = 1.0`
- `tumor_growth_start_mcs = 100`
- `tumor_target_volume = 33.0`
- `tumor_lambda_volume = 10.0`
- `tumor_target_surface = 90.0`
- `tumor_lambda_surface = 2.0`
- `vascular_target_volume = 60.0`
- `vascular_lambda_volume = 13.0`
- `vascular_target_surface = 150.0`
- `vascular_lambda_surface = 3.0`
- `tumor_doubling_volume = 54.0`
- `vascular_doubling_volume = 80.0`

### `Simulation/Angiogenesis.xml`
This file contains the CompuCell3D XML configuration.
It keeps the paper-style:
- cell types,
- contact energies,
- `Volume` and `Surface` plugins,
- `VEGF1` / `VEGF2` fields,
- diffusion, uptake, secretion, and constant-concentration field sources.

It now uses a smaller development lattice:
- `<Dimensions x="100" y="100" z="1"/>`
- `<Steps>2000</Steps>`

### `Simulation/AngiogenesisSteppables.py`
This file contains the runtime rules.
Key staged-control regions:
- `ConstraintInitializerSteppable.start`
  - respects `enable_initial_vascular_strip`, `enable_initial_sprouts`, and `enable_initial_tumor_mass`
- `GrowthSteppable.step`
  - respects `enable_type_switching`, `enable_tumor_growth`, and `enable_vascular_growth`
- `MitosisSteppable.step`
  - respects `enable_mitosis`
- `MonitoringSteppable`
  - keeps Python-side metrics and CSV export available across presets

## Monitoring and CSV output
Python-side monitoring remains enabled by default.
The CSV is written to:
- `assignment3/Angiogenesis/results/angiogenesis_metrics.csv`

Tracked metrics include:
- `normal_cells`
- `hypoxic_cells`
- `necrotic_cells`
- `tumor_like_cells`
- `active_neovascular_cells`
- `vascular_cells`
- `inactive_neovascular_cells`
- `avg_tumor_volume`
- `avg_tumor_target_volume`
- `avg_tumor_volume_growth_rate`
- `total_tumor_volume_growth_rate`
- `mean_tumor_oxygen`
- `mean_tumor_vegf1`
- `mean_tumor_vegf2`
- vascular / endothelial-compatible aliases for plotting

Monitoring toggles are in `Simulation/AngiogenesisConfig.py`:
- `enable_python_monitoring`
- `monitor_frequency`
- `monitor_to_console`
- `monitor_to_csv`
- `monitor_include_growth_rates`
- `monitor_include_field_means`
- `monitor_include_vascular_metrics`

## Plotting the monitoring CSV
Use:
- `assignment3/scripts/plot_angiogenesis_metrics.py`

Example:
- `python assignment3/scripts/plot_angiogenesis_metrics.py --input assignment3/Angiogenesis/results/angiogenesis_metrics.csv`

Plots are written by default to:
- `assignment3/Angiogenesis/results/plots/`

## 2D adaptation notes
This is **not** a direct copy of the paper's 3D project.
The following choices were made intentionally:
- the lattice remains 2D,
- initialization is generated in Python instead of loading the paper PIF,
- the model keeps the paper cell mechanics and thresholds, but applies them on a 2D domain,
- the model can now be enabled in steps again, from simple to more complex.

## Quick validation checklist
When the model is working as expected, you should see:
1. oxygen highest near the boundary vessels,
2. low-oxygen tumor regions become `Hypoxic` in the hypoxia/sprouting/full presets,
3. the most deprived cells turn `Necrotic`,
4. `VEGF2` rises around hypoxic regions,
5. `InactiveNeovascular` / `ActiveNeovascular` cells enlarge and divide in the fuller presets.
