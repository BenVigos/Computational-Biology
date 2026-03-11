# Paper-aligned 2D angiogenesis model

This model lives in `assignment3/Angiogenesis`.
It adapts the reference implementation in `assignment3/Paper files/Simulation` so the **cell mechanics, thresholds, growth laws, and field names follow the paper model**, while the lattice remains **2D** (`z = 1`).

## What was changed
The old MVP-style model used the simplified cell set
`Tumor / Hypoxic / Endothelial / BloodVessel`.

The current model now follows the paper cell taxonomy:
- `Normal`
- `Hypoxic`
- `Necrotic`
- `ActiveNeovascular`
- `Vascular`
- `InactiveNeovascular`

It also now uses the paper field layout:
- `Oxygen`
- `VEGF1`
- `VEGF2`

## What the simulation does
At startup the model creates:
- a vertical `Vascular` strip on the left side,
- several `InactiveNeovascular` sprouts attached to that vessel,
- a central tumor mass seeded as `Normal` cells.

During the run:
- oxygen diffuses from `Vascular` and `ActiveNeovascular` cells,
- `Normal` tumor cells become `Hypoxic` when oxygen drops below the paper nutrient threshold,
- tumor cells become `Necrotic` when oxygen drops below the paper necrotic threshold,
- `Hypoxic` cells secrete `VEGF2`,
- vascular-like cells secrete `VEGF1`,
- `InactiveNeovascular` and `ActiveNeovascular` cells grow using the same VEGF-driven rules as the paper steppables,
- tumor and vascular cells divide at the same paper doubling volumes,
- Python monitoring records growth and field metrics to CSV.

## Files and where parameters are set
### `Simulation/AngiogenesisConfig.py`
This is the main Python-side parameter file.
It contains the paper-derived values that were previously hard-coded in the paper steppables:

- phenotype thresholds:
  - `nutrient_thresh = 5.0`
  - `necrotic_thresh = 1.0`
- delay before phenotype switching / growth:
  - `tumor_growth_start_mcs = 100`
- tumor mechanics:
  - `tumor_target_volume = 33.0`
  - `tumor_lambda_volume = 10.0`
  - `tumor_target_surface = 90.0`
  - `tumor_lambda_surface = 2.0`
- vascular mechanics:
  - `vascular_target_volume = 60.0`
  - `vascular_lambda_volume = 13.0`
  - `vascular_target_surface = 150.0`
  - `vascular_lambda_surface = 3.0`
- tumor growth law:
  - `0.04 * oxygen / (10 + oxygen)` for volume
  - `0.12 * oxygen / (10 + oxygen)` for surface
- vascular growth law:
  - `0.06 * VEGF2 / (0.5 + VEGF2)` for volume
  - `0.15 * VEGF2 / (0.5 + VEGF2)` for surface
- vascular activation / crowding thresholds:
  - `vascular_vegf_activation_threshold = 0.5`
  - `inactive_neighbor_area_limit = 70.0`
  - `active_neighbor_area_limit = 50.0`
- mitosis thresholds:
  - `tumor_doubling_volume = 54.0`
  - `vascular_doubling_volume = 80.0`

This file also contains:
- the 2D geometry used for initialization,
- monitoring toggles,
- runtime feature toggles like `enable_tumor_growth`, `enable_vascular_growth`, and `enable_mitosis`.

### `Simulation/Angiogenesis.xml`
This file contains the CompuCell3D XML configuration.
It now mirrors the paper model more closely by using:
- the paper cell types,
- the paper `Contact` energies,
- `Volume` and `Surface` plugins,
- `VEGF1` and `VEGF2` chemotaxis entries,
- paper-inspired diffusion / secretion / uptake settings.

Important XML regions:
- `<Plugin Name="CellType">`: defines all simulation cell types.
- `<Plugin Name="Contact">`: paper contact energies.
- `<Plugin Name="Chemotaxis">`: paper chemotaxis couplings.
- `<Steppable Type="DiffusionSolverFE">`: field diffusion, secretion, uptake, and constant-concentration sources.

### `Simulation/AngiogenesisSteppables.py`
This file contains the runtime rules.
Key regions:
- `ConstraintInitializerSteppable.start`:
  - creates the 2D vessel, sprouts, and tumor geometry.
- `GrowthSteppable.start`:
  - assigns paper target volume / surface parameters to each cell type.
- `GrowthSteppable.step`:
  - implements the paper growth and phenotype-switching logic.
- `MitosisSteppable.step`:
  - divides `Normal`, `Hypoxic`, `ActiveNeovascular`, and `InactiveNeovascular` cells at the paper thresholds.
- `MitosisSteppable.update_attributes`:
  - resets parent/child mechanics after division using the paper values.
- `MonitoringSteppable`:
  - samples tumor, necrotic, and vascular metrics in Python and optionally writes them to CSV.

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
- the lattice remains `180 x 180 x 1`,
- initialization is generated in Python instead of loading the paper PIF,
- the model keeps the paper cell mechanics and thresholds, but applies them on a 2D domain,
- the README and monitoring were expanded so parameter changes are easy to inspect from Python.

## Quick validation checklist
When the model is working as expected, you should see:
1. oxygen highest near the left-side vascular strip,
2. low-oxygen tumor regions become `Hypoxic`,
3. the most deprived cells turn `Necrotic`,
4. `VEGF2` rises around hypoxic regions,
5. `InactiveNeovascular` / `ActiveNeovascular` cells enlarge and divide when exposed to sufficient VEGF.
