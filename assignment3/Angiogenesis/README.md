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

## How to toggle features
Open `Simulation/AngiogenesisConfig.py` and change:
- `SELECTED_PRESET = "mvp"` for the simplest validation model.
- `SELECTED_PRESET = "hypoxia_hif"` to add HIF-1α dynamics.
- `SELECTED_PRESET = "sprouting"` to add endothelial growth + chemotaxis.
- `SELECTED_PRESET = "full"` to also enable mitosis and apoptosis.

You can also edit any individual boolean in `ModelConfig` if you want a custom combination.

## Suggested validation order
1. `mvp`: check oxygen stays highest near the vessel and hypoxia appears in the tumor core.
2. `hypoxia_hif`: verify `mean_HIF` increases in the log when oxygen is low.
3. `sprouting`: check the endothelial cells respond to VEGF gradients.
4. `full`: check division and low-oxygen shrinkage are stable.

## Notes
- Diffusion/decay coefficients are defined in `Simulation/Angiogenesis.xml` because CompuCell3D reads them from XML.
- Most reaction, switching, growth, and secretion logic is in `Simulation/AngiogenesisSteppables.py` so it can be toggled easily.
- The chemotaxis toggle is implemented in Python, but if the local CompuCell3D build exposes a different chemotaxis API, the rest of the model will still run and chemotaxis will stay off.

