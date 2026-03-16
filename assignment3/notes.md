# Notes


## Model development and equations
We should look at:

### 1. Tumor dynamics
* tumor cell count: growth is slow, then accelerates, eventually limited by oxygen (MM dynamics)
* hypoxic fraction: tumor grows, less oxygne: hypoxic cells increase, secrete VEGF2: new vessels appear, oxygen improves: hypoxia decrease
* necrotic fraction
* tumor volume

### 2. Vascular response
* endothelial cell count: VEGF chemotaxis activates sprouting, hence increase population.
* active sprout cells
* vascular volume

### 3. Signaling fields
* oxygen: tumor consumes oxygen, oxygen drops, angiogenesis restores supply
* VEGF2: rise with hypoxia

### 4. System-level dynamics
* tumor growth rate: faster when oxygen is available, active angiogenesis. Then limited
* oxygen recovery after angiogenesis

## Plots Results and Comparison

### Introduction
- Will go over figures and videos in this assignment
- Discuss why they were chosen and what they show
- Highlight key features from the results

### Conceptual overview
- Understanding what is being modelled is crucial to interpreting the results
- Figure 1: Schematic of the tumor-angiogenesis model
  - Shows the key components and interactions in the model
  - Tumor cells, hypoxia, VEGF, endothelial cells
  - Helps to understand the biological processes being simulated

### Results
- Methodology:
    - Isolate each of the processes present (tumor growth, angiogenesis and impact of hif1a) 
    - helps understand their individual contributions
    - That's why we have 3 different figures

- Figure 1: Comparison of tumor growth in simulations without angiogenesis and in the full model including angiogenesis

    - Panel (a): Temporal evolution of tumor size
      - assumes spherical tumor geometry: effective tumor radius
      - mean radial distance of tumor cells from the tumor centroid allows us to understand how the tumor changes spatially over time
      - Deviation in the two measures can indicate tumor is not fully spherical or doesn't have uniform cell density
      
    - Panel (b): Tumor cell composition over time
      - Provides insights into number of normoxic, hypoxic, and necrotic cells
      - Shows how angiogenesis affects the balance between proliferative and oxygen-limited states

    - In bowth panels we see that angiogenesis allows for more sustained tumor growth.

    - Panel (c): Hypoxic fraction
      - Answers the question: does angiogenesis reduce the degree of oxygen stress experienced by the tumor during progression?
      - We see that angiogenesis reduces the hypoxic fraction, indicating improved oxygen supply to the tumor

- Figure 2: Comparison of endothelial activation given a fully hypoxic non-growing tumor with and without HIF-1$\alpha$ signaling enabled

    - Panel (a): Time evolution of VEGF and HIF-1$\alpha$-associated signaling
      - Shows how the inclusion of HIF-1$\alpha$ alters the molecular pro-angiogenic response (i.e., VEGF production)
      - We see no significant change in the production of VEGF.
    
    - Panel (b): Numbers of active and inactive neovascular cells
      - Provides a direct measure of how signaling differences translate into endothelial activation and sprout dynamics
      - We see that with HIF-1$\alpha$, there are more active neovascular cells, indicating enhanced angiogenesis
    
    - Panel (c): Tumor oxygen levels together with the minimum distance between the tumor and nearby sprouts
      - Links biochemical signaling to the spatial approach of the chemotaxis of the vascular network toward the tumor mass
      - We see that with HIF-1$\alpha$, tumor oxygen levels improve over time as sprouts approach, demonstrating the functional impact of enhanced angiogenesis

- Figure 3: Integrated comparison of tumor-only growth and the full tumor--angiogenesis model
    - Goal: to directly compare the tumor dynamics and signaling environment in the presence versus absence of angiogenesis, highlighting how vascular support reshapes tumor progression
    - Contains aggregate metrics from both previous figures to provide a comprehensive view of the differences between the two conditions

    - Panel (a): Temporal evolution of total tumor-like cells together with the normoxic, hypoxic, and necrotic subpopulations
      - Allows direct assessment of how angiogenesis reshapes tumor composition
      - We see that angiogenesis leads to a larger total tumor cell population, with a higher proportion of normoxic cells.
      
    - Panel (b): VEGF and HIF-1$\alpha$ signaling dynamics
      - Allows us to directly couple the signaling fields to the quantitative makeup of the tumor population

    - Panel (c): Tumor oxygen levels together with the minimum tumor-to-vessel distance
      - Demonstrates how vascular proximity is associated with changes in oxygen availability during tumor progression
      - We see that in the angiogenesis model, tumor oxygen levels improve over time as the minimum distance to vessels decreases, confirming that angiogenesis effectively enhances oxygen delivery to the tumor
      - In contrast, in the non-angiogenesis model, tumor oxygen levels remain low and the minimum distance to vessels does not change (until very late), illustrating the lack of vascular support
      - The late part of the figures can also give us insights into some of the finite-size effects of the model.

## Video descriptions
- Becaue an image is worth a thousand words, a video is worth a million. The videos provide a dynamic visualization of the processes we are simulating, allowing us to see how the tumor and vasculature evolve over time in a way that static figures cannot capture.
- Video 1: Tumor growth without angiogenesis
- Video 2: Angiogenesis without tumor growth
- Video 3: Tumor growth with angiogenesis





## Old stuff
- General model: tumour that begins normoxic, outgrows its oxygen supply,
enters a hypoxic crisis that drives VEGF2 secretion, recruits new vessels via
angiogenesis, and recovers.

### Plot 1 - tumor dynamics

1. Cell counts

***Metrics***
```
normal_cells   = count of cells with type NORMAL
hypoxic_cells  = count of cells with type HYPOXIC
necrotic_cells = count of cells with type NECROTIC
tumor_like_cells = sum of all three
```

Cell count is the most direct readout of tumour growth. The three subtypes tells us
whether growth is healthy (normal cells dominating), stressed (hypoxic cells
rising), or dying (necrotic cells appearing). The ratio between them at any moment
summarises the tumour's oxygen status. (KEEP)

2. Tumor phenotype fractions and VGEF2

```
hypoxic_fraction  = hypoxic_cells  / tumor_like_cells
necrotic_fraction = necrotic_cells / tumor_like_cells
normal fraction = 1 − hypoxic_fraction − necrotic_fraction
mean_tumor_vegf2  = mean of VEGF2 field sampled at each tumour cell's centre of mass
```

Causal loop of our model.  VEGF2 on the right axis shows the signalling coupling: hypoxia
drives VEGF2 secretion, which in turn drives angiogenesis. (KEEP)

3. Per-phenotype cell volumes

```
avg_normal_volume   = mean(cell.volume for cell in normal_cells)
avg_hypoxic_volume  = mean(cell.volume for cell in hypoxic_cells)
avg_necrotic_volume = mean(cell.volume for cell in necrotic_cells)
```

In a CPM, cells grow by increasing `targetVolume` and the Metropolis dynamics
drive actual volume toward the target. Average volume tells you whether the
growth mechanics are working correctly for each phenotype. Necrotic cells should
shrink toward zero (in our case: hollow-core); normal cells should oscillate around the target (sawtooth path is because of mitosis); hypoxic
cells should be slightly below target (reduced growth rate).

4. Growth pressure and HIF1a

```
avg_tumor_target_volume = mean(cell.targetVolume for all tumour cells)
avg_tumor_volume        = mean(cell.volume for all tumour cells)
growth_pressure         = target - actual  (shaded region)
mean_tumor_hif1a        = mean(cell.dict["hif1a"] for all tumour cells)
```

The gap between target and actual volume is the mechanical drive for growth in the
CPM: a cell with targetVolume > volume exerts outward pressure on neighbours. During the hypoxic crisis: gas are created because cells that were normally growing have their volume reset or gorwth rate reduced (necrotic / hypoxic)
HIF1a peak is expected before growth pressure peaks, because HIF drives VGEF which drives angiogenesis that restores 02 and allows growth.

### Plot 2 - Endothelial Response

1. Neovascular populations
```
active_neovascular_cells   = count of ACTIVENEOVASCULAR cells
inactive_neovascular_cells = count of INACTIVENEOVASCULAR cells
```
Active = tip cells (chemotaxing toward VEGF2, growing).
Inactive = stalk cells (lower growth rate, following the tip).


2. Per-type neovascular volumes

```
avg_active_neovascular_volume   = mean(cell.volume for active cells)
avg_inactive_neovascular_volume = mean(cell.volume for inactive cells)
```

Active tip cells should be growing (volume > target) because they are receiving
VEGF2 above threshold and have `_vascular_growth_step` adding to `targetVolume`.
Inactive stalk cells should stay near `vascular_target_volume = 40` if no VEGF2
is driving them. The divergence between the two lines should validate the differential
growth rule.

3. Endothelial Growth Rate

```
avg_endothelial_volume_growth_rate = _mean_delta_rate(
    current_endothelial_volumes, previous_endothelial_volumes, delta_mcs
)
= (sum(current) - sum(previous)) / (delta_mcs * n_current_cells)
```

This is the per-cell average rate of volume change for neovascular cells. A
positive rate means cells are growing on average; a negative rate means they
are shrinking. The early period reveals initialisation artefacts; the main
period shows how actively sprouts are growing in response to VEGF2.

### Plot 3 - Signaling Fields

1. Oxygen Field
```
mean_tumor_oxygen = mean(Oxygen[x,y,z] at each tumour cell's COM)
```

2. VGEF field
```
mean_tumor_vegf2 = mean(VEGF2[x,y,z] at each tumour cell's COM)
mean_tumor_vegf1 = mean(VEGF1[x,y,z] at each tumour cell's COM)
```
VEGF2 is secreted by hypoxic cells.
VEGF1 is secreted by vascular/neovascular cells and drives
a secondary chemotaxis signal. (DO WE NEED IT?)

MAYBE COMBINE TWO FIELDS

### Plot 4 - System Level Dynamics

1. Tumor volume growth rate
```
total_tumor_volume_growth_rate = (
    sum(current cell volumes) - sum(previous cell volumes)
) / delta_mcs
```
This is the most direct readout of whether the tumour is expanding, stagnant, or
shrinking. Because it uses total volume sums, it correctly captures growth from
both cell enlargement and cell division. The zero-crossings mark biological phase
transitions: onset of hypoxic stasis, angiogenic recovery, and re-acceleration.
Target growth rate (green dashed) near zero throughout: confirms this column is
tracking target-volume changes, which are minimal because necrotic cells reset
to zero and living cells reset after division

### Plot 5 - HIF1A gene network
1. HIF1a vs 02
```
mean_tumor_hif1a = mean(cell.dict["hif1a"] for all tumour cells)
```
2.  VEGF drive and VEF2 field
```
mean_tumor_vegf_drive = mean(cell.dict["vegf_drive"] for all tumour cells)
vegf_drive = V_basal + (V_max - V_basal) * f(H)
f(H)       = H^n / (K_H^n + H^n)    # Hill function
```

How HIF1 boosts VEGF.

### Plot 6 - Spatial Dynamics
1. Tumor Morphology
```
tumor_effective_radius = sqrt(total_tumor_volume / π)
# treats the tumour as a circle of the same area

tumor_mean_radius = mean(euclidean_distance(cell_COM, tumor_centroid)
                         for all tumour cells)
# average distance of cells from the weighted centre of mass

tumor_radial_cv = std(radial_distances) / mean(radial_distances)
# coefficient of variation: 0 = all cells at same radius, high = irregular

tumor_aspect_ratio = x_span / y_span
# bounding box ratio: 1.0 = square, <1 = taller than wide
```

Effective radius tells you how large the tumour is as a pure size metri: tumor stops expanding when oxygen is limited.
Mean radius tells you where cells actually are relative to the centre, if
mean radius approaches effective radius (ratio ~1.0 instead of the expected
~0.67 for a filled disk), cells are concentrated at the periphery and the
interior is becoming hollow. Radial CV measures boundary roughness and
invasiveness: irregularity. Aspect ratio measure absence of directional growth bias.

2. Vascular Proximity

```
# For each tumour cell, find the nearest neovascular or parent vessel cell
dists_to_any_vessel = [min_distance(tumour_cell, all_vessel_cells)
                       for tumour_cell in tumour_cells]
mean_dist_to_vessel = mean(dists_to_any_vessel)
min_dist_to_vessel  = min(dists_to_any_vessel)

# Separately for neovascular sprouts only:
mean_dist_to_sprout = mean(min_distance(tumour_cell, neovascular_cells)
                           for tumour_cell in tumour_cells)
```

This is the spatial validation of angiogenesis. Falling distances confirm that
vessels are genuinely approaching the tumour in space, not just appearing in
counts. The gap between mean-to-any-vessel and mean-to-sprout shows how much
of the proximity improvement comes from new sprouts vs inherited geometry from
parent walls.

3. Radial Architecture
```
normal_mean_radius   = mean(distance(normal_cell,   tumour_centroid))
hypoxic_mean_radius  = mean(distance(hypoxic_cell,  tumour_centroid))
necrotic_mean_radius = mean(distance(necrotic_cell, tumour_centroid))
tumor_effective_radius  # repeated for reference

hypoxic_core_offset = distance(tumour_centroid, centroid_of(hypoxic + necrotic cells))

hypoxic_inner_fraction  = fraction of hypoxic cells within 0.45 * effective_radius
necrotic_inner_fraction = fraction of necrotic cells within 0.45 * effective_radius
```

In a biologically correct tumour, oxygen diffuses inward from peripheral vessels,
creating a gradient: normoxic cells at the outside, hypoxic in the middle,
necrotic at the core. This predicts: `necrotic_radius < hypoxic_radius < normal_radius`.
Measuring whether this ordering holds is the spatial validation of your O₂-phenotype coupling.
Core offset shows whether the tumour is growing symmetrically or has a directional bias.


Compared to model without angiogenesis... both tumor experience identical hypoxic crises driven by the same O₂ depletion and HIF1a accumulation. Angiogenesis is the only difference, and its presence produces nearly 2× more oxygen delivery (crisis is not prevented by angiogenesis but only vascularised condition can escape it), 3× more tumour cells, and allows the HIF1a network to downregulate as intended. Without angiogenesis, the GRN remains locked in a hypoxic steady state: HIF1alpha elevated, VEGF2 high, growth impaired — exactly the regulatory phenotype the model predicts for an unrescued tumour. This directly demonstrates that GRN is not just running in parallel with the simulation but is coupled to the biological outcome in a meaningful way.


*Anna's ideas*
- compare full config with tumor growth only (but hif1a)
- keep all panels from plot 1, since relvant to explain the full mechanism
- discard plot 2, or keep only the first panel: mitosis always produces active daughters, so nearly all neovascular cells become and remain active (shouldnt be the case)
- merge VEGF and Oxygen field in plot 3, noting that VEGF1 is not relevant
- discard plot 4
- keep plot 5, interesting to compare to the one without angiogenesis. It shows HIF1a dynamics working
- keep some spatial metrics from 6, especially those that show that are tumor has hollow structure and that necrotic_radius < hypoxic_radius < normal_radius for necrotic "core". And vessels distance
- limitation: no decoupling in metrics collection?, boundaries effects? FiPy vs CC3D solver? pharmacokinetics and dynamics extensions?

### Results and Discussion

