# Assignment 3 – Detailed Project Plan (Ideas 1 & 2)

> **Course context / constraints** (from `assignment3/assignment.md`)
>
> - Must involve **coupled PDEs** based **reaction–diffusion** system(s) *and* **cell-based model(s)**.
> - Should combine course concepts (enzyme kinetics, metabolic networks, biotech, GRNs, nonlinear dynamics, evolutionary dynamics, reaction–diffusion, cell-based modelling).
>
> **Time budget:** ~1.5 weeks, team of 3.

---

## How to read this document
We propose **two feasible neuroscience-flavoured modeling directions**. Both satisfy the “coupled PDEs + cell-based” requirement, but they emphasize different phenomena:

- **Idea 1 (Neuromodulator learning wave):** a *chemical* diffuses (reaction–diffusion PDE). It gates plasticity in discrete neuron/synapse agents.
- **Idea 2 (Excitable tissue waves):** the *activity itself* is modeled as an excitable-medium PDE (FitzHugh–Nagumo/Barkley). Optional extra PDEs for metabolism or modulators provide “coupling”. Discrete agents can represent specialized cells (e.g., glia, neurons) that locally source/sink fields.

We recommend implementing **one** as the main deliverable (depending on what the group prefers), and keeping the other as a contingency/extension.

---

# IDEA 1 — Reaction–diffusion neuromodulator + synaptic plasticity (“learning wave”)

## 1. Motivation
Brains don’t just transmit signals; they **adapt** their connectivity based on experience. Many learning processes are modulated by spatially distributed chemicals such as **dopamine** or **nitric oxide** that can act beyond a single synapse.

This idea aims to show a simple, visual story:

1. A stimulus activates a small group of neurons.
2. A **neuromodulator pulse** is released (e.g., “reward”).
3. The neuromodulator **diffuses**, creating a wave/gradient.
4. Synapses exposed to sufficient neuromodulator undergo potentiation.
5. After repeated pairing, the network preferentially routes activity along the “learned” path.

This is attractive for the assignment because:

- The neuromodulator is naturally a **reaction–diffusion PDE**.
- Neurons/synapses are naturally **discrete agents**.
- Enzyme kinetics fits cleanly via **Michaelis–Menten clearance** of neuromodulator.
- GRN modelling fits as slow adaptation (e.g., receptor expression) with Hill functions.

## 2. Goal (what we will demonstrate)
A minimal set of deliverables:

- **Coupled PDE(s):** 2D reaction–diffusion for neuromodulator concentration (optionally add a second PDE such as “degrading enzyme” or “metabolic resource”).
- **Cell-based component:** discrete neurons on a 2D lattice or embedded graph; synapses with plastic weights.
- **Learning outcome:** show that a reward pulse **selectively strengthens** synapses in the region it reaches, producing a learned route from an input zone to an output zone.

### Concrete “demo scenario”
- A 2D grid domain with two neuron clusters: **Input** (left) and **Output** (right).
- Baseline: random sparse connections / or local nearest-neighbour coupling.
- Training: repeatedly stimulate Input; release neuromodulator near Output at a short delay (mimicking reward).
- Result: after training, stimulating Input produces stronger propagation to Output.

## 3. Theoretical background

### 3.1 Reaction–diffusion PDE for neuromodulator
Let \(M(x,y,t)\) be extracellular neuromodulator concentration (e.g., dopamine proxy).

$$
\frac{\partial M}{\partial t} = D_M \nabla^2 M + S_M(x,y,t) - \underbrace{k_M M}_{\text{linear decay}} - \underbrace{\frac{V_{\max} M}{K_M + M}}_{\text{enzymatic clearance (Michaelis–Menten)}}
$$

- \(D_M\): diffusion coefficient.
- \(S_M\): source term. Nonzero near “reward sites” or active neurons.
- \(k_M\): first-order loss (uptake/decay).
- \(V_{\max}, K_M\): enzymatic clearance parameters (enzyme kinetics / biotech lever).

**Boundary conditions:** no-flux (Neumann) boundaries are realistic for a local tissue patch:
$$
\nabla M \cdot \vec{n} = 0 \quad \text{on the boundary}
$$

### 3.2 Cell-based neurons (rate model) and synapses
Represent neurons as discrete agents \(i = 1,\dots,N\). Each has an activity \(a_i(t)\) (firing rate proxy).

A simple, stable rate dynamics:
$$
\tau_a \frac{d a_i}{dt} = -a_i + \phi\left( \sum_j w_{ij} a_j + I_i(t) \right)
$$

- \(w_{ij}\): synaptic weight from neuron \(j\) to \(i\).
- \(I_i\): external input (stimulus pulses).
- \(\phi(\cdot)\): nonlinearity (e.g., sigmoid or ReLU clipped), giving nonlinear dynamics.

### 3.3 Plasticity rule gated by neuromodulator (learning)
Let neuron \(i\) occupy location \((x_i,y_i)\). The neuron “senses” local modulator \(M_i(t) = M(x_i,y_i,t)\).

A minimal 3-factor Hebbian rule:
$$
\frac{d w_{ij}}{dt} = \eta\, g(M_i)\, a_i a_j - \lambda_w w_{ij}
$$

- \(\eta\): learning rate.
- \(\lambda_w\): weight decay / regularization.
- \(g(M)\): gating function. A Hill function fits the GRN-style nonlinear thresholding:
$$
 g(M) = \frac{M^n}{K_g^n + M^n}
$$

This yields:
- If there is no neuromodulator, \(g(M)\approx 0\) and weights barely change.
- If neuromodulator is high, \(g(M)\approx 1\) and correlated activity strengthens weights.

### 3.4 Optional GRN (slow receptor expression)
To explicitly include “gene regulatory network modelling”: define a slow variable \(R_i(t)\) representing receptor / plasticity protein expression.

$$
\frac{dR_i}{dt} = \alpha_R \frac{M_i^n}{K_R^n + M_i^n} - \gamma_R R_i
$$

Then modulate plasticity by \(R_i\):
$$
\frac{dw_{ij}}{dt} = \eta\, R_i\, a_i a_j - \lambda_w w_{ij}
$$

This creates a biologically plausible separation of timescales: fast modulator diffusion; slower gene/protein response; slower synaptic change.

### 3.5 Coupling summary (what depends on what)
- PDE \(M(x,y,t)\) depends on sources \(S_M\) which depend on agents (active neurons or reward location).
- Agents read back \(M_i\) from the PDE field to gate plasticity.

This is exactly the “continuous field \(\leftrightarrow\) discrete cell” coupling pattern emphasized in the course.

## 4. Key assumptions (kept simple on purpose)
- Rate neurons (not full spiking) are enough to show learning and routing.
- A single modulator captures reward-gated plasticity.
- Diffusion on a 2D sheet is an acceptable cortical-tissue abstraction.
- Plasticity is controlled by local modulator sensed at each neuron’s position.

## 5. Expected outputs / figures
- Heatmaps of \(M(x,y,t)\) showing diffusion from reward pulse.
- Connectivity or adjacency plot before vs after learning.
- Activity propagation (space-time plot) showing stronger Input \(\rightarrow\) Output after training.
- Control experiments: with \(\eta=0\) (no learning) or with fast clearance (high \(V_{\max}\)).

## 6. Where the “other concepts” fit
- **Enzyme kinetics:** Michaelis–Menten clearance term of \(M\); modulate \(V_{\max}\) as “drug/enzyme therapy”.
- **Metabolic networks:** add an “energy” variable \(E\) per neuron; activity consumes \(E\), limiting plasticity.
- **Biotechnology:** simulate optogenetic stimulation input \(I_i(t)\) + pharmacology changing clearance.
- **Nonlinear dynamics:** sigmoidal \(\phi\), Hill gating \(g\), potential bistability in receptor module.
- **Evolutionary dynamics (optional):** optimize \(D_M, V_{\max}, \eta\) to maximize routing performance.
- **Cell-based modelling:** neurons are discrete; potentially add glia agents that uptake \(M\).

---

# IDEA 2 — Excitable media PDE (FitzHugh–Nagumo/Barkley) + coupling (metabolism/modulator) + cell agents

## 1. Motivation
The cortex (and many biological tissues) can support **travelling waves** of activity. Excitable-media models reproduce key qualitative features of action-potential propagation (threshold, wavefront, refractory period) without the complexity of full Hodgkin–Huxley.

This idea aims to show:
- Signal transmission in space via waves.
- How chemical/metabolic constraints modulate propagation.
- (Optional) simple wave-based computation: wave collision/annihilation, gating via refractory regions.

## 2. Goal (what we will demonstrate)
A minimal set of deliverables:

- **Coupled PDEs:** Excitable system \((u,v)\) in 2D; plus at least one additional coupled PDE (e.g., metabolic resource \(E\) or neuromodulator \(M\)) to satisfy “coupled PDEs” explicitly.
- **Cell-based component:** discrete agents controlling local parameters or acting as sinks/sources (e.g., glia that absorb activator or release inhibitor).
- **Outcome:** show robust wave propagation from a stimulus point, and how coupling (energy depletion, modulator inhibition, or glial uptake) can block or shape waves.

## 3. Theoretical background

### 3.1 FitzHugh–Nagumo reaction–diffusion system
A standard 2-field excitable model:

$$
\frac{\partial u}{\partial t} = D_u \nabla^2 u + u - \frac{u^3}{3} - v + I(x,y,t)
$$
$$
\frac{\partial v}{\partial t} = D_v \nabla^2 v + \varepsilon (u + a - b v)
$$

- \(u\): activator (membrane potential proxy).
- \(v\): recovery variable (refractoriness proxy).
- \(I(x,y,t)\): external stimulus current (pulse input).
- \(a,b,\varepsilon\): parameters controlling excitability.

Often \(D_v\) is small (or zero) compared to \(D_u\).

**Boundary conditions:** no-flux is again a reasonable tissue patch approximation.

### 3.2 Coupled metabolic/resource PDE (adds physiology + coupling)
Define an energy/resource field \(E(x,y,t)\) (e.g., ATP availability proxy):

$$
\frac{\partial E}{\partial t} = D_E \nabla^2 E + k_{prod}(E_0 - E) - k_{cons}\, h(u)\, E
$$

- \(k_{prod}\): recovery/production rate toward baseline \(E_0\).
- \(k_{cons}\): consumption rate, depends on activity via \(h(u)\) (e.g., \(h(u)=u^2\) or a Hill function).

Couple energy back into excitability by modulating input gain or threshold:

$$
\frac{\partial u}{\partial t} = D_u \nabla^2 u + \big(1+\alpha_E E\big)\left(u - \frac{u^3}{3}\right) - v + I
$$

or change parameter \(a\rightarrow a(E)\). This makes the PDEs *explicitly coupled*.

### 3.3 Cell-based layer (agents as local heterogeneity)
Introduce discrete agent types:

- **Neuron agents:** define locations that receive stimulus pulses \(I(x,y,t)\), or locally have different excitability parameters.
- **Glia agents:** uptake activator or release inhibitor locally (source/sink terms):

$$
\frac{\partial u}{\partial t} = \dots - \sum_i \kappa_i\, \delta(x-x_i,y-y_i)\, u
$$

In implementation, replace \(\delta\) with “affects one grid cell” or a small Gaussian kernel.

### 3.4 Optional enzyme kinetics / biotech intervention
Add a diffusing inhibitor/drug \(M(x,y,t)\) that reduces excitability:

$$
\frac{\partial M}{\partial t} = D_M \nabla^2 M + S_M - \frac{V_{\max} M}{K_M + M}
$$

Then couple to \(u\) by lowering stimulus gain or increasing recovery:
$$
 I_{eff} = \frac{I}{1 + \alpha_M M}
$$

This yields a direct biotech narrative: “apply drug/enzyme that blocks waves”.

## 4. Key assumptions
- Excitable media captures macroscopic wave phenomena.
- We don’t claim exact action-potential biophysics; we claim qualitative wave mechanics.
- Agents represent localized heterogeneity (neurons/glia) rather than full cell morphologies.

## 5. Expected outputs / figures
- Time snapshots of \(u(x,y,t)\) showing wavefronts.
- Phase-plane plot at a point showing excitable dynamics (\(u\) vs \(v\)).
- Demonstrate “propagation failure” when metabolic resource is depleted or inhibitor applied.
- Compare homogeneous vs heterogeneous agent layouts.

## 6. Where the “other concepts” fit
- **Nonlinear dynamics:** excitable nullclines, thresholds, refractory effects.
- **Reaction–diffusion:** intrinsic to FHN + optional metabolic/drug field.
- **Metabolic networks:** resource field limits excitability.
- **Enzyme kinetics / biotechnology:** MM clearance of inhibitor/drug; simulate interventions.
- **GRNs:** optionally add slow \(R\) gene-expression field that sets excitability parameter \(a\).
- **Evolutionary dynamics (optional):** optimize parameters for fastest reliable transmission or minimal energy use.
- **Cell-based modelling:** agent layer provides sinks/sources and heterogeneity.

---

# Possible extensions (computation-focused, fits Idea 1 & 2)

## Extension A (fits **Idea 1**): Reward-gated routing = simple reinforcement learning in space

**Motivation:** The “neuromodulator learning wave” can be reframed as a minimal **reinforcement learning** mechanism: only when a **reward** arrives (high neuromodulator) do synapses that were recently “eligible” get strengthened. This clarifies *what computation is being learned*: reliable routing from an input region to a goal/output region.

### Minimal extra state: eligibility trace
For each synapse \(w_{ij}\), define an eligibility \(e_{ij}(t)\) capturing recent pre/post co-activity:

$$
\tau_e \frac{d e_{ij}}{dt} = -e_{ij} + a_i a_j
$$

Then update weights using a 3-factor rule gated by local neuromodulator \(M_i\):

$$
\frac{d w_{ij}}{dt} = \eta\, g(M_i)\, e_{ij} - \lambda_w w_{ij}
$$

- “**Credit assignment**” emerges: only synapses that were active **before** reward (high \(e\)) are modified **when** reward arrives (high \(M\)).

### Experimental demo
- Place two potential targets **B** (rewarded) and **C** (not rewarded).
- Stimulate **Input A** repeatedly.
- Deliver reward (modulator source \(S_M\)) only when activity reaches B (or externally at B).
- Show that over training, activity from A preferentially routes to B more often/stronger than to C.

### Figures to add
- Before/after activity routing probability A→B vs A→C.
- Weight heatmap or network graph showing strengthened pathways.
- Modulator field \(M(x,y,t)\) snapshot during reward.

**Course concept hooks:** nonlinear learning dynamics, GRN-style Hill gating \(g(M)\), enzyme kinetics via MM clearance of \(M\), plus the same coupled PDE ↔ agent coupling.

---

## Extension B (fits **Idea 2**): Wave-based logic gates (excitable media computing)

**Motivation:** Excitable-media PDEs can implement simple **logic** because wavefronts interact nonlinearly:
- collisions can annihilate or reinforce
- refractory regions prevent passage (diode-like behavior)
- geometry can emulate AND/OR junctions

This gives a direct “brain-like computation” message: **information can be processed by spatiotemporal patterns**, not only by wired digital circuits.

### Implementation sketch
Use the same excitable PDE (e.g., FitzHugh–Nagumo):

- Define two input ports A and B, and one output port O.
- Encode bit=1 as a stimulus pulse \(I(x,y,t)\) at the port at \(t=t_0\). Encode bit=0 as no pulse.

**AND gate (typical design):** Only if two waves arrive simultaneously at a junction does the combined excitation exceed threshold and propagate to output.

This is easiest if we introduce a spatially varying excitability parameter \(a(x,y)\) or threshold region at the junction (heterogeneity), potentially implemented by agent “gate cells” or a fixed mask.

### Optional coupling to satisfy ‘coupled PDEs’ even more clearly
Add a diffusing inhibitory/drug field \(M(x,y,t)\) (or metabolic field \(E(x,y,t)\)) that dynamically changes excitability in parts of the gate:

$$
\frac{\partial M}{\partial t} = D_M \nabla^2 M + S_M - k_M M
$$

and couple via \(I_{eff}=I/(1+\alpha_M M)\) or \(a \to a+\alpha_M M\).

This lets us “train/tune” a gate by adjusting inhibitor distribution (biotech narrative) or by slowly adapting \(a(x,y)\) using a plasticity-like rule.

### Figures to add
- A panel with 4 rows (00, 01, 10, 11) showing wave snapshots and whether O fires (a truth table).
- Parameter sensitivity: show how changing \(a\) or inhibitor \(M\) toggles gate reliability.

**Course concept hooks:** strong nonlinear dynamics, coupled PDEs (u,v plus M/E), reaction–diffusion, and cell-based agents providing localized sinks/sources (glia-like control of propagation).

---

# Choosing between Idea 1 vs Idea 2 (practical criteria)
- Choose **Idea 1** if you want the headline to be **learning / plasticity** and “reward-gated” adaptation.
- Choose **Idea 2** if you want the headline to be **signal propagation / computation by waves** with very nice spatiotemporal visuals.

Both can be extended with the same “concept modules” (enzyme clearance, metabolism, gene expression, agent sinks/sources).

---

# Implementation plan & timeline (1.5-week sprint)

We’ll assume the sprint is **10 days** (including weekends) and split the work into three parallel tracks. This matches the rough plan implied by `assignment3/assignment.md` (finalize topic early; then equations; then simulations+results).

## Roles (consistent with `assignment3/contributions.md`)
- **Anna:** Model development & equations (derivations, assumptions, parameter choices).
- **Konstantinos:** Figures and tables (visualization pipeline, animations, publication-ready plots).
- **Tim:** Results and implications (experiment design, ablations, interpretation, writing).

> Note: Everyone codes; the split mainly defines *ownership* so work progresses in parallel.

## Day-by-day timeline

### Day 1 (Today): finalize model choice + repo scaffolding
- **All:** pick Idea 1 *or* Idea 2 as primary, decide minimal deliverables.
- **Anna:** produce “equations & assumptions” markdown draft (copy into `assignment3/equations.md`).
- **Konstantinos:** create plotting templates: heatmaps, quiver/streamplots, time-series, panel layouts.
- **Tim:** define 2–3 experimental questions and success metrics.

### Day 2–3: get the coupled PDE core running + one validation figure
- **Anna:** parameter sanity ranges; boundary conditions; stability considerations for explicit scheme.
- **Konstantinos:** implement 2D Laplacian + no-flux boundaries; produce first GIF/frames.
- **Tim:** set up experiment scripts/notebook to reproduce runs and compare settings.

Validation targets:
- Idea 1: diffusion from point/patch source creates expected gradient.
- Idea 2: pulse stimulus creates a travelling wave (no blow-up).

### Day 4–5: add cell-based coupling (agents) + verify coupling works
- **Anna:** define agent-to-field coupling terms (sources/sinks), plasticity/energy coupling.
- **Konstantinos:** visualization overlays (agent positions on top of field heatmap).
- **Tim:** run ablations (no coupling vs coupling) and summarize differences.

Validation targets:
- Idea 1: modulator reaches agents and gates plasticity (weights change only where \(M\) high).
- Idea 2: agent sinks/sources reshape wave patterns.

### Day 6–7: run “main story” experiments + polish the narrative
- **Anna:** finalize theoretical background text and parameter table.
- **Konstantinos:** generate final panels: before/after learning (Idea 1) or wave propagation and modulation (Idea 2).
- **Tim:** write results/discussion, limitations, and biological interpretation.

### Day 8–9: robustness checks + optional extra concept
- Add one extra concept if time allows:
  - Idea 1: enzymatic clearance term (MM) vs linear decay comparison.
  - Idea 2: metabolic resource field \(E\) coupling and energy depletion experiment.
  - Optional “evolutionary” sweep: random search to optimize a metric.

### Day 10: writeup integration + final QA
- **All:** merge sections into final report structure.
- **Konstantinos:** ensure figures have consistent style, captions, units, and legends.
- **Tim:** ensure claims are supported; add “what we learned” and future work.

---

# Technical plan (lightweight, realistic)

## Numerics (suggested)
- 2D finite differences on a grid (NumPy).
- Time integration: explicit Euler for reaction terms; for diffusion, explicit (small dt) or semi-implicit if needed.
- No-flux boundaries via mirrored edges / Neumann stencil.

## Reproducibility
- Single `scripts/` entry point per experiment.
- Save parameters and random seeds alongside outputs.
- Save figures into `assignment3/results/`.

---

# References (starting list)

## Reaction–diffusion & excitable media
- Murray, J. D. (2002/2003). *Mathematical Biology I/II*. Springer. (reaction–diffusion fundamentals)
- FitzHugh, R. (1961). Impulses and physiological states in theoretical models of nerve membrane. *Biophysical Journal*.
- Nagumo, J., Arimoto, S., & Yoshizawa, S. (1962). An active pulse transmission line simulating nerve axon. *Proceedings of the IRE*.

## Neural plasticity & neuromodulation
- Dayan, P., & Abbott, L. F. (2001). *Theoretical Neuroscience*. (rate models, learning rules)
- Gerstner, W., Kistler, W. M., Naud, R., & Paninski, L. (2014). *Neuronal Dynamics*. (plasticity + dynamical systems)

## Enzyme kinetics / Hill functions / GRNs
- Michaelis, L., & Menten, M. L. (1913). Die Kinetik der Invertinwirkung. (Michaelis–Menten kinetics)
- Alon, U. (2019). *An Introduction to Systems Biology: Design Principles of Biological Circuits*. (Hill-function GRN motifs)

## Coupled PDE + agent-based coupling patterns
- Fellermann, H., et al. (2017). Spatio-temporal modeling of the origin of life. *Artificial Life*. (continuous fields coupled to discrete agents; delta-function style coupling)

---

## Next step (decision point)
At the next group meeting, decide:

1. Primary direction: **Idea 1** *or* **Idea 2**.
2. Minimum coupled PDE set: (a) just modulator + enzyme field, or (b) FHN + metabolic resource field.
3. Primary output: (a) “learning / rewiring” animation, or (b) “signal wave + modulation” animation.

