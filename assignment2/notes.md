# Notes

## Ben Notes (not sure where they go yet, but here for now):
In mechanism 2 we have a lot of misfolded proteins (which we don't model) 
and a lot of unspliced RNA (which we do model) because of the large amount of unspliced mRNA. 
 Even though we don't model the misfolded proteins directly, we can infer that they contribute to the 
aggressive phenotype observed in Patient Beta.)

In Viterbi algorithm the 1st state is uncertain. If we wanted to be more rigorous, we could run the algorithm multiple times with 
different initial states (stochastic) and see if we consistently get the same final state sequence. 
However, given the clear dominance of exons in Patient Alpha and introns in Patient Beta, we can be reasonably confident in our 
conclusions about the underlying mechanisms. Another option would be to use the reverse Viterbi algorithm, which starts from the end 
of the sequence and works backward, to see if it converges on the same state sequence.

## Equations and Model

### 1. Biological System Overview

- Two-gene regulatory network:
  - Gene A → tumor suppressor
  - Gene B → oncogene
- Proteins mutually regulate each other.
- Cancer disrupts this regulation via:
  - Mechanism I: Transcriptional Hijack
  - Mechanism II: Splicing Sabotage

Goal:
Identify mechanism from RNA data and model resulting dynamics.

---

### 2. Mechanism Identification – Hidden Markov Model

Why HMM?
- Observed: RNA sequences
- Hidden: Regulatory mechanism
- Need probabilistic classification framework

We apply the Viterbi algorithm to:
- Compute most likely hidden state sequence
- Infer which disruption mechanism is active

Assumptions:
- The system transitions between discrete regulatory states
- Emission probabilities correctly represent sequence likelihood
- Observations are conditionally independent given the state

Modeling decision:
Mechanism classification must precede dynamic modeling.

---

### 3. Mechanism I → Deterministic ODE Model

Biology:
- Protein A loses transcriptional repression of Gene B
- No defect in splicing
- mRNA fully processed

Why ODEs?
- Dynamics are continuous
- Molecular noise is not dominant
- System approximated by average concentrations

Core structure:
- mRNA production - degradation
- Protein translation - degradation
- Regulation via Hill functions (nonlinear threshold behavior)

Key assumption:
- Concentrations are large enough for deterministic approximation
- Binding interactions reach quasi-steady-state
- Splicing is efficient and not rate-limiting
- Transcriptional repression is completely disabled under hijack

Interpretation:
Predict steady-state shifts toward oncogene dominance.

---

### 4. Mechanism II → Stochastic SDEVelo Model

Biology:
- Protein A disrupts splicing
- Pre-mRNA accumulates
- Aggressive tumor phenotype

Why SDE model?
- Splicing is stochastic
- High transcriptional variability
- Need to track:
  - Unspliced RNA (U)
  - Spliced RNA (S)
  - Protein (P)

Structure:
- Transcription - regulated splicing + noise
- Translation - degradation

Noise modeled via Wiener processes.

Key assumptions:
- Splicing rate depends on protein concentration
- Transcription occurs in bursts (sigmoidal activation)
- Molecular fluctuations significantly affect dynamics
- Noise is Gaussian and temporally uncorrelated

Interpretation:
System may exhibit variability-driven instability and rapid shifts.

---

### 5. Downstream Metabolic Model

Equations describe:
- Resource growth and consumption
- Enzyme decay and activation

Why include it?
- Connect molecular dysregulation to cellular phenotype
- Study long-term stability

Assumptions:
- Resource-enzyme interaction follows mass-action kinetics
- Parameters remain constant over time
- No external resource input beyond modeled growth term

Stability analysis shows:
- One unstable saddle
- One oscillatory center
- Suggests persistent metabolic cycling

---

### 6. Overall Modeling Logic

1. Classify mechanism using HMM.
2. Choose framework based on biology and noise dominance:
   - ODE → deterministic transcription defect
   - SDE → stochastic splicing defect
3. Use Hill functions for nonlinear regulation.
4. Analyze stability to predict long-term cell fate.

Next video:
Parameter tables, simulations, and phase portraits.

## Figures and tables

Figure 1: Illustrates the hijacking in mechanism I. Usually there is Inhibition of Gene B by Protein A, and Activation of Gene A by Protein B. However, the hijacking caused by the studied cancer means there is no inhibition of Gene B by Protein A.

Figure 2: This Phase Space of Protein A and B and mRNA A and B aims to illustrate the effect of the hijacking. Protein B activates production of Protein A, which reduces the amount of Protein B. This feedback loop ensures that cell division by Protein B remains controlled. With hijacking (b) we see that protein A and B both increase, which causes uncontrolled cell division.

Figure 3 shows the Time Evolution of the Protein and mRNA concentrations. In a healthy state the increase in Protein B causes a slightly delayed increase in Protein A (seen by the shift of the peak to the right), which then stops mRNA B production and decreases Protein B concentration.

In the Hijacked state the concentrations both reach their steady state equillibrium at circa 2.3 molar.

Figure 4 shows the Translation and Transcription rates. These follow a very similar shape to Figure 2 as they are the first derivative of Figure 2.

Healthy states results were validated against Polynikis A, Hogan SJ, di Bernardo M. Comparing different ODE modelling approaches for gene regulatory networks. *J Theor Biol*. 2009;261(4):511-530. doi:10.1016/j.jtbi.2009.07.040


Figure 5. Aims to illustrate the behaviour of Mechanism II, modelled using SDEVelo, compared to the healthy state, modelled using ODE. The SDE has clear stochastic behaviour but follows a somewhat similar trajectory to the Healthy ODE.

Figure 6. Aims to show the effect of Mechanism II, by showing the concentrations of unspliced and spliced RNA and protein over time. This shows that, as we saw in Figure 5, the concentration of Protein B is kept below 0.5 molar. However, we can see a massive buildup of unspliced RNA for Gene B.

Figure 7 shows the Phase Portrait with two distinct fixed points. At (0,0) we have a saddle point that is attractive in the E direction, but repulsive in the R direction. Around (1.1, 1.8) we see a stable fixed point: a stable spiral. Nullclines at dE/dt and dR/dt show no change in that axis.


## Results and implications

Discussion Draft 
Patient Alpha: Transcriptional Hijack and Uncontrolled Proliferation
The application of the Viterbi algorithm to Patient Alpha's RNA-seq data revealed 
a sequence composed entirely of exons. This indicates that mRNA splicing proceeds 
normally; however, the regulatory network is disrupted at the transcriptional level.
We diagnosed this as Mechanism I (Transcriptional Hijack), where Protein A fails to 
inhibit the transcription of Gene B.

In a healthy cell, this two-gene network acts as a strict negative feedback loop, maintaining protein and mRNA concentrations at safe, 
controlled levels (which we also compared with the homeostatic baseline from [insert paper name]). However, because the tumor cell 
bypassed Protein A's inhibition, we observe a runaway process in our ODE model. Starting from initial concentrations, mRNA B and 
Protein B surge, eventually reaching a pathologically high steady state. Consequently, Protein A—which is promoted by Protein B 
is also overexpressed in a futile attempt to restore balance.

Cellular Fate of Patient Alpha: What do these behaviors imply for the cellular fate? The continuous, high-level expression of Protein 
B (the PROLIFERATOR) locks the cell into a state of unchecked, continuous division. The cellular fate here is classic, uncontrolled 
tumor growth, driven by the sheer overexpression of oncogenes, which is typical of early-stage or standard proliferative cancers.

Patient Beta: Splicing Sabotage and Aggressive Toxicity
For Patient Beta, the Viterbi algorithm returned a sequence dominated by introns, pointing to a severe disruption in RNA maturation. 
We identified this as Mechanism II (Splicing Sabotage). By applying the stochastic SDEVelo model, we successfully captured the massive 
buildup of unspliced pre-mRNA.

Because the pre-mRNA cannot mature properly, the translation into functional Protein B is severely hampered, likely resulting in 
misfolded or non-functional products. Consequently, functional Protein B levels dampen over time. Since Protein B is required to 
promote Protein A, the Guardian protein also stabilizes at much lower levels. Overall, the phase space trajectories appear highly 
chaotic and noisy due to biological stochasticity, standing in stark contrast to the smooth, deterministic curves of Patient Alpha's 
ODE model.
    
Cellular Fate of Patient Beta:
Despite having lower levels of functional Protein B, the cellular fate here is an extremely aggressive and fast-growing tumor. 
The massive accumulation of unspliced RNA creates a highly toxic intracellular environment. This "intron retention" and RNA buildup 
likely sequesters critical cellular machinery or triggers alternative, highly malignant pathways, leading to the severe clinical 
phenotype observed in this patient.

Downstream Metabolic Effects (Bonus)
Finally, to understand how the tumor sustains this growth, we modeled the downstream metabolic interactions between a cellular resource
(R) and a growth-promoting enzyme (E). The stability analysis of this system returned two fixed points. The origin (0,0) was identified
as a saddle point, meaning total extinction is unstable.
More importantly, the coexistence equilibrium $(R^*, E^*)$ yielded a Jacobian matrix with purely imaginary eigenvalues, classifying it 
mathematically as a "center". Biologically, this means the metabolic system does not settle into a static state, but rather enters a 
regime of sustained, stable oscillations. This cyclic dynamic allows the tumor to continuously and efficiently harvest metabolic
resources without completely exhausting its host environment, perfectly fueling its relentless growth.


Coupled protein, slide 126