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

## Figures and tables

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