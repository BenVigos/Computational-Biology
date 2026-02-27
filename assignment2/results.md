---
geometry:
- margin=1in
header-includes:
  - \usepackage{float}
  - \usepackage{subcaption}
---

# Assignment 2: Figures and tables

## Patient Alpha

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{plots/alpha_4_hill_functions.png}
\caption{Activation and Inhibition Hill Functions}
\end{figure}

\begin{figure}[h]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_3_phase_space.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_3_phase_space_hijack.png}}
\caption{Phase Portrait of Proteins and mRNAs concentrations for Patient Alpha}
\end{figure}


\begin{figure}[h]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_1_proteins_mrna_time.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_1_proteins_mrna_time_hijack.png}}
\caption{Concentration of Proteins and mRNAs overtime for Patient Alpha}
\end{figure}

\begin{figure}[h]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_5_rates.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_5_rates_hijack.png}}
\caption{Translation and Transcription Rates overtime for Patient Alpha}
\end{figure}

## Patient Beta

![Phase Portrait of the Protein concentrations of Patient Beta](plots/sdevelo_phase_portrait.png)

![Comparison of the (average) concentration of protein over time from ODE versus 10 SDEVelo simulations with 95% confidence interval](plots/sdevelo_concentrations_over_time.png)

<!-- ![Average concentration of (un)spliced RNA and protein over time from 10 SDEVelo simulations with 95% confidence interval](plots/sdevelo_versus_ode.png) -->

![Stream plot with equilibrium point and nullclines for interaction between metabolite R, and gorwth-promoting enzyme E](plots\bonus_phase_portrait.png) 
