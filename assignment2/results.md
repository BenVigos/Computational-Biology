---
geometry:
- margin=1in
- paperwidth=210mm
- paperheight=210mm
header-includes:
  - \usepackage{float}
  - \usepackage{subcaption}
  - \usepackage{graphicx}
---
# Assignment 2: Figures and tables

## Patient Alpha

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{plots/alpha_4_hill_functions.png}
\caption{Activation and Inhibition Hill Functions}
\end{figure}

\begin{figure}[H]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_3_phase_space.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_3_phase_space_hijack.png}}
\caption{Phase Portrait of Proteins and mRNAs concentrations for Patient Alpha}
\end{figure}

\begin{figure}[H]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_1_proteins_mrna_time.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_1_proteins_mrna_time_hijack.png}}
\caption{Concentration of Proteins and mRNAs overtime for Patient Alpha}
\end{figure}

\begin{figure}[H]
\centering
\subfloat[Healthy State]{\includegraphics[width=0.45\textwidth]{plots/alpha_5_rates.png}}
\hfill
\subfloat[Hijacked State]{\includegraphics[width=0.45\textwidth]{plots/alpha_5_rates_hijack.png}}
\caption{Translation and Transcription Rates overtime for Patient Alpha}
\end{figure}

\newpage
## Patient Beta

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{plots/sdevelo_phase_portrait.png}
\caption{Phase Portrait of the Protein concentrations of Patient Beta}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[height=0.8\textwidth]{plots/sdevelo_concentrations_over_time.png}
\caption{Comparison of the (average) concentration of protein over time from ODE versus 10 SDEVelo simulations with 95\% confidence interval}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{plots/bonus_phase_portrait.png}
\caption{Stream plot with equilibrium point and nullclines for interaction between metabolite R, and growth-promoting enzyme E}
\end{figure}