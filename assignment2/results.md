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

Most likely sequence of hidden states for AGCGC: ['?', 'E', 'E', 'E', 'E']

Table: Viterbi probabilities for Patient Alpha's observed sequence AGCGC, where E=Exon and I=Intron

|         | t=0     | t=1      | t=2      | t=3      |          t=4 |
|---------|--------:|---------:|---------:|---------:|-------------:|
| Exon    | 0.125   | 0.028125 | 0.006328 | 0.001424 | **0.000320** |
| Intron  | 0.200   | 0.008000 | 0.000960 | 0.000038 |     0.000021 |


Table: Backpointer indices for Patient Alpha's observed sequence AGCGC, where E=Exon and I=Intron

|         | t=0 | t=1 | t=2 | t=3 | t=4 |
|---------|:---:|:---:|:---:|:---:|:---:|
| Exon    |  ?  |  **E**  |  **E**  |  **E** |  **E** |
| Intron  |  ?  |  I  |  I  |  I  |  E  |

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

Most likely sequence of hidden states for AUUAU: ['?', 'I', 'I', 'I', 'I']

Table: Viterbi probabilities for Patient Beta's observed sequence AUUAU, where E=Exon and I=Intron

|         | t=0     | t=1      | t=2      | t=3      |          t=4 |
|---------|--------:|---------:|---------:|---------:|-------------:|
| Exon    | 0.125   | 0.028125 | 0.006328 | 0.001424 |     0.000328 |
| Intron  | 0.200   | 0.064000 | 0.020480 | 0.006554 | **0.002097** |


Table: Backpointer indices for Patient Beta's observed sequence AUUAU, where E=Exon and I=Intron

|         | t=0 | t=1 | t=2 | t=3 | t=4 |
|---------|:---:|:---:|:---:|:---:|:---:|
| Exon    |  ?  |  E  |  E  |  E  |  I  |
| Intron  |  ?  | **I** | **I** | **I** | **I** |

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