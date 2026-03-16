---
geometry:
- margin=0.5in
- paperwidth=210mm
- paperheight=210mm
header-includes:
  - \usepackage{float}
  - \usepackage{graphicx}
---
# Assignment 3: Figures and videos

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{angiogenesis-3612748292.jpg}
\caption{Schematic representation of the tumor microenvironment and angiogenesis process. The figure illustrates the interactions between tumor cells, endothelial cells, and the extracellular matrix, highlighting key signaling pathways such as VEGF and HIF-1$\alpha$ that drive angiogenic sprouting and vascular remodeling in response to hypoxic conditions within the tumor mass. Image Source:Tocris Bioscience. Angiogenesis in cancer research product guide, 
3 Ed., (2015),  https://www.tocris.com/literature/product-guides/angiogenesis-in-cancer }
\label{fig:angiogenesis-schematic}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{results/comparison_figures/figure_A_tumor_vs_angiogenesis.png}
\caption{Figure A. Comparison of tumor growth in simulations without angiogenesis and including angiogenesis. Panel (a) shows the temporal evolution of tumor size, reported as both the effective tumor radius and the mean radial distance of tumor cells from the tumor centroid. Panel (b) shows the corresponding tumor cell composition over time, separating total tumor-like cells into normoxic, hypoxic, and necrotic subpopulations. Panel (c) shows the hypoxic fraction, highlighting how the emergence of angiogenesis changes the degree of oxygen stress experienced by the tumor.}
\label{fig:comparison-a}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{results/comparison_figures/figure_B_endothelial_activation.png}
\caption{Figure B. Comparison of endothelial activation in simulations using a constant fully hypoxic tumor with and without HIF-1$\alpha$ signaling enabled. Panel (a) shows the time evolution of VEGF and HIF-1$\alpha$-associated signaling, illustrating how the inclusion of HIF-1$\alpha$ alters (or does not alter) the molecular pro-angiogenic response. Panel (b) shows the numbers of active and inactive neovascular cells, providing a direct measure of how signaling differences translate into endothelial activation. Panel (c) shows tumor oxygen levels together with the minimum distance between the tumor and nearby sprouts.}
\label{fig:comparison-b}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{results/comparison_figures/figure_C_combined_tumor_dynamics.png}
\caption{Figure C. Integrated comparison of tumor-only growth and the full tumor--angiogenesis model. Panel (a) shows the temporal evolution of total tumor-like cells together with the normoxic, hypoxic, and necrotic subpopulations. Panel (b) shows VEGF and HIF-1$\alpha$ signaling dynamics. Panel (c) shows tumor oxygen levels together with the minimum tumor-to-vessel.}
\label{fig:comparison-c}
\end{figure}

# Supplementary Videos
4 supplementary videos are included to visually illustrate the dynamic processes involved in each of the 3 performed simulations.

- Video 1: Tumor growth without angiogenesis
- Video 2: Angiogenesis without tumor growth
- Video 3: Tumor growth with angiogenesis
