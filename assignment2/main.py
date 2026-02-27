import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.interpolate import griddata
from matplotlib.lines import Line2D
import pandas as pd


def get_letter_index(letter):
    """
    Helper function to convert a letter (A, C, G, U) to an index (0, 1, 2, 3).
    :param letter: A character representing a nucleotide.
    :return: An integer index corresponding to the nucleotide.
    """
    index = {"A": 0, "U": 1, "G": 2, "C": 3}
    return index.get(letter, ValueError)


def get_exon_intron_state(state):
    """
    Helper function to convert a state index to a string representation (Exon or Intron).
    :param state: An integer index representing the state.
    :return: A string 'Exon' or 'Intron' corresponding to the state.
    """
    index = {0: "E", 1: "I"}
    return index.get(state, "?")


def viterbi_algorithm(sequence, A, B, P):
    """
    Implements the Viterbi algorithm to find the most likely sequence of hidden states.
    :param sequence: The observed sequence (e.g., 'AGCGC').
    :param A: Transition probability matrix.
    :param B: Emission probability matrix.
    :param P: Initial state probabilities.
    :return: The most likely sequence of hidden states.
    """
    n_states = len(P)
    T = len(sequence)

    # Initialize the Viterbi matrix and backpointer
    V = np.zeros((n_states, T))
    backpointer = np.zeros((n_states, T), dtype=int) - 1

    # Initialization step
    for s in range(n_states):
        V[s, 0] = P[s] * B[s, get_letter_index(sequence[0])]

    # Recursion step
    for t in range(1, T):
        for s in range(n_states):
            max_prob = -1
            max_state = 0
            for s_prev in range(n_states):
                prob = (
                    V[s_prev, t - 1]
                    * A[s_prev, s]
                    * B[s, get_letter_index(sequence[t])]
                )
                if prob > max_prob:
                    max_prob = prob
                    max_state = s_prev
            V[s, t] = max_prob
            backpointer[s, t] = max_state

    # Termination step
    best_path_prob = -1
    best_last_state = 0
    for s in range(n_states):
        if V[s, T - 1] > best_path_prob:
            best_path_prob = V[s, T - 1]
            best_last_state = s

    # Backtrack to find the best path
    best_path = [best_last_state]
    for t in range(T - 2, 0, -1):
        best_path.append(backpointer[best_path[-1], t])

    best_path.append(-1)
    best_path.reverse()

    for i, s in enumerate(best_path):
        best_path[i] = get_exon_intron_state(s)

    backpointer_symbolic = np.empty(backpointer.shape, dtype=object)

    for i, r in enumerate(backpointer):
        for j, s in enumerate(r):
            backpointer_symbolic[i, j] = get_exon_intron_state(s)


    return best_path, V, backpointer_symbolic


def hill_activation(P, theta, n):
    """Calculate the Hill activation function.

    :param P: Concentration of the regulator protein.
    :param theta: Threshold concentration for activation.
    :param n: Hill coefficient.
    """
    return P**n / (theta**n + P**n)


def hill_inhibition(P, theta, n):
    """Calculate the Hill inhibition function.
    :param P: Concentration of the regulator protein.
    :param theta: Threshold concentration for inhibition.
    :param n: Hill coefficient."""
    return 1 - hill_activation(P, theta, n)


def patient_alpha_ode(t, y, hijack=None):
    """Calculate the derivatives of mRNA and protein for Patient Alpha.
    :param t: Time variable (not used in this ODE but required by solve_ivp).
    :param y: A list or array containing the current values of [ra, rb, P_a, P_b].
    """
    ra, rb, P_a, P_b = y
    dra_dt = m_a * hill_activation(P_b, theta_b, n_b) - gamma_a * ra
    if hijack:
        drb_dt = m_b * 1 - gamma_b * rb  # no inhibition on rb
    else:
        drb_dt = m_b * hill_inhibition(P_a, theta_a, n_a) - gamma_b * rb
    dPa_dt = k_pa * ra - delta_pa * P_a
    dPb_dt = k_pb * rb - delta_pb * P_b
    return [dra_dt, drb_dt, dPa_dt, dPb_dt]

def concentration_overtime(sol, t, hijack=False):
    """Plot the time evolution of mRNA and protein concentrations.
    :param sol: Solution object.
    :param t: Time points corresponding to the solution.
    :param hijack: Boolean indicating whether to plot the hijacked scenario or not.
    """
    ra, rb, pA, pB = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    # pA, pB over time
    plt.figure(figsize=(8, 6))
    plt.plot(t, pA, label='Protein A', color = "tab:blue", linewidth=1.5)
    plt.plot(t, pB, label='Protein B', color = "tab:red", linewidth=1.5)
    plt.plot(t, ra, label='mRNA A',linewidth=1.5, color = "tab:blue", linestyle="--")
    plt.plot(t, rb, label='mRNA B', linewidth=1.5, color = "tab:red", linestyle="--")
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("Concentration [M]", fontsize=12)
    plt.title("Time Evolution of Protein and mRNA Concentrations", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("plots/alpha_1_proteins_mrna_time_hijack.png", dpi=300)
    else:
        plt.savefig("plots/alpha_1_proteins_mrna_time.png", dpi=300)
    plt.clf()

def phase_plots(sol, t, hijack=False):
    """"Plot the phase space of protein concentrations and mRNA concentrations.
    :param sol: Solution object.
    :param t: Time points corresponding to the solution.
    :param hijack: Boolean indicating whether to plot the hijacked scenario or not.
    """
    ra, rb, pA, pB = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    # phase plot pA vs pB and mRNA_A vs mRNA_B
    plt.figure(figsize=(8, 6))
    plt.plot(pA, pB, color="tab:orange", linewidth=2, label="Protein Trajectory", alpha = 0.8)
    plt.scatter([pA[0]], [pB[0]], color="tab:orange", zorder=5, label="Start Proteins (t=0)")

    plt.plot(ra, rb, color="tab:blue", linewidth=2, label="mRNA Trajectory", alpha = 0.6)
    plt.scatter([ra[0]], [rb[0]], color="tab:blue", zorder=5, label="Start mRNA (t=0)")
    plt.xlabel("$P_a, r_a$ Concentration [M]", fontsize=12)
    plt.ylabel("$P_b, r_b$ Concentration [M]", fontsize=12)
    plt.title("Phase Space: Protein A vs Protein B and mRNA A vs mRNA B", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.scatter(
        [ra[-1]],
        [rb[-1]],
        color="tab:orange",
        zorder=5,
        s=100,
        marker="*",
        label="End (mRNA Equilibrium)")
        plt.scatter(
        [pA[-1]],
        [pB[-1]],
        color="tab:blue",
        zorder=5,
        s=100,
        marker="*",
        label="End (Proteins Equilibrium)")
        plt.savefig("plots/alpha_3_phase_space_hijack.png", dpi=300)     
    else:
        plt.savefig("plots/alpha_3_phase_space.png", dpi=300)
    plt.clf()

def hill_plots(sol, t, hijack=False):
    """Plot the Hill activation and inhibition functions for the given parameters.

    :param sol: solution object containing the trajectories of mRNA and protein concentrations.
    :param t: time points corresponding to the solution.
    :param hijack: Boolean indicating whether to plot the hijacked scenario (no inhibition) or not.
    """
    ra, rb, pA, pB = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    # hill Functions
    P_vals = np.linspace(0, 1.5, 100)
    act = hill_activation(P_vals, theta_a, n_a)
    inh = hill_inhibition(P_vals, theta_a, n_a)

    plt.figure(figsize=(6, 6))
    plt.plot(
        P_vals,
        act,
        label=f"Activation ($\\theta$={theta_a}, $n$={n_a})",
        color="tab:red",
        linewidth=2,
    )
    plt.plot(
            P_vals,
            np.ones_like(P_vals),
            label=f"Hijacked, no inhibition",
            color="tab:blue", linestyle="--",
            linewidth=2,
        )

    plt.plot(
            P_vals,
            inh,
            label=f"Inhibition ($\\theta$={theta_a}, $n$={n_a})",
            color="tab:blue",
            linewidth=2,
        )
    plt.axvline(x=theta_a, color="grey", linestyle=":", label="Threshold ($\\theta$)")
    plt.xlabel("Regulator Protein Concentration [M]", fontsize=12)
    plt.ylabel("Hill Function", fontsize=12)
    plt.title("Hill Functions for Activation & Inhibition", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("plots/alpha_4_hill_functions.png", dpi=300)
    plt.clf()

def plot_rates(sol, t, hijack=False):
    """Plot the transcription and translation rates over time based on the solution trajectories.
    :param sol: solution object containing the trajectories of mRNA and protein concentrations.
    :param t: time points corresponding to the solution.
    :param hijack: Boolean indicating whether to plot the hijacked scenario (no inhibition) or not.
    """
    ra, rb, pA, pB = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    # transcription rates over time (m_a, m_b, hill functions)
    transcription_rate_A = m_a * hill_activation(pB, theta_b, n_b)
    if hijack:
        transcription_rate_B = np.full_like(t, m_b)  # 1*m_b, no inhibition on B
    else:
        transcription_rate_B = m_b * hill_inhibition(pA, theta_a, n_a)

    # translation rates over time (k_pa, k_pb, mRNA concentrations)
    translation_rate_A = k_pa * ra
    translation_rate_B = k_pb * rb

    plt.figure(figsize=(8, 6))
    plt.plot(t, transcription_rate_A, label='Gene A Transcription Rate', linewidth=1.5, color="tab:blue", linestyle='--')
    plt.plot(t, transcription_rate_B, label='Gene B Transcription Rate', linewidth=1.5, color ="tab:red", linestyle ='--')
    plt.plot(t, translation_rate_A, label='Protein A Translation Rate', linewidth=1.5, color="tab:blue")
    plt.plot(t, translation_rate_B, label='Protein B Translation Rate', linewidth=1.5, color="tab:red")
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("Rate [M/s]", fontsize=12)
    plt.title("Transcription and Translation Rates Over Time", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("plots/alpha_5_rates_hijack.png", dpi=300)
    else:
        plt.savefig("plots/alpha_5_rates.png", dpi=300)
    plt.clf()

def solve_sde_velo(dt, steps):
    """Solve the SDEVelo model
    :param dt: Time step for the simulation.
    :param steps: Number of time steps to simulate.
    """
    t = np.linspace(0, steps * dt, steps)

    # empty arrays:
    u = np.zeros((2, steps)) # unspliced mRNA for gene A and B
    s = np.zeros((2, steps)) # spliced mRNA for gene A and B
    p = np.zeros((2, steps)) # protein for gene A and B

    # initial concentrations:
    u[:, 0] = 0.8
    s[:, 0] = 0.8
    p[:, 0] = 0.8

    for i in range(steps - 1):
        alpha = c / (1 + np.exp(b * (t[i] - a)))
        activation_of_A = hill_activation(p[1, i], theta[1], n[1])
        inhibition_of_B = hill_inhibition(p[0, i], theta[0], n[0])

        beta_star = np.array([beta[0] * activation_of_A, beta[1] * inhibition_of_B])

        dW1 = np.random.normal(0, np.sqrt(dt), 2)
        dW2 = np.random.normal(0, np.sqrt(dt), 2)

        du = (alpha - beta_star * u[:, i]) * dt + sigma_1 * dW1
        u[:, i + 1] = np.maximum(0, u[:, i] + du)

        ds = (beta_star * u[:, i] - gamma * s[:, i]) * dt + sigma_2 * dW2
        s[:, i + 1] = np.maximum(0, s[:, i] + ds)

        dp = (k_P * s[:, i] - delta * p[:, i]) * dt
        p[:, i + 1] = np.maximum(0, p[:, i] + dp)

    return u, s, p, t


def sdevelo_phase_portrait(
    sol_sde, 
    label_sde, 
    ode_func=None, 
    sol_ode=None, 
    label_ode="ODE Trajectory", 
    savefig=None, 
    title="", 
    ode_args=(False), 
    padding=0.1, 
    grid_size=20,
    show_vector_field=True
):
    """Plot the phase portrait of the SDEVelo solution, optionally overlaying the ODE trajectory and vector field.
    :param sol_sde: Solution object from the SDEVelo simulation containing u, s, p, t.
    :param label_sde: Label for the SDEVelo trajectory in the legend.
    :param ode_func: Function defining the ODE system (for vector field calculation).
    :param sol_ode: Solution object from the ODE simulation to overlay on the phase portrait.
    :param label_ode: Label for the ODE trajectory in the legend.
    :param savefig: Filename to save the plot, if not None.
    :param title: Title of the plot.
    :param ode_args: Additional arguments to pass to the ODE function when calculating the vector field.
    :param padding: Fractional padding to add around the min/max of the data for vector field limits.
    :param grid_size: Number of points along each axis for the vector field grid.
    :param show_vector_field: Boolean indicating whether to show the vector field or not.
    """
    u, s, p, t = sol_sde
    pA_sde, pB_sde = np.ravel(p[0]), np.ravel(p[1])

    pA_min, pA_max = pA_sde.min(), pA_sde.max()
    pB_min, pB_max = pB_sde.min(), pB_sde.max()

    if sol_ode is not None:
        pA_ode, pB_ode = sol_ode.y[2], sol_ode.y[3]
        pA_min, pA_max = min(pA_min, pA_ode.min()), max(pA_max, pA_ode.max())
        pB_min, pB_max = min(pB_min, pB_ode.min()), max(pB_max, pB_ode.max())

    pA_range = pA_max - pA_min if (pA_max - pA_min) > 0 else 1.0
    pB_range = pB_max - pB_min if (pB_max - pB_min) > 0 else 1.0
    
    fig, ax = plt.subplots(figsize=(8, 6))

    if show_vector_field and ode_func is not None:
        pA_grid = np.linspace(max(0, pA_min - padding * pA_range), pA_max + padding * pA_range, grid_size)
        pB_grid = np.linspace(max(0, pB_min - padding * pB_range), pB_max + padding * pB_range, grid_size)
        PA, PB = np.meshgrid(pA_grid, pB_grid)
        
        U = np.zeros_like(PA)
        V = np.zeros_like(PB)
        
        def mrna_roots(vars_r, pa_val, pb_val):
            ra, rb = vars_r
            derivs = ode_func(0, [ra, rb, pa_val, pb_val], *ode_args)
            return [derivs[0], derivs[1]] 
        
        ra_guess, rb_guess = 0.8, 0.8
        
        for i in range(grid_size):
            for j in range(grid_size):
                pa_val = PA[i, j]
                pb_val = PB[i, j]
                
                try:
                    ra_ss, rb_ss = fsolve(mrna_roots, [ra_guess, rb_guess], args=(pa_val, pb_val))
                except Exception:
                    ra_ss, rb_ss = ra_guess, rb_guess
                    
                derivatives = ode_func(0, [ra_ss, rb_ss, pa_val, pb_val], *ode_args) 
                
                U[i, j] = float(np.squeeze(derivatives[2]))
                V[i, j] = float(np.squeeze(derivatives[3]))
                
        N = np.sqrt(U**2 + V**2)
        N[N == 0] = 1.0
        U_norm = U / N
        V_norm = V / N

        ax.quiver(PA, PB, U_norm, V_norm, color='lightgray', alpha=0.8, pivot='mid')

    if sol_ode is not None:
        ax.plot(pA_ode, pB_ode, label=label_ode, color="tab:orange", zorder=3, linestyle="--", alpha=0.6)
        ax.scatter(pA_ode[0], pB_ode[0], color="tab:orange", marker="o", zorder=5)
        ax.scatter(pA_ode[-1], pB_ode[-1], color="tab:orange", marker="v", zorder=5)

    ax.plot(pA_sde, pB_sde, label=label_sde, color="tab:purple", zorder=4, alpha=0.6)
    ax.scatter(pA_sde[0], pB_sde[0], color="tab:purple", marker="o", zorder=5, label="Initial state (SDE)")

    ax.set_xlabel("Protein A Concentration [M]")
    ax.set_ylabel("Protein B Concentration [M]")
    ax.grid(True, alpha=0.3)

    handles, labels = ax.get_legend_handles_labels()
    
    if show_vector_field and ode_func is not None:
        ode_proxy = Line2D([0], [0], color='gray', marker=r'$\rightarrow$', markersize=8, linestyle='None')
        handles.append(ode_proxy)
        labels.append("ODE Vector Field")
        
    if sol_ode is not None:
        ode_start_proxy = Line2D([0], [0], color='tab:orange', marker='o', markersize=6, linestyle='None')
        handles.extend([ode_start_proxy])
        labels.extend(["Initial state (ODE)"])

    ax.legend(handles=handles, labels=labels, loc='best', fontsize='small')
    ax.set_title(title)

    if savefig is not None:
        plt.savefig(f"plots/{savefig}", bbox_inches="tight")
    
    plt.close(fig)

def sdevelo_multiplot(n_runs, dt, steps, sol_ode, title="", savefig=None):
    """Run multiple simulations of the SDEVelo model, calculate mean and confidence intervals."""
    t = np.linspace(0, steps * dt, steps)
    all_results = [solve_sde_velo(dt, steps) for _ in range(n_runs)]
    protein_runs = np.array([res[2] for res in all_results])
    mean_p = np.mean(protein_runs, axis=0)
    err_p = 1.96 * (np.std(protein_runs, axis=0) / np.sqrt(n_runs))

    plt.figure(figsize=(8, 6))
    plt.plot(t, mean_p[0], color="tab:blue")
    plt.fill_between(
        t, mean_p[0] - err_p[0], mean_p[0] + err_p[0], color="tab:blue", alpha=0.2
    )
    plt.plot(t, sol_ode.y[2], color="tab:blue", linestyle="--")

    plt.plot(t, mean_p[1], color="tab:red")
    plt.fill_between(
        t, mean_p[1] - err_p[1], mean_p[1] + err_p[1], color="tab:red", alpha=0.2
    )
    plt.plot(t, sol_ode.y[3], color="tab:red", linestyle="--")

    species_handles = [
        Line2D([0], [0], color="tab:blue", lw=2, label="Protein A"),
        Line2D([0], [0], color="tab:red", lw=2, label="Protein B"),
    ]
    model_handles = [
        Line2D([0], [0], color="black", lw=2, linestyle="-", label="SDEVelo"),
        Line2D([0], [0], color="black", lw=2, linestyle="--", label="ODE"),
    ]

    ax = plt.gca()

    plt.subplots_adjust(right=0.78)

    leg_species = ax.legend(
        handles=species_handles,
        title="Species",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    ax.add_artist(leg_species)

    ax.legend(
        handles=model_handles,
        title="Model",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.80),
        borderaxespad=0.0,
    )
    plt.title(f"{title} ($n_{{\\text{{SDEVelo}}}}={n_runs}$)")
    plt.xlabel("Time [s]")
    plt.ylabel("Concentration [M]")
    plt.grid(True, alpha=0.3)

    if savefig is not None:
        plt.savefig(f"plots/{savefig}")
    plt.clf()


def plot_sdevelo_concentrations(
    n_runs: int, dt: float, steps: int, title: str = "", savefig=None
):
    """Run multiple simulations of the SDEVelo model, calculate mean and confidence intervals for U, S, P, and plot them.
    :param n_runs: Number of independent SDEVelo simulations to run.
    :param dt: Time step for each simulation.
    :param steps: Number of time steps to simulate.
    :param title: Title for the plot.
    :param savefig: Filename to save the plot, if not None.
    """
    all_u = np.zeros((n_runs, 2, steps))
    all_s = np.zeros((n_runs, 2, steps))
    all_p = np.zeros((n_runs, 2, steps))

    for i in range(n_runs):
        u, s, p, t = solve_sde_velo(dt, steps)
        all_u[i] = u
        all_s[i] = s
        all_p[i] = p

    # 2. Calculate Means and Standard Errors
    mean_u = np.mean(all_u, axis=0)
    mean_s = np.mean(all_s, axis=0)
    mean_p = np.mean(all_p, axis=0)

    # 95% Confidence Interval (1.96 * Standard Error)
    err_u = 1.96 * (np.std(all_u, axis=0) / np.sqrt(n_runs))
    err_s = 1.96 * (np.std(all_s, axis=0) / np.sqrt(n_runs))
    err_p = 1.96 * (np.std(all_p, axis=0) / np.sqrt(n_runs))

    # 3. Plotting
    fig, ax = plt.subplots(figsize=(8, 6))

    # Define styling: Colors for genes, linestyles for molecule types
    styles = [
        {"mean": mean_u, "err": err_u, "ls": ":", "alpha_line": 0.8, "alpha_fill": 0.1},
        {
            "mean": mean_s,
            "err": err_s,
            "ls": "--",
            "alpha_line": 0.8,
            "alpha_fill": 0.1,
        },
        {"mean": mean_p, "err": err_p, "ls": "-", "alpha_line": 1.0, "alpha_fill": 0.2},
    ]

    colors = ["tab:blue", "tab:red"]  # Gene A, Gene B

    # Loop through genes and styles to plot lines and fill_betweens
    for gene_idx, color in enumerate(colors):
        for style in styles:
            # Plot mean line
            ax.plot(
                t,
                style["mean"][gene_idx],
                color=color,
                linestyle=style["ls"],
                alpha=style["alpha_line"],
            )
            # Plot shaded error region
            ax.fill_between(
                t,
                style["mean"][gene_idx] - style["err"][gene_idx],
                style["mean"][gene_idx] + style["err"][gene_idx],
                color=color,
                alpha=style["alpha_fill"],
            )

    # 4. Create Custom Legends
    gene_handles = [
        Line2D([0], [0], color="tab:blue", lw=2, label="Gene A"),
        Line2D([0], [0], color="tab:red", lw=2, label="Gene B"),
    ]

    molecule_handles = [
        Line2D([0], [0], color="black", lw=2, linestyle=":", label=r"Unspliced ($U$)"),
        Line2D([0], [0], color="black", lw=2, linestyle="--", label=r"Spliced ($S$)"),
        Line2D([0], [0], color="black", lw=2, linestyle="-", label=r"Protein ($P$)"),
    ]

    # Adjust layout and add legends
    # plt.subplots_adjust(right=0.75)

    leg_genes = ax.legend(
        handles=gene_handles,
        title="Gene",
        loc="upper center",
        # bbox_to_anchor=(1.02, 1.00),
        # borderaxespad=0.0,
    )
    ax.add_artist(leg_genes)

    ax.legend(
        handles=molecule_handles,
        title="Molecule Type",
        loc="upper left",
        # bbox_to_anchor=(1.02, 0.75),
        # borderaxespad=0.0,
    )

    plt.xlabel("Time [s]")
    plt.ylabel("Concentration [M]")
    plt.title(rf"{title} ($n={n_runs})$")
    plt.grid(True, alpha=0.3)

    if savefig is not None:
        plt.savefig(f"plots/{savefig}")
    plt.clf()

def phase_portrait(system,
                    x0=None,
                    xlim=None,
                    ylim=None,
                    dx=5,
                    dy=5,
                    grid_points=100,
                    figsize=(6, 6),
                    fontsize=12,
                    streamplot_kwargs=None,
                    plot_nullclines=True,
                    title=None,
                    xlabel=None,
                    ylabel=None,
                    nullcline_labels=None,
                    fixedpoint_label='stable fixed point',
                    show_legend=True):
    """Create a phase portrait for a 2D system.

    Parameters
    ----------
    system : callable
        Function that accepts a length-2 array-like [x, y] (or two arrays) and
        returns [dx/dt, dy/dt]. The function should support numpy arrays.
    x0 : array-like, optional
        Initial guess for finding a fixed point with fsolve. Default [0, 0].
    xlim, ylim : tuple, optional
        Axis limits as (min, max). If None they are set relative to the fixed
        point using dx/dy.
    dx, dy : float
        Range around the fixed point to build the plotting window when xlim/
        ylim are not provided.
    grid_points : int
        Number of points along each axis for the vector field grid.
    figsize : tuple
        Figure size passed to plt.subplots.
    fontsize : int
        Font size for labels and title.
    streamplot_kwargs : dict, optional
        Additional keyword args forwarded to ax.streamplot.
    plot_nullclines : bool
        If True, draw contour lines where dx/dt=0 and dy/dt=0.

    Returns
    -------
    fig, ax
    """
    # Prepare default kwargs
    if streamplot_kwargs is None:
        streamplot_kwargs = {'density': 1.0, 'color': 'gray'}

    # Find fixed point
    x0 = np.asarray(x0) if x0 is not None else np.array([0.0, 0.0])
    try:
        fp = fsolve(system, x0)
    except Exception:
        # fallback to zero if fsolve fails
        fp = np.array([0.0, 0.0])

    # Determine plotting ranges
    if xlim is None:
        x_min, x_max = fp[0] - dx, fp[0] + dx
    else:
        x_min, x_max = xlim
    if ylim is None:
        y_min, y_max = fp[1] - dy, fp[1] + dy
    else:
        y_min, y_max = ylim

    x = np.linspace(x_min, x_max, int(grid_points))
    y = np.linspace(y_min, y_max, int(grid_points))
    X, Y = np.meshgrid(x, y)

    # Evaluate vector field. Accept systems that return lists or arrays.
    d = system([X, Y])
    # ensure we have two arrays
    dX = np.asarray(d[0])
    dY = np.asarray(d[1])

    fig, ax = plt.subplots(figsize=figsize)

    # Streamplot
    ax.streamplot(X, Y, dX, dY, **streamplot_kwargs)

    # Plot fixed point
    ax.plot(fp[0], fp[1], marker='x', color='k', markersize=8, markeredgewidth=2)
    if fixedpoint_label:
        ax.annotate(fixedpoint_label, xy=(fp[0], fp[1]), xytext=(5, 5), textcoords='offset points', fontsize=fontsize)

    # Nullclines via contour at level 0
    handles = []
    labels = []
    if plot_nullclines:
        try:
            c1 = ax.contour(X, Y, dX, levels=[0], colors='C0', linestyles='--', linewidths=1)
            c2 = ax.contour(X, Y, dY, levels=[0], colors='C1', linestyles='-.', linewidths=1)
            # create legend handles
            # use provided labels if any
            if nullcline_labels is None:
                nc_labels = ('dx/dt = 0', 'dy/dt = 0')
            else:
                # ensure we have at least two labels
                nc_labels = tuple(nullcline_labels)
                if len(nc_labels) < 2:
                    defaults = ('dx/dt = 0', 'dy/dt = 0')
                    nc_labels = tuple(list(nc_labels) + list(defaults[len(nc_labels):]))
            handles.append(Line2D([0], [0], color='C0', linestyle='--'))
            labels.append(nc_labels[0])
            handles.append(Line2D([0], [0], color='C1', linestyle='-.'))
            labels.append(nc_labels[1])
        except Exception:
            # if contouring fails (e.g., non-finite values), skip nullclines
            pass

    # axis labels and formatting
    ax.set_xlabel(xlabel if xlabel is not None else 'x', fontsize=fontsize)
    ax.set_ylabel(ylabel if ylabel is not None else 'y', fontsize=fontsize)
    if title is not None:
        ax.set_title(title, fontsize=fontsize)
    else:
        ax.set_title('Phase portrait', fontsize=fontsize)

    # Add a legend combining nullcline proxies (if any)
    if show_legend and handles:
        # include fixed point in legend too if requested
        fp_label = fixedpoint_label if fixedpoint_label else None
        if fp_label:
            handles = [Line2D([0], [0], color='k', marker='x', linestyle='')] + handles
            labels = [fp_label] + labels
        ax.legend(handles, labels, fontsize=fontsize)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    return fig, ax


if __name__ == "__main__":
    np.random.seed(3)

    os.makedirs("plots/", exist_ok=True)

    ### Task 1: Patient alpha ###

    # (a) Identify the Regulatory Mechanism: Analyze the sequence AGCGC to infer which regulatory
    # mechanism (I or II) is active in this cell line.
    # * Hint: The classification of the total RNA-sequence data should be carried out using the Viterbi
    # algorithm. Use the HMM parameters specified in Tables 1, 2, and 3.

    sequence = "AGCGC"

    B = np.array(
        [[0.25, 0.25, 0.25, 0.25], [0.4, 0.4, 0.05, 0.15]]
    )  # Emission probabilities for the HMM states

    A = np.array(
        [[0.9, 0.1], [0.2, 0.8]]
    )  # Transition probabilities for the HMM states

    P = [0.5, 0.5]  # Initial state probabilities

    path, V, backpointer = viterbi_algorithm(sequence, A, B, P)
    print("For patient Alpha:")
    print(f"Most likely sequence of hidden states for {sequence}:", path)
    print("Viterbi Matrix (probabilities):\n", pd.DataFrame(V, index=["Exon", "Intron"], columns=[f"t={i}" for i in range(len(sequence))]).map(lambda x: f"{x:.6f}"))
    print("Backpointer Matrix (state indices):\n", pd.DataFrame(backpointer, index=["Exon", "Intron"], columns=[f"t={i}" for i in range(len(sequence))]))
    print("\n")

    # (b) Based on the mechanism you identified for patient alpha, choose the appropriate mathematical
    # framework to model the gene regulation dynamics ODE model vs SDEVelo model. Also, the choice of
    # model formulation may depend on the observed mRNA status. Parameter values for these models are
    # given in Tables 4 and 5, respectively. Ensure to present relevant graphs/visualizations.

    # Trascription problem: use ODE (assumption: simple two-gene network)
    m_a, m_b = 2.35, 2.35  # unit
    gamma_a, gamma_b = 1, 1  # s-1
    k_pa, k_pb = 1.0, 1.0  # s-1
    theta_a, theta_b = 0.21, 0.21  # M
    n_a, n_b = (
        3,
        3,
    )  # hill coefficient (assumption: use hill and not PWL activation because of low n coefficients)
    delta_pa, delta_pb = 1.0, 1.0  # s-1
    mRNA_a, mRNA_b = 0.8, 0.8  # M
    P_a, P_b = 0.8, 0.8  # M

    t_eval = np.linspace(0, 100, 1000)
    y0 = [mRNA_a, mRNA_b, P_a, P_b]

    sol_alpha = solve_ivp(patient_alpha_ode, [0, 100], y0, t_eval=t_eval, args=(True,))
    sol_alpha_regular = solve_ivp(
        patient_alpha_ode, [0, 100], y0, t_eval=t_eval, args=(False,)
    )

    # regular
    concentration_overtime(sol_alpha_regular, t_eval, hijack=False)
    phase_plots(sol_alpha_regular, t_eval, hijack=False)
    hill_plots(sol_alpha_regular, t_eval, hijack=False)
    plot_rates(sol_alpha_regular, t_eval, hijack=False)

    # hijack 
    concentration_overtime(sol_alpha, t_eval, hijack=True)
    phase_plots(sol_alpha, t_eval, hijack=True)
    hill_plots(sol_alpha, t_eval, hijack=True)
    plot_rates(sol_alpha, t_eval, hijack=True)

    ## patient beta

    ### Task 2: Analysis of Patient Sample Beta ###

    # Next analyze “Patient Beta”, a cell line derived from a highly aggressive and fast-growing tumor. The
    # RNA-seq data for the PROLIFERATOR gene in this sample is dominated by the sequence AUUAU.
    # Repeat the full analysis pipeline (a) and (b) from Patient Alpha for this new sample.

    sequence = "AUUAU"

    path, V, backpointer = viterbi_algorithm(sequence, A, B, P)
    print("For patient Beta:")
    print(f"Most likely sequence of hidden states for {sequence}:", path)
    print("Viterbi Matrix (probabilities):\n", pd.DataFrame(V, index=["Exon", "Intron"], columns=[f"t={i}" for i in range(len(sequence))]).map(lambda x: f"{x:.6f}"))
    print("Backpointer Matrix (state indices):\n", pd.DataFrame(backpointer, index=["Exon", "Intron"], columns=[f"t={i}" for i in range(len(sequence))]))

    a = np.array([1.0, 0.25])
    b = np.array([0.0005, 0.0005])
    c = np.array([2.0, 0.5])
    beta = np.array([2.35, 2.35])
    gamma = np.array([1.0, 1.0])
    n = np.array([3, 3])
    theta = np.array([0.21, 0.21])
    k_P = np.array([1.0, 1.0])
    m = np.array([2.35, 2.35])
    delta = np.array([1.0, 1.0])
    sigma_1 = np.array([0.05, 0.05])
    sigma_2 = np.array([0.05, 0.05])

    dt, t_max = 0.01, 100
    steps = int(t_max / dt)
    t_eval = np.linspace(0, t_max, steps)
    sol_ode = solve_ivp(patient_alpha_ode, [0, t_max], y0, t_eval=t_eval, args=(False,))

    sdevelo_phase_portrait(
        solve_sde_velo(dt, steps),
        "Patient Beta (SDE)",
        # ode_func=patient_alpha_ode,
        sol_ode=sol_ode,
        label_ode="Healthy State (ODE)",
        # ode_args=(False,),
        grid_size=30,
        title="Phase Portrait of Protein concentrations",
        savefig="sdevelo_phase_portrait",
    )

    # hijacked_ode = solve_ivp(patient_alpha_ode, [0, t_max], y0, t_eval=t_eval, args=(True,))
    # sdevelo_multiplot(
    #     10,
    #     dt,
    #     steps,
    #     hijacked_ode,
    #     title="Protein concentration over time",
    #     savefig="sdevelo_versus_ode",
    # )

    plot_sdevelo_concentrations(
        10,
        dt,
        steps,
        title="Concentration of (un)spliced RNA and protein over time",
        savefig="sdevelo_concentrations_over_time",
    )

    ### Task 3: Comparative Analysis & Diagnosis ###

    # Synthesize your findings from Patient Alpha and Patient Beta to deliver a diagnosis.
    # Compare the protein phase portraits for Patient Alpha and Patient Beta. Describe the dynamic
    # behaviors observed in each plot. What do these different behaviors imply about the cellular "fate" of
    # each sample? Justify your conclusions based on the trajectories in your plots.


    # BONUS
    # a) Explain what each term in the equations (αR, −βRE, −γE, δRE) represents biologically in the 
    # context of a cellular resource and a growth-promoting enzyme. (In your Model/Eqn video) 
    # b) Carry out a stability analysis  for, 𝛼 = 2, 𝛽 = 1.1, 𝛾 = 1 and 𝛿 = 0.9, R(0) = 1 and E(0) =0.5 
    # (And discuss in your results and discussions video) 
    # c) Plot the stream plot with equilibrium point and nullclines (And explain in your Figures video) 

    alpha, beta, gamma, delta = 2, 1.1, 1, 0.9
    a0, b0 = 1, 0.5
    def dadt(a, b): # a is R, b is E
        return alpha*a - beta*a*b 

    def dbdt(a, b):
        return -gamma*b + delta*a*b 

    def nca(alpha, beta):
        return alpha/beta

    def ncb(gamma, delta):
        return gamma/delta

    def system(vars):
        a, b = vars[0], vars[1]
        return [dadt(a, b), dbdt(a, b)]
    
    fig, ax = phase_portrait(system,
                             x0=[a0, b0],
                             dx=5,
                             dy=5,
                             grid_points=200,
                             figsize=(6, 6),
                             fontsize=12,
                             streamplot_kwargs={'density': 1.5, 'color': 'gray'},
                             xlabel="R",
                             ylabel="E")
    plt.savefig("plots/bonus_phase_portrait.png", dpi=300)