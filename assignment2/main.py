import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.lines import Line2D


np.random.seed(3)

### Functions for Task 1, 2, 3 ###


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

    return best_path


def rna_vel(t, u, s, beta, gamma, c, a, b):
    """Calculate the RNA velocity based on the given parameters.
    dudt = alpha - b * u"""

    alpha = c / (1 + np.exp(b * (t - a)))

    dudt = alpha - beta * u

    dsdt = beta * u - gamma * s


### Task 1: Patient alpha ###

# (a) Identify the Regulatory Mechanism: Analyze the sequence AGCGC to infer which regulatory
# mechanism (I or II) is active in this cell line.
# * Hint: The classification of the total RNA-sequence data should be carried out using the Viterbi
# algorithm. Use the HMM parameters specified in Tables 1, 2, and 3.

sequence = "AGCGC"

B = np.array(
    [[0.25, 0.25, 0.25, 0.25], [0.4, 0.4, 0.05, 0.15]]
)  # Emission probabilities for the HMM states

A = np.array([[0.9, 0.1], [0.2, 0.8]])  # Transition probabilities for the HMM states

P = [0.5, 0.5]  # Initial state probabilities


path = viterbi_algorithm(sequence, A, B, P)
print(f"Most likely sequence of hidden states for {sequence}:", path)

# (b) Based on the mechanism you identified for patient alpha, choose the appropriate mathematical
# framework to model the gene regulation dynamics ODE model vs SDEVelo model. Also, the choice of
# model formulation may depend on the observed mRNA status. Parameter values for these models are
# given in Tables 4 and 5, respectively. Ensure to present relevant graphs/visualizations.

# Trascription problem: use ODE (assumption: simple two-gene network)
m_a, m_b = 2.35, 2.35 #s-1 # assumption: no steady state since same rates
gamma_a, gamma_b = 1,1  #s-1
k_pa, k_pb = 1.0, 1.0 #s-1 
theta_a, theta_b = 0.21, 0.21 #M
n_a, n_b = 3, 3 # hill coefficient (assumption: use hill and not PWL activation because of low n coefficients)
delta_pa, delta_pb = 1.0, 1.0 #s-1
mRNA_a, mRNA_b = 0.8, 0.8 #M
P_a, P_b = 0.8, 0.8 #M


def hill_activation(P, theta, n):
    """Calculate the Hill activation function."""
    return P**n / (theta**n + P**n)

def hill_inhibition(P, theta, n):
    """Calculate the Hill inhibition function."""
    return 1- hill_activation(P, theta, n)

def patient_alpha_ode(t, y, hijack = None):
    """Calculate the derivatives of mRNA and protein for Patient Alpha."""
    ra, rb, P_a, P_b = y
    dra_dt = m_a * hill_activation(P_b, theta_b, n_b) - gamma_a * ra
    if hijack:
        drb_dt = m_b * 1 - gamma_b * rb # no inhibition on rb
    else:
        drb_dt = m_b * hill_inhibition(P_a, theta_a, n_a) - gamma_b * rb
    dPa_dt = k_pa * ra - delta_pa * P_a
    dPb_dt = k_pb * rb - delta_pb * P_b
    return [dra_dt, drb_dt, dPa_dt, dPb_dt]

t_eval = np.linspace(0, 100, 1000)
y0 = [mRNA_a, mRNA_b, P_a, P_b]

sol_alpha = solve_ivp(patient_alpha_ode, [0, 100], y0, t_eval=t_eval, args=(True,))
sol_alpha_regular = solve_ivp(patient_alpha_ode, [0, 100], y0, t_eval=t_eval, args=(False,))

def plots_alpha(sol, t, hijack=None):
    ra, rb, pA, pB = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    
    # pA, pB over time
    plt.figure(figsize=(8, 6))
    plt.plot(t, pA, label='Protein A (Guardian)', linewidth=2)
    plt.plot(t, pB, label='Protein B (Proliferator)', linewidth=2)
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Concentration [M]', fontsize=12)
    plt.title('Time Evolution of Protein Concentrations', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_1_proteins_time_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_1_proteins_time.png", dpi=300)
    plt.show()
    
    # mRNAs over time
    plt.figure(figsize=(8, 6))
    plt.plot(t, ra, label='mRNA A',linewidth=2, linestyle='--')
    plt.plot(t, rb, label='mRNA B', linewidth=2, linestyle='--')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Concentration [M]', fontsize=12)
    plt.title('Time Evolution of mRNA Concentrations', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_2_mRNAs_time_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_2_mRNAs_time.png", dpi=300)
    plt.show()
    
    # phase plot pA vs pB and mRNA_A vs mRNA_B
    plt.figure(figsize=(8, 6))
    plt.plot(pA, pB, color='black', linewidth=2, label='Protein Trajectory')
    plt.scatter([pA[0]], [pB[0]], color='black', zorder=5, label='Start Proteins (t=0)')
    plt.scatter([pA[-1]], [pB[-1]], color='yellow', zorder=5, s=100, marker='*', label='End (Proteins Equilibrium)')
    plt.plot(ra, rb, color='blue', linewidth=2, label='mRNA Trajectory')
    plt.scatter([ra[0]], [rb[0]], color='blue', zorder=5, label='Start mRNA (t=0)')
    plt.scatter([ra[-1]], [rb[-1]], color='orange', zorder=5, s=100, marker='*', label='End (mRNA Equilibrium)')
    plt.xlabel('$P_a, r_a$ Concentration [M]', fontsize=12)
    plt.ylabel('$P_b, r_b$ Concentration [M]', fontsize=12)
    plt.title('Phase Space: Protein A vs Protein B and mRNA A vs mRNA B', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_3_phase_space_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_3_phase_space.png", dpi=300)
    plt.show()
    
    # hill Functions
    P_vals = np.linspace(0, 1.5, 100)
    act = hill_activation(P_vals, theta_a, n_a)
    inh = hill_inhibition(P_vals, theta_a, n_a)
    
    plt.figure(figsize=(8, 6))
    plt.plot(P_vals, act, label=f'Activation ($\\theta$={theta_a}, n={n_a})', color='green', linewidth=2)
    if hijack: 
        plt.plot(P_vals, np.ones_like(P_vals), label=f'Hijacked, no inhibition', color='orange', linewidth=2)
    else:
        plt.plot(P_vals, inh, label=f'Inhibition ($\\theta$={theta_a}, n={n_a})', color='orange', linewidth=2)
    plt.axvline(x=theta_a, color='grey', linestyle=':', label='Threshold ($\\theta$)')
    plt.xlabel('Regulator Protein Concentration [M]', fontsize=12)
    plt.ylabel('Hill Function', fontsize=12)
    plt.title('Hill Functions for Activation & Inhibition', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_4_hill_functions_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_4_hill_functions.png", dpi=300)
    plt.show()
    
    # transcription rates over time (m_a, m_b, hill functions)
    transcription_rate_A = m_a * hill_activation(pB, theta_b, n_b)
    if hijack:
        transcription_rate_B = np.full_like(t, m_b) # 1*m_b, no inhibition on B
    else:
        transcription_rate_B = m_b * hill_inhibition(pA, theta_a, n_a)
    
    plt.figure(figsize=(8, 6))
    plt.plot(t, transcription_rate_A, label='Gene A Transcription Rate $m_a * H_{\\theta_b, n_b}(P_B)$', linewidth=2)
    plt.plot(t, transcription_rate_B, label='Gene B Transcription Rate $m_b * H_{\\theta_a, n_a}(P_A)$', linewidth=2)
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Rate [M/s]', fontsize=12)
    plt.title('Transcription Rates Over Time', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_5_transcription_rates_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_5_transcription_rates.png", dpi=300)
    plt.show()
    
    # translation rates over time (k_pa, k_pb, mRNA concentrations)
    translation_rate_A = k_pa * ra
    translation_rate_B = k_pb * rb
    
    plt.figure(figsize=(8, 6))
    plt.plot(t, translation_rate_A, label='Protein A Translation Rate $k_{pa} * r_a$', linewidth=2)
    plt.plot(t, translation_rate_B, label='Protein B Translation Rate $k_{pb} * r_b$', linewidth=2)
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Rate [M/s]', fontsize=12)
    plt.title('Translation Rates Over Time', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if hijack:
        plt.savefig("alpha_6_translation_rates_hijack.png", dpi=300)
    else:
        plt.savefig("alpha_6_translation_rates.png", dpi=300)
    plt.show()

plots_alpha(sol_alpha, t_eval, hijack=True)
plots_alpha(sol_alpha_regular, t_eval, hijack=False)

### Task 2: Analysis of Patient Sample Beta ###

# Next analyze “Patient Beta”, a cell line derived from a highly aggressive and fast-growing tumor. The
# RNA-seq data for the PROLIFERATOR gene in this sample is dominated by the sequence AUUAU.
# Repeat the full analysis pipeline (a) and (b) from Patient Alpha for this new sample.

sequence = "AUUAU"

path = viterbi_algorithm(sequence, A, B, P)
print(f"Most likely sequence of hidden states for {sequence}:", path)

a = np.array([1.0, 0.25])
b = np.array([0.0005, 0.0005])
c = np.array([2.0, 0.5])
beta = np.array([2.35, 2.35])
gamma = np.array([1.0, 1.0])
n = np.array([3, 3])
theta = np.array([0.21, 0.21])
k_p = np.array([1.0, 1.0])
m = np.array([2.35, 2.35])
delta = np.array([1.0, 1.0])
sigma_1 = np.array([0.05, 0.05])
sigma_2 = np.array([0.05, 0.05])


def solve_sde_velo(dt, steps):
    t = np.linspace(0, steps * dt, steps)

    # empty arrays:
    u = np.zeros((2, steps))
    s = np.zeros((2, steps))
    p = np.zeros((2, steps))

    # initial concentrations:
    u[:, 0] = 0.8
    s[:, 0] = 0.8
    p[:, 0] = 0.8

    for i in range(steps - 1):
        inhibition_B = hill_inhibition(p[0, i], theta[1], n[1])
        current_beta = np.array([beta[0], beta[1] * inhibition_B])

        for g in range(2):
            dW1 = np.random.normal(0, np.sqrt(dt))
            dW2 = np.random.normal(0, np.sqrt(dt))

            du = (c[g] - current_beta[g] * u[g, i]) * dt + sigma_1[g] * dW1
            u[g, i + 1] = max(0, u[g, i] + du)

            ds = (current_beta[g] * u[g, i] - gamma[g] * s[g, i]) * dt + sigma_2[g] * dW2
            s[g, i + 1] = max(0, s[g, i] + ds)

            dp = (k_p[g] * s[g, i] - delta[g] * p[g, i]) * dt
            p[g, i + 1] = max(0, p[g, i] + dp)

    return u, s, p, t

ra, rb, pA, pB = sol_alpha.y[0], sol_alpha.y[1], sol_alpha.y[2], sol_alpha.y[3]
u, s, p, t = solve_sde_velo(0.01, 1000)

plt.plot(p[0], p[1], label="Patient Beta", color="tab:red")
plt.axvline(theta[0], color="gray", linestyle="--", label="Binding Threshold (θ)")
plt.axhline(theta[1], color="gray", linestyle="--")
plt.plot(sol_alpha.y[2], sol_alpha.y[3], label="Patient Alpha", color="tab:blue")

plt.xlabel("Protein A Concentration [M]")
plt.ylabel("Protein B Concentration [M]")
plt.legend()
plt.grid(True, which="both", linestyle=":", alpha=0.5)

plt.scatter(p[0, -1], p[1, -1], color="black", zorder=5)
plt.annotate(
    "Final State", (p[0, -1], p[1, -1]), textcoords="offset points", xytext=(10, -10)
)

plt.savefig("A_vs_B_Phase")
plt.show()


n_runs = 10
all_results = [solve_sde_velo(0.01, 1000) for _ in range(n_runs)]
protein_runs = np.array([res[2] for res in all_results])
mean_p = np.mean(protein_runs, axis=0)
err_p = 1.96 * (np.std(protein_runs, axis=0) / np.sqrt(n_runs))

plt.plot(t, mean_p[0], color="tab:blue")
plt.fill_between(
    t, mean_p[0] - err_p[0], mean_p[0] + err_p[0], color="tab:blue", alpha=0.2
)
plt.plot(t, sol_alpha.y[2], color="tab:blue", linestyle="--")

plt.plot(t, mean_p[1], color="tab:red")
plt.fill_between(
    t, mean_p[1] - err_p[1], mean_p[1] + err_p[1], color="tab:red", alpha=0.2
)
plt.plot(t, sol_alpha.y[3], color="tab:red", linestyle="--")


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
    handles=species_handles, title="Species",
    loc="upper left", bbox_to_anchor=(1.02, 1.00), borderaxespad=0.0
)
ax.add_artist(leg_species)

ax.legend(
    handles=model_handles, title="Model",
    loc="upper left", bbox_to_anchor=(1.02, 0.80), borderaxespad=0.0
)


plt.xlabel("Time [s]")
plt.ylabel("Concentration [M]")

plt.savefig("ODE_vs_SDE")
plt.show()

### Task 3: Comparative Analysis & Diagnosis ###

# Synthesize your findings from Patient Alpha and Patient Beta to deliver a diagnosis.
# Compare the protein phase portraits for Patient Alpha and Patient Beta. Describe the dynamic
# behaviors observed in each plot. What do these different behaviors imply about the cellular "fate" of
# each sample? Justify your conclusions based on the trajectories in your plots.

# %%
