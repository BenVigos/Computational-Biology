import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.lines import Line2D


np.random.seed(3)

### Functions for Task 1, 2, 3 ###


def get_letter_index(letter):
    """
    Helper function to convert a letter (A, C, G, T) to an index (0, 1, 2, 3).
    :param letter: A character representing a nucleotide.
    :return: An integer index corresponding to the nucleotide.
    """
    index = {"A": 0, "U": 1, "G": 2, "C": 3}
    return index.get(letter)


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

# Trascription problem: use ODE
m_a, m_b = 2.35, 2.35  # s-1
gamma_a, gamma_b = 1, 1  # s-1
k_pa, k_pb = 1.0, 1.0  # s-1
theta_a, theta_b = 0.21, 0.21  # M
n_a, n_b = 3, 3  # hill coefficient
delta_pa, delta_pb = 1.0, 1.0  # s-1
mRNA_a, mRNA_b = 0.8, 0.8  # M
P_a, P_b = 0.8, 0.8  # M


def hill_activation(P, theta, n):
    """Calculate the Hill activation function."""
    return P**n / (theta**n + P**n)


def hill_inhibition(P, theta, n):
    """Calculate the Hill inhibition function."""
    return 1 - hill_activation(P, theta, n)


def patient_alpha_ode(t, y):
    """Calculate the derivatives of mRNA and protein for Patient Alpha."""
    ra, rb, P_a, P_b = y
    dra_dt = m_a * hill_inhibition(P_b, theta_b, n_b) - gamma_a * ra
    drb_dt = m_b * 1 - gamma_b * rb  # no inhibition on rb
    dPa_dt = k_pa * ra - delta_pa * P_a
    dPb_dt = k_pb * rb - delta_pb * P_b
    return [dra_dt, drb_dt, dPa_dt, dPb_dt]


t = np.linspace(0, 50, 1000)
y0 = [mRNA_a, mRNA_b, P_a, P_b]
sol_alpha = solve_ivp(patient_alpha_ode, [0, 50], y0, t_eval=t)

pA_alpha = sol_alpha.y[2]
pB_alpha = sol_alpha.y[3]

plt.figure(figsize=(8, 6))
plt.plot(t, pA_alpha, label="pA(t)", color="red", linewidth=2)
plt.xlabel("Tempo [s]", fontsize=12)
plt.ylabel("Concentrazione Proteina A [pA]", fontsize=12)
plt.title("Time Evolution of pA Concentration", fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig("pA over time")
plt.clf()

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


u, s, p, t = solve_sde_velo(0.01, 1000)

plt.plot(p[0], p[1], label="Patient Beta", color="tab:red")
plt.axvline(theta[0], color="gray", linestyle="--", label="Binding Threshold (θ)")
plt.axhline(theta[1], color="gray", linestyle="--")
plt.plot(pA_alpha, pB_alpha, label="Patient Alpha", color="tab:blue")

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
plt.plot(t, pA_alpha, color="tab:blue", linestyle="--")

plt.plot(t, mean_p[1], color="tab:red")
plt.fill_between(
    t, mean_p[1] - err_p[1], mean_p[1] + err_p[1], color="tab:red", alpha=0.2
)
plt.plot(t, pB_alpha, color="tab:red", linestyle="--")


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
