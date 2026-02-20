import numpy as np
import matplotlib.pyplot as plt


### Functions for Task 1, 2, 3 ###

def get_letter_index(letter):
    """
    Helper function to convert a letter (A, C, G, T) to an index (0, 1, 2, 3).
    :param letter: A character representing a nucleotide.
    :return: An integer index corresponding to the nucleotide.
    """
    match letter:
        case 'A':
            return 0
        case 'U':
            return 1
        case 'G':
            return 2
        case 'C':
            return 3

def get_exon_intron_state(state):
    """
    Helper function to convert a state index to a string representation (Exon or Intron).
    :param state: An integer index representing the state.
    :return: A string 'Exon' or 'Intron' corresponding to the state.
    """
    match state:
        case 0:
            return 'E'
        case 1:
            return 'I'
        case _:
            return '?'

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
    backpointer = np.zeros((n_states, T), dtype=int)-1

    # Initialization step
    for s in range(n_states):
        V[s, 0] = P[s] * B[s, get_letter_index(sequence[s])]

    # Recursion step
    for t in range(1, T):
        for s in range(n_states):
            max_prob = -1
            max_state = 0
            for s_prev in range(n_states):
                prob = V[s_prev, t-1] * A[s_prev, s] * B[s, get_letter_index(sequence[s])]
                if prob > max_prob:
                    max_prob = prob
                    max_state = s_prev
            V[s, t] = max_prob
            backpointer[s, t] = max_state

    # Termination step
    best_path_prob = -1
    best_last_state = 0
    for s in range(n_states):
        if V[s, T-1] > best_path_prob:
            best_path_prob = V[s, T-1]
            best_last_state = s

    # Backtrack to find the best path
    best_path = [best_last_state]
    for t in range(T-2, 0, -1):
        best_path.append( backpointer[best_path[-1], t])

    best_path.append(-1)
    best_path.reverse()

    for i, s in enumerate(best_path):
        best_path[i] = get_exon_intron_state(s)

    return best_path

def rna_vel(t, u, s, beta, gamma, c, a, b):
    """Calculate the RNA velocity based on the given parameters.
    dudt = alpha - b * u"""

    alpha = c / (1 + np.exp(b*(t-a)))

    dudt = alpha - beta * u

    dsdt = beta * u - gamma * s



### Task 1: Patient alpha ###

# (a) Identify the Regulatory Mechanism: Analyze the sequence AGCGC to infer which regulatory
# mechanism (I or II) is active in this cell line.
# * Hint: The classification of the total RNA-sequence data should be carried out using the Viterbi
# algorithm. Use the HMM parameters specified in Tables 1, 2, and 3.

sequence = 'AGCGC'

B = np.array([[0.25, 0.25, 0.25, 0.25],
     [0.4, 0.4, 0.05, 0.15]]) # Emission probabilities for the HMM states

A = np.array([[0.9, 0.1],
     [0.2, 0.8]]) # Transition probabilities for the HMM states

P = [0.5, 0.5] # Initial state probabilities


path = viterbi_algorithm(sequence, A, B, P)
print("Most likely sequence of hidden states:", path)






### Task 2: Analysis of Patient Sample Beta ###

# Next analyze “Patient Beta”, a cell line derived from a highly aggressive and fast-growing tumor. The
# RNA-seq data for the PROLIFERATOR gene in this sample is dominated by the sequence AUUAU.
# Repeat the full analysis pipeline (a) and (b) from Patient Alpha for this new sample.

sequence = 'AUUAU'

path = viterbi_algorithm(sequence, A, B, P)
print("Most likely sequence of hidden states:", path)


### Task 3: Comparative Analysis & Diagnosis ###

# Synthesize your findings from Patient Alpha and Patient Beta to deliver a diagnosis.
# Compare the protein phase portraits for Patient Alpha and Patient Beta. Describe the dynamic
# behaviors observed in each plot. What do these different behaviors imply about the cellular "fate" of
# each sample? Justify your conclusions based on the trajectories in your plots.