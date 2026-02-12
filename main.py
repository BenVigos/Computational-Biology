import matplotlib.pyplot as plt
import numpy as np

def main():
    print("Hello from computational-biology!")


if __name__ == "__main__":
    main()

    # Transition matrix A between letters [I, N, T, U].
    A = [[0.1, 0.4, 0.2, 0.3], [0.1, 0.1, 0.6, 0.2], [0.25, 0.25, 0.25, 0.25], [0.3, 0.3, 0.1, 0.3]]

    # Initial selection probabilities for [I, N, T, U]
    p_init = [0.3, 0.3, 0, 0.4]

    import random
    from typing import List, Tuple

    LETTERS = ['I', 'N', 'T', 'U']

    def _sample_index(dist: List[float]) -> int:
        """Sample an index from a (possibly non-normalized) discrete distribution."""
        total = sum(dist)
        if total <= 0:
            raise ValueError("Distribution must have positive total probability")
        r = random.random() * total
        cum = 0.0
        for i, p in enumerate(dist):
            cum += p
            if r <= cum:
                return i
        return len(dist) - 1

    def generate_word(A: List[List[float]], p_init: List[float], letters: List[str]) -> str:
        """Generate a single 4-letter word from the Markov chain.

        Starts by sampling the initial letter from p_init, then performs three transitions
        using rows of A (A[current_index] gives next-state distribution).
        """
        idx = _sample_index(p_init)
        word_chars = [letters[idx]]
        for _ in range(3):
            row = A[idx]
            idx = _sample_index(row)
            word_chars.append(letters[idx])
        return ''.join(word_chars)

    def find_until_target(A: List[List[float]], p_init: List[float], target: str, letters: List[str], max_iters: int = 1_000_000) -> Tuple[int, str]:
        """Keep generating 4-letter words until `target` is generated or max_iters reached.

        Returns (attempts, word) when target found. Raises RuntimeError on exhaustion.
        """
        for attempt in range(1, max_iters + 1):
            w = generate_word(A, p_init, letters)
            if w == target:
                return attempt, w
        raise RuntimeError(f"Target '{target}' not found in {max_iters} attempts")

    # Run the search for the target word 'UNIT'

    att = []

    for i in range(1000):

        attempts, found = find_until_target(A, p_init, 'UNIT', LETTERS)
        att.append(attempts)

    print(f"Found target word '{found}' after {attempts} generated word(s).")
    plt.hist(att)
    plt.show()
    print(f"Average attempts: {np.mean(att):.2f}, Std Dev: {np.std(att):.2f}")