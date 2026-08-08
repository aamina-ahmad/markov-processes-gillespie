"""Run the submitted two-state Gillespie experiment and save the comparison plot."""

from pathlib import Path

import matplotlib.pyplot as plt

from src.gillespie import analytical_solution, estimate_probabilities


ALPHA = 0.3
BETA = 0.5
T_END = 15
TRIALS = 10000
N_BINS = 300
STARTING_STATE = 1


def main():
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    time_grid, p1_est, p2_est = estimate_probabilities(
        ALPHA,
        BETA,
        T_END,
        TRIALS,
        N_BINS,
        starting_state=STARTING_STATE,
    )

    t, p1_exact, p2_exact = analytical_solution(
        ALPHA,
        BETA,
        T_END,
        N_BINS,
        p10=1.0 if STARTING_STATE == 1 else 0.0,
    )

    pi_1 = BETA / (ALPHA + BETA)
    pi_2 = ALPHA / (ALPHA + BETA)

    plt.figure(figsize=(9, 5))
    plt.plot(time_grid, p1_est, label="Simulated p1(t)")
    plt.plot(time_grid, p2_est, label="Simulated p2(t)")
    plt.plot(t, p1_exact, "--", label="Analytical p1(t)")
    plt.plot(t, p2_exact, "--", label="Analytical p2(t)")
    plt.axhline(pi_1, linestyle=":", label=f"pi1 = {pi_1:.4f}")
    plt.axhline(pi_2, linestyle=":", label=f"pi2 = {pi_2:.4f}")
    plt.xlabel("Time")
    plt.ylabel("Probability")
    plt.title(
        "Gillespie Simulation vs Analytical Solution "
        f"(alpha={ALPHA}, beta={BETA})"
    )
    plt.legend()
    plt.tight_layout()

    output = results_dir / "gillespie_vs_analytical.png"
    plt.savefig(output, dpi=160)
    plt.close()

    print(f"Saved {output}")
    print(f"Analytical stationary probabilities: pi1={pi_1:.4f}, pi2={pi_2:.4f}")
    print(
        "Final simulated probabilities: "
        f"p1={p1_est[-1]:.4f}, p2={p2_est[-1]:.4f}"
    )


if __name__ == "__main__":
    main()
