import numpy as np


# simulates a two-state Markov process using the Gillespie algorithm
# returns (times, states) as lists
def gillespie(alpha, beta, t_end, starting_state=1):
    t = 0
    current_state = starting_state
    times = [0]
    states = [current_state]

    # total rate out of each state
    rates = {1: alpha, 2: beta}

    while t <= t_end:
        # draw waiting time exponentially distributed with rate out of current state
        rate = rates[current_state]
        tau = np.random.exponential(1.0 / rate)
        t += tau

        # only one possible transition from each state, so always jump to the other state
        current_state = 2 if current_state == 1 else 1

        times.append(t)
        states.append(current_state)

    return times, states


# estimates p_1(t) and p_2(t) by averaging over many Gillespie trajectories
# returns (time_grid, p1_estimates, p2_estimates)
def estimate_probabilities(alpha, beta, t_end, trials, n_bins, starting_state=1):
    time_grid = np.linspace(0, t_end, n_bins)
    p1_counts = np.zeros(n_bins)
    p2_counts = np.zeros(n_bins)

    for trial in range(trials):
        np.random.seed(trial)  # seed per trial for reproducibility
        times, states = gillespie(
            alpha, beta, t_end, starting_state=starting_state
        )

        # track state at each time step
        state = starting_state
        swap_index = 0

        for i in range(len(time_grid)):
            t = time_grid[i]

            # process all jumps that occurred before or at t
            while (
                swap_index < len(times) - 1
                and times[swap_index + 1] <= t
            ):
                swap_index += 1
                state = states[swap_index]

            if state == 1:
                p1_counts[i] += 1
            else:
                p2_counts[i] += 1

    # normalise counts to get probability estimates
    p1_estimates = p1_counts / trials
    p2_estimates = p2_counts / trials

    return time_grid, p1_estimates, p2_estimates


# analytical solution from part (e)
# computes p_1(t) and p_2(t)
def analytical_solution(alpha, beta, t_end, n_bins, p10):
    # stationary probabilities (eigenvector of W for eigenvalue 0)
    pi_1 = beta / (alpha + beta)
    pi_2 = alpha / (alpha + beta)

    p20 = 1 - p10  # initial probability of state 2
    t = np.linspace(0, t_end, n_bins)

    # decaying exponential with relaxation rate (alpha + beta)
    p1 = pi_1 + (p10 - pi_1) * np.exp(-(alpha + beta) * t)
    p2 = pi_2 + (p20 - pi_2) * np.exp(-(alpha + beta) * t)

    return t, p1, p2
