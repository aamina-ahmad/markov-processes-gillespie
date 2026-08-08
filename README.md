# Markov Processes and Gillespie Algorithm

Analytical and stochastic study of a two-state continuous-time Markov process, combining master equations, steady-state analysis and Gillespie simulation.

## Overview

The model contains two states with transition rates:

```text
State 1 → State 2: α
State 2 → State 1: β
```

where `α > 0` and `β > 0`.

The state probabilities are `p1(t)` and `p2(t)`, with:

```text
p1(t) + p2(t) = 1
```

The project derives the analytical behaviour of the process and then verifies it using stochastic Gillespie simulations.

---

## Master Equation

The probability dynamics are:

```text
dp1/dt = βp2 - αp1
dp2/dt = αp1 - βp2
```

The generator matrix is:

```text
     [ -α    β ]
W =  [  α   -β ]
```

Its key properties are:

- each column sums to zero, conserving total probability
- off-diagonal entries are non-negative
- one eigenvalue is zero
- the second eigenvalue is negative

The eigenvalues are:

```text
λ1 = 0
λ2 = -(α + β)
```

---

## Stationary Distribution

The stationary distribution satisfies:

```text
Wπ = 0
π1 + π2 = 1
```

which gives:

```text
π1 = β / (α + β)
π2 = α / (α + β)
```

For the experiment:

```text
α = 0.3
β = 0.5
```

so:

```text
π1 = 0.625
π2 = 0.375
```

The system therefore approaches:

```text
π = (0.625, 0.375)
```

as time increases.

---

## Time-Dependent Solution

The analytical probabilities are:

```text
p1(t) = π1 + (p1(0) - π1) exp(-(α + β)t)

p2(t) = π2 + (p2(0) - π2) exp(-(α + β)t)
```

The transient component decays exponentially at rate:

```text
α + β
```

giving the relaxation time:

```text
τ = 1 / (α + β)
```

For `α = 0.3` and `β = 0.5`:

```text
τ = 1.25
```

---

## Detailed Balance and Convergence

The stationary distribution satisfies detailed balance:

```text
π1 α = π2 β
```

The project also uses Kullback-Leibler divergence to show that the distribution converges toward the stationary state.

The KL divergence is non-increasing over time:

```text
dD(p(t) || π) / dt ≤ 0
```

so it acts as a Lyapunov function for the dynamics.

---

## Gillespie Simulation

The Gillespie algorithm simulates individual random trajectories of the continuous-time Markov process.

While the system is in a state, the waiting time before the next transition is sampled from an exponential distribution:

```python
tau = np.random.exponential(1.0 / rate)
```

Because the model contains only two states, each transition moves the system to the opposite state.

The experiment uses:

```text
alpha = 0.3
beta = 0.5
t_end = 15
trajectories = 10000
time bins = 300
starting state = 1
```

The probabilities `p1(t)` and `p2(t)` are estimated by averaging across the 10,000 simulated trajectories.

---

## Results

The analytical stationary probabilities are:

```text
π1 = 0.6250
π2 = 0.3750
```

The original Gillespie experiment produced simulated values of approximately:

```text
p1 = 0.6247
p2 = 0.3753
```

The simulated probability curves closely match the analytical solution and converge to the predicted stationary distribution.

This demonstrates that averaging many stochastic Gillespie trajectories recovers the deterministic probability evolution described by the master equation.

---

## Repository Structure

```text
markov-processes-gillespie/
│
├── README.md
├── requirements.txt
├── .gitignore
├── run_simulation.py
│
├── src/
│   ├── __init__.py
│   └── gillespie.py
│
└── results/
    └── README.md
```

`src/gillespie.py` contains:

- the Gillespie trajectory simulator
- probability estimation across trajectories
- the analytical solution

`run_simulation.py` runs the experiment and generates the comparison plot.

---

## Running the Project

Clone the repository:

```bash
git clone https://github.com/aamina-ahmad/markov-processes-gillespie.git
```

Move into the project:

```bash
cd markov-processes-gillespie
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run:

```bash
python run_simulation.py
```

On Windows:

```bash
py -m pip install -r requirements.txt
py run_simulation.py
```

The generated figure is saved to:

```text
results/gillespie_vs_analytical.png
```

---

## Technologies

- Python
- NumPy
- Matplotlib
- Continuous-time Markov chains
- Gillespie simulation
- Eigenvalue analysis
- Stochastic processes

---

## Quantitative Relevance

Continuous-time Markov processes appear in areas including:

- regime-switching models
- credit-rating transitions
- default modelling
- queueing systems
- reliability modelling
- event-driven stochastic systems

The project demonstrates both analytical and simulation-based approaches to stochastic modelling.

---

## Key Takeaways

- The generator matrix conserves total probability.
- The zero eigenvalue determines the stationary distribution.
- The negative eigenvalue determines the rate of convergence.
- The system satisfies detailed balance.
- KL divergence provides a convergence argument.
- Gillespie simulation closely reproduces the analytical solution.
- Averaging random trajectories recovers the deterministic probability dynamics.

---

## Disclaimer

This repository is an educational and portfolio project focused on stochastic processes and simulation.

All trajectories are synthetically generated from the specified transition rates.