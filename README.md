# Markov Processes and Gillespie Algorithm

A two-state continuous-time Markov process analysed both **analytically** and through **stochastic simulation with the Gillespie algorithm**.

The project develops the master equation, generator matrix, stationary distribution, eigenstructure, detailed-balance condition, Lyapunov argument using Kullback-Leibler divergence, symmetric reformulation, time-dependent solution and relaxation time, before comparing the analytical probabilities with Gillespie simulation.

---

## Overview

The system has two states with transition rates

$$
W(1 \to 2) = \alpha,
\qquad
W(2 \to 1) = \beta,
$$

where

$$
\alpha > 0,
\qquad
\beta > 0.
$$

The state probabilities are

$$
p(t)
=
\begin{pmatrix}
p_1(t) \\
p_2(t)
\end{pmatrix}.
$$

The project studies how these probabilities evolve over time, what stationary distribution they approach, and whether a stochastic Gillespie simulation reproduces the analytical solution.

---

## Master Equation

For state 1,

$$
\frac{dp_1}{dt}
=
\beta p_2 - \alpha p_1.
$$

For state 2,

$$
\frac{dp_2}{dt}
=
\alpha p_1 - \beta p_2.
$$

These can be written in matrix form as

$$
\frac{dp}{dt}
=
Wp(t),
$$

with generator matrix

$$
W
=
\begin{pmatrix}
-\alpha & \beta \\
\alpha & -\beta
\end{pmatrix}.
$$

Important properties of this generator are:

- each column sums to zero, so total probability is conserved
- off-diagonal entries are non-negative
- one eigenvalue is zero
- the second eigenvalue is strictly negative

---

## Eigenvalues

The characteristic equation is

$$
\det(W-\lambda I)
=
\lambda(\lambda+\alpha+\beta).
$$

Therefore,

$$
\lambda_1 = 0,
\qquad
\lambda_2 = -(\alpha+\beta).
$$

Because $\alpha,\beta>0$, the second eigenvalue is negative.

The zero eigenvalue corresponds to the stationary distribution, while the negative eigenvalue governs relaxation toward equilibrium.

---

## Stationary Distribution

The stationary state $\pi$ satisfies

$$
W\pi = 0
$$

together with

$$
\pi_1+\pi_2=1.
$$

Solving gives

$$
\pi_1
=
\frac{\beta}{\alpha+\beta},
$$

and

$$
\pi_2
=
\frac{\alpha}{\alpha+\beta}.
$$

For the parameters used in the simulation,

$$
\alpha=0.3,
\qquad
\beta=0.5,
$$

so

$$
\pi_1 = 0.625,
\qquad
\pi_2 = 0.375.
$$

---

## Left and Right Eigenvectors

For $\lambda_1=0$, the right eigenvector is the stationary distribution

$$
\pi
=
\begin{pmatrix}
\frac{\beta}{\alpha+\beta} \\
\frac{\alpha}{\alpha+\beta}
\end{pmatrix}.
$$

The corresponding left eigenvector is

$$
u^T=(1,1),
$$

which reflects conservation of total probability.

For

$$
\lambda_2=-(\alpha+\beta),
$$

a corresponding right eigenvector is

$$
\begin{pmatrix}
-1 \\
1
\end{pmatrix}.
$$

---

## Detailed Balance

Detailed balance requires

$$
\pi_1 W(1\to2)
=
\pi_2 W(2\to1).
$$

Substituting the rates and stationary probabilities gives

$$
\frac{\beta}{\alpha+\beta}\alpha
=
\frac{\alpha}{\alpha+\beta}\beta.
$$

Therefore the two-state process satisfies detailed balance.

---

## KL Divergence and Convergence

The Kullback-Leibler divergence between the time-dependent distribution and the stationary distribution is

$$
D(p(t)\|\pi)
=
\sum_{n=1}^{2}
p_n(t)
\log
\left(
\frac{p_n(t)}{\pi_n}
\right).
$$

Using the master equations and detailed balance, the project derives

$$
\frac{d}{dt}D(p(t)\|\pi)
=
-
(p_1\alpha-p_2\beta)
\log
\left(
\frac{p_1\alpha}{p_2\beta}
\right)
\le 0.
$$

Therefore, $D(p(t)\|\pi)$ is non-increasing and acts as a Lyapunov function.

Equality occurs only when

$$
p_1\alpha=p_2\beta,
$$

which is the stationary condition.

Hence,

$$
\lim_{t\to\infty}p(t)=\pi.
$$

---

## Symmetric Reformulation

Using the stationary probabilities, the generator can be transformed into the symmetric matrix

$$
M
=
\begin{pmatrix}
-\alpha & \sqrt{\alpha\beta} \\
\sqrt{\alpha\beta} & -\beta
\end{pmatrix}.
$$

This can be written as a similarity transform

$$
M=D^{-1}WD,
$$

where $D$ is diagonal with entries based on the square roots of the stationary probabilities.

Since $M$ is real and symmetric, it is diagonalisable by the spectral theorem.

---

## Time-Dependent Solution

The general solution is formed from the two eigenmodes.

The resulting probabilities are

$$
p_1(t)
=
\pi_1
+
\left(
p_1(0)-\pi_1
\right)
e^{-(\alpha+\beta)t},
$$

and

$$
p_2(t)
=
\pi_2
+
\left(
p_2(0)-\pi_2
\right)
e^{-(\alpha+\beta)t}.
$$

The probabilities therefore approach the stationary distribution exponentially.

---

## Relaxation Time

The decay rate is

$$
\alpha+\beta.
$$

The relaxation timescale is therefore

$$
\tau_{\text{relax}}
=
\frac{1}{\alpha+\beta}.
$$

Larger transition rates lead to faster convergence toward the stationary distribution.

---

## Gillespie Algorithm

The Gillespie algorithm simulates individual trajectories of the continuous-time Markov process.

When the system is in its current state, the waiting time until the next jump is exponentially distributed.

For this two-state process:

- in state 1, the outgoing rate is $\alpha$
- in state 2, the outgoing rate is $\beta$
- after each waiting time, the system must jump to the other state

The submitted experiment uses:

```text
alpha = 0.3
beta = 0.5
t_end = 15
trajectories = 10000
time bins = 300
starting state = 1
```

The simulated probabilities $p_1(t)$ and $p_2(t)$ are estimated by averaging across the independent Gillespie trajectories.

---

## Submitted Gillespie Implementation

The core simulation logic is retained in the repository:

```python
def gillespie(alpha, beta, t_end, starting_state=1):
    t = 0
    current_state = starting_state
    times = [0]
    states = [current_state]

    rates = {1: alpha, 2: beta}

    while t <= t_end:
        rate = rates[current_state]
        tau = np.random.exponential(1.0 / rate)
        t += tau

        current_state = 2 if current_state == 1 else 1

        times.append(t)
        states.append(current_state)

    return times, states
```

Probability estimates are then formed by evaluating the state of each trajectory across a common time grid and averaging across all runs.

---

## Analytical vs Simulated Results

For

$$
\alpha=0.3,
\qquad
\beta=0.5,
$$

the analytical stationary values are

$$
\pi_1=0.6250,
\qquad
\pi_2=0.3750.
$$

In the submitted experiment, the simulated probabilities converged to approximately

$$
p_1=0.6247,
\qquad
p_2=0.3753.
$$

The simulated curves closely overlap the analytical time-dependent solution, showing that the Gillespie algorithm accurately reproduces the theoretical continuous-time Markov dynamics.

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

### `src/gillespie.py`

Contains the core functions from the submitted computational work:

- Gillespie trajectory simulation
- probability estimation across many trajectories
- analytical time-dependent solution

### `run_simulation.py`

Provides a reproducible portfolio runner using the same main parameter values as the submitted experiment and saves a comparison plot.

### `results/`

Stores generated output figures.

---

## Running the Project

Clone the repository:

```bash
git clone https://github.com/aamina-ahmad/markov-processes-gillespie.git
```

Move into the project directory:

```bash
cd markov-processes-gillespie
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the simulation:

```bash
python run_simulation.py
```

On Windows:

```bash
py -m pip install -r requirements.txt
py run_simulation.py
```

The comparison plot will be written to:

```text
results/gillespie_vs_analytical.png
```

---

## Technologies

- Python
- NumPy
- Matplotlib
- continuous-time Markov chains
- stochastic simulation
- Gillespie algorithm
- eigenvalue analysis
- analytical ODE solutions

---

## Quantitative Relevance

Continuous-time Markov models appear across quantitative modelling wherever systems move randomly between discrete states.

Examples include:

- credit-rating transitions
- regime-switching models
- queueing systems
- reliability models
- order-state models
- population processes
- chemical reaction networks
- event-driven stochastic systems

The project demonstrates both sides of the modelling problem:

1. deriving the probability dynamics analytically
2. reproducing those dynamics through pathwise stochastic simulation

---

## Key Takeaways

1. The two-state generator conserves total probability and has one zero eigenvalue.
2. The zero eigenvalue determines the unique stationary distribution.
3. Detailed balance holds for the two-state process.
4. KL divergence decreases toward equilibrium and acts as a Lyapunov function.
5. The non-zero eigenvalue determines exponential relaxation.
6. The relaxation timescale is $1/(\alpha+\beta)$.
7. Gillespie simulation reproduces the analytical time-dependent probabilities closely.
8. Averaging many stochastic trajectories recovers the deterministic probability evolution predicted by the master equation.

---

## Disclaimer

This repository is an educational and portfolio project focused on stochastic processes and simulation.

The simulation uses synthetic trajectories generated from the specified transition rates.
