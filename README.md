# Markov Processes and Gillespie Algorithm

A two-state continuous-time Markov process analysed both **analytically** and through **stochastic simulation using the Gillespie algorithm**.

The project develops the master equation, generator matrix, stationary distribution, eigenstructure, detailed balance, a Kullback-Leibler divergence convergence argument, symmetric reformulation, the time-dependent solution and relaxation time, before comparing the analytical probabilities with Gillespie simulation.

---

## Overview

The system has two states with transition rates

$$
W(1 \to 2) = \alpha,
\qquad
W(2 \to 1) = \beta
$$

where

$$
\alpha > 0,
\qquad
\beta > 0
$$

The state-probability vector is

$$
p(t)
=
\begin{pmatrix}
p_1(t) \\
p_2(t)
\end{pmatrix}
$$

The project studies how these probabilities evolve over time, the stationary distribution they approach, and whether stochastic Gillespie simulation reproduces the analytical solution.

---

## Master Equation

For a continuous-time Markov process, the probability of occupying each state changes according to the transition rates.

For state 1,

$$
\frac{dp_1}{dt}
=
\beta p_2 - \alpha p_1
$$

For state 2,

$$
\frac{dp_2}{dt}
=
\alpha p_1 - \beta p_2
$$

These equations can be written in matrix form as

$$
\frac{dp}{dt}
=
Wp(t)
$$

with generator matrix

$$
W
=
\begin{pmatrix}
-\alpha & \beta \\
\alpha & -\beta
\end{pmatrix}
$$

Important properties of this generator are:

- each column sums to zero, so total probability is conserved
- off-diagonal entries are non-negative
- one eigenvalue is zero
- the second eigenvalue is strictly negative

---

## Probability Conservation

Because each column of $W$ sums to zero,

$$
\frac{d}{dt}
\left(
p_1(t)+p_2(t)
\right)
=
0
$$

Therefore,

$$
p_1(t)+p_2(t)=1
$$

for all $t$, provided the initial probabilities sum to one.

This reflects the fact that the system must always be in either state 1 or state 2.

---

## Eigenvalues

The characteristic equation is

$$
\det(W-\lambda I)
=
0
$$

For this generator,

$$
\det(W-\lambda I)
=
\lambda(\lambda+\alpha+\beta)
$$

Therefore,

$$
\lambda_1 = 0,
\qquad
\lambda_2 = -(\alpha+\beta)
$$

Because $\alpha,\beta>0$,

$$
\lambda_2<0
$$

The zero eigenvalue corresponds to the stationary distribution, while the negative eigenvalue governs relaxation toward equilibrium.

---

## Stationary Distribution

The stationary state $\pi$ satisfies

$$
W\pi = 0
$$

with

$$
\pi
=
\begin{pmatrix}
\pi_1 \\
\pi_2
\end{pmatrix}
$$

and the normalisation condition

$$
\pi_1+\pi_2=1
$$

From

$$
-\alpha\pi_1+\beta\pi_2=0
$$

we obtain

$$
\pi_2
=
\frac{\alpha}{\beta}\pi_1
$$

Using the normalisation condition gives

$$
\pi_1
=
\frac{\beta}{\alpha+\beta}
$$

and

$$
\pi_2
=
\frac{\alpha}{\alpha+\beta}
$$

Therefore,

$$
\pi
=
\begin{pmatrix}
\dfrac{\beta}{\alpha+\beta} \\
\dfrac{\alpha}{\alpha+\beta}
\end{pmatrix}
$$

---

## Parameter Values Used in the Simulation

The computational experiment uses

$$
\alpha=0.3,
\qquad
\beta=0.5
$$

Therefore,

$$
\pi_1
=
\frac{0.5}{0.3+0.5}
=
0.625
$$

and

$$
\pi_2
=
\frac{0.3}{0.3+0.5}
=
0.375
$$

So the analytical stationary distribution is

$$
\pi
=
\begin{pmatrix}
0.625 \\
0.375
\end{pmatrix}
$$

---

## Left and Right Eigenvectors

For a right eigenvector $v$,

$$
Wv=\lambda v
$$

For a left eigenvector $u^T$,

$$
u^T W=\lambda u^T
$$

### Right Eigenvector for the Zero Eigenvalue

For

$$
\lambda_1=0
$$

the right eigenvector is the stationary distribution

$$
\pi
=
\begin{pmatrix}
\dfrac{\beta}{\alpha+\beta} \\
\dfrac{\alpha}{\alpha+\beta}
\end{pmatrix}
$$

because

$$
W\pi=0
$$

### Left Eigenvector for the Zero Eigenvalue

The corresponding left eigenvector is

$$
u^T=(1,1)
$$

since

$$
(1,1)W=(0,0)
$$

This is directly related to conservation of total probability.

### Right Eigenvector for the Negative Eigenvalue

For

$$
\lambda_2=-(\alpha+\beta)
$$

a corresponding right eigenvector is

$$
v_2
=
\begin{pmatrix}
-1 \\
1
\end{pmatrix}
$$

---

## Detailed Balance

Detailed balance requires

$$
\pi_1 W(1\to2)
=
\pi_2 W(2\to1)
$$

Substituting the transition rates gives

$$
\pi_1\alpha
=
\pi_2\beta
$$

Using the stationary probabilities,

$$
\frac{\beta}{\alpha+\beta}\alpha
=
\frac{\alpha}{\alpha+\beta}\beta
$$

Both sides equal

$$
\frac{\alpha\beta}{\alpha+\beta}
$$

so detailed balance is satisfied.

For this two-state process, there is only one probability flux between the two states, and at equilibrium the two directional fluxes balance exactly.

---

## Kullback-Leibler Divergence

The Kullback-Leibler divergence between the time-dependent distribution $p(t)$ and the stationary distribution $\pi$ is

$$
D(p(t)\|\pi)
=
\sum_{n=1}^{2}
p_n(t)
\log
\left(
\frac{p_n(t)}{\pi_n}
\right)
$$

For the two-state system,

$$
D(p(t)\|\pi)
=
p_1(t)
\log
\left(
\frac{p_1(t)}{\pi_1}
\right)
+
p_2(t)
\log
\left(
\frac{p_2(t)}{\pi_2}
\right)
$$

Using the master equations and detailed balance, the derivative can be written as

$$
\frac{d}{dt}D(p(t)\|\pi)
=
-
(p_1\alpha-p_2\beta)
\log
\left(
\frac{p_1\alpha}{p_2\beta}
\right)
$$

Using the identity

$$
(e^x-e^y)(x-y)\ge0
$$

gives

$$
(p_1\alpha-p_2\beta)
\log
\left(
\frac{p_1\alpha}{p_2\beta}
\right)
\ge0
$$

and therefore

$$
\frac{d}{dt}D(p(t)\|\pi)
\le0
$$

So the KL divergence is non-increasing with time.

Equality occurs only when

$$
p_1\alpha=p_2\beta
$$

which is precisely the stationary condition.

Therefore,

$$
\lim_{t\to\infty}p(t)=\pi
$$

and the KL divergence acts as a Lyapunov function for the dynamics.

---

## Symmetric Reformulation

Using the stationary probabilities, the generator can be transformed into the symmetric matrix

$$
M
=
\begin{pmatrix}
-\alpha & \sqrt{\alpha\beta} \\
\sqrt{\alpha\beta} & -\beta
\end{pmatrix}
$$

This matrix is related to $W$ by the similarity transform

$$
M=D^{-1}WD
$$

where $D$ is diagonal with entries related to the square roots of the stationary probabilities.

Because $M$ is real and symmetric,

$$
M_{12}=M_{21}
$$

and the spectral theorem guarantees that it is diagonalisable.

If

$$
q(t)=D^{-1}p(t)
$$

then

$$
\frac{dq}{dt}
=
Mq(t)
$$

which gives an equivalent symmetric representation of the dynamics.

---

## Diagonalisation

The two eigenvalues are

$$
\lambda_1=0
$$

and

$$
\lambda_2=-(\alpha+\beta)
$$

with corresponding eigenvectors

$$
v_1
=
\begin{pmatrix}
\pi_1 \\
\pi_2
\end{pmatrix}
$$

and

$$
v_2
=
\begin{pmatrix}
-1 \\
1
\end{pmatrix}
$$

The eigenvector matrix can therefore be written as

$$
P
=
\begin{pmatrix}
\pi_1 & -1 \\
\pi_2 & 1
\end{pmatrix}
$$

with diagonal eigenvalue matrix

$$
\Lambda
=
\begin{pmatrix}
0 & 0 \\
0 & -(\alpha+\beta)
\end{pmatrix}
$$

---

## Time-Dependent Solution

Because the generator is diagonalisable, the general solution is a linear combination of the two eigenmodes.

The resulting probabilities are

$$
p_1(t)
=
\pi_1
+
\left(
p_1(0)-\pi_1
\right)
e^{-(\alpha+\beta)t}
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
e^{-(\alpha+\beta)t}
$$

The transient component therefore decays exponentially at rate

$$
\alpha+\beta
$$

As

$$
t\to\infty
$$

the exponential term vanishes and

$$
p_1(t)\to\pi_1
$$

and

$$
p_2(t)\to\pi_2
$$

---

## Relaxation Time

The non-zero eigenvalue is

$$
\lambda_2=-(\alpha+\beta)
$$

so the characteristic relaxation time is

$$
\tau_{\text{relax}}
=
\frac{1}{|\lambda_2|}
$$

Therefore,

$$
\tau_{\text{relax}}
=
\frac{1}{\alpha+\beta}
$$

For the simulation parameters,

$$
\alpha+\beta=0.8
$$

so

$$
\tau_{\text{relax}}
=
\frac{1}{0.8}
=
1.25
$$

Larger transition rates therefore lead to faster convergence toward the stationary distribution.

---

## Gillespie Algorithm

The Gillespie algorithm provides a stochastic simulation of the continuous-time Markov process.

At any point in time, the system occupies one state.

The waiting time until the next jump is exponentially distributed.

If the system is in state 1, the rate is

$$
\alpha
$$

and if it is in state 2, the rate is

$$
\beta
$$

The waiting time is generated using

```python
tau = np.random.exponential(1.0 / rate)
```

Because the process has only two states, the next state is always the opposite state.

---

## Gillespie Simulation Parameters

The simulation uses:

```text
alpha = 0.3
beta = 0.5
t_end = 15
trajectories = 10000
time bins = 300
starting state = 1
```

The probabilities $p_1(t)$ and $p_2(t)$ are estimated by averaging over 10,000 independent Gillespie trajectories.

---

## Gillespie Implementation

The core simulation function is:

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

The probability-estimation routine then evaluates each trajectory on a common time grid and counts the proportion of trajectories occupying each state.

---

## Probability Estimation

For each independent trajectory:

1. a Gillespie path is simulated
2. the state is evaluated across a common time grid
3. the state-1 and state-2 counts are recorded
4. counts are divided by the number of trajectories

The estimated probabilities are therefore

$$
\hat{p}_1(t)
=
\frac{\text{number of trajectories in state 1 at time }t}
{\text{total number of trajectories}}
$$

and

$$
\hat{p}_2(t)
=
\frac{\text{number of trajectories in state 2 at time }t}
{\text{total number of trajectories}}
$$

---

## Analytical vs Simulated Results

For

$$
\alpha=0.3,
\qquad
\beta=0.5
$$

the analytical stationary probabilities are

$$
\pi_1=0.6250
$$

and

$$
\pi_2=0.3750
$$

In the submitted experiment, the simulated probabilities converged approximately to

$$
p_1=0.6247
$$

and

$$
p_2=0.3753
$$

The simulated trajectories therefore converge very closely to the theoretical stationary distribution.

The analytical time-dependent curves also closely overlap the Gillespie probability estimates.

This confirms that the stochastic simulation reproduces the dynamics predicted by the master equation.

---

## Project Workflow

```mermaid
flowchart LR
    A[Define alpha and beta] --> B[Construct generator matrix W]

    B --> C1[Analytical solution]
    B --> C2[Gillespie simulation]

    C1 --> D1[Stationary distribution]
    C1 --> D2[Eigenvalues and eigenvectors]
    C1 --> D3[Time-dependent probabilities]
    C1 --> D4[Relaxation time]

    C2 --> E[Simulate many stochastic trajectories]
    E --> F[Estimate p1(t) and p2(t)]

    D3 --> G[Compare analytical and simulated probabilities]
    F --> G

    G --> H[Verify convergence to stationary distribution]
```

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

Contains the computational functions for:

- Gillespie trajectory simulation
- probability estimation across independent trajectories
- the analytical time-dependent solution

### `run_simulation.py`

Runs the main experiment using the project parameters and generates the analytical-versus-simulation comparison plot.

### `results/`

Stores generated figures.

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

Install the required packages:

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

The generated plot will be saved as:

```text
results/gillespie_vs_analytical.png
```

---

## Expected Output

The program prints the analytical stationary probabilities and the final simulated probability estimates.

The analytical values are:

```text
pi1 = 0.6250
pi2 = 0.3750
```

The generated figure compares:

- simulated $p_1(t)$
- simulated $p_2(t)$
- analytical $p_1(t)$
- analytical $p_2(t)$
- stationary probability levels

---

## Technologies

The project uses:

- Python
- NumPy
- Matplotlib
- continuous-time Markov chains
- stochastic simulation
- Gillespie algorithm
- linear ODE systems
- eigenvalue analysis

---

## Quantitative Relevance

Continuous-time Markov models are useful whenever a system moves randomly between discrete states.

Applications include:

- credit-rating transitions
- regime-switching models
- default-state modelling
- queueing systems
- reliability models
- event-driven systems
- population processes
- chemical reaction networks

The project demonstrates both analytical and computational approaches to stochastic modelling.

The analytical model describes the evolution of state probabilities, while the Gillespie algorithm generates individual random paths whose average behaviour recovers those probabilities.

---

## Analytical and Simulation Perspectives

The project illustrates two complementary views of the same stochastic process.

### Distribution-Level View

The master equation describes how the probabilities evolve:

$$
\frac{dp}{dt}=Wp(t)
$$

This gives a deterministic evolution for the probability distribution.

### Pathwise View

The Gillespie algorithm generates individual random sample paths.

Each path jumps unpredictably between the two states, but averaging over many independent trajectories produces the distribution predicted by the master equation.

This provides a direct numerical demonstration of the relationship between microscopic stochastic trajectories and macroscopic probability evolution.

---

## Key Takeaways

1. The generator matrix conserves total probability.
2. The system has one zero eigenvalue and one strictly negative eigenvalue.
3. The zero eigenvalue determines the stationary distribution.
4. The stationary distribution satisfies detailed balance.
5. KL divergence decreases toward equilibrium and acts as a Lyapunov function.
6. The transient probability mode decays exponentially.
7. The relaxation time is

$$
\frac{1}{\alpha+\beta}
$$

8. Gillespie simulation reproduces the analytical time-dependent probabilities.
9. Averaging many stochastic trajectories recovers the deterministic probability evolution predicted by the master equation.

---

## Disclaimer

This repository is an educational and portfolio project focused on stochastic processes and simulation.

All trajectories are synthetically generated from the specified transition rates.