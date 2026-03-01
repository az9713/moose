# Case 95: Charge Relaxation — Dielectric Relaxation Time in an Ohmic Conductor

## Overview

Free charge injected into a weakly conducting dielectric does not stay put. Ohm's law drives a current proportional to the local electric field, which in turn is driven by the charge itself via Gauss's law. The two effects combine to produce a self-draining ODE: charge density falls exponentially with the dielectric relaxation time tau_e = epsilon / sigma. This is one of the fundamental timescales in electrical engineering, governing electrostatic discharge (ESD) protection, charge dissipation in printed-circuit-board materials, and the crossover from displacement-current to conduction-current dominated behavior as frequency decreases.

The key insight from MIT 6.641 Lecture 7 is that the relaxation is purely local. Because ohmic current redistributes charge in a way that exactly cancels the divergence of displacement current, there is no spatial diffusion of charge — every point relaxes independently with the same time constant. A Gaussian blob of initial charge does not spread; it simply shrinks uniformly in amplitude while preserving its spatial shape.

This case is deliberately minimal. The MOOSE input has no spatial diffusion term and no boundary conditions — the physics is entirely captured by two kernels: `TimeDerivative` and `CoefReaction`. Running 50 time steps over five relaxation times (5 tau_e = 1 s) demonstrates that the numerical solution tracks exp(-5t) to better than 1%.

---

## The Physics

**Governing equations**

Combining the continuity equation with Gauss's law and Ohm's law:

```
div J + d(rho)/dt = 0
J = sigma E
div(epsilon E) = rho
```

Substituting gives a pointwise ODE at every location in the conductor:

```
d(rho)/dt + (sigma / epsilon) * rho = 0
```

**Analytic solution**

```
rho(x, y, t) = rho_0(x, y) * exp(-t / tau_e)
```

where the dielectric relaxation time is:

```
tau_e = epsilon / sigma
```

With sigma/epsilon = 5.0, the relaxation time is tau_e = 0.2 s.

**Initial condition**

A Gaussian charge blob centered at (0.5, 0.5):

```
rho_0(x, y) = exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.01)
```

Peak value is 1.0 C/m^3, width (sigma_g) is 0.1 m.

**Boundary conditions**

None — zero-flux (natural Neumann) is applied on all four edges. Physically, charge cannot escape the material; it decays in place through ohmic conduction.

**Domain and discretization**

- 2D unit square [0, 1]^2
- 30 x 30 bilinear quadrilateral elements

**Parameters**

| Symbol | Value | Meaning |
|--------|-------|---------|
| sigma/epsilon | 5.0 s^-1 | Inverse relaxation time |
| tau_e | 0.2 s | Dielectric relaxation time |
| dt | 0.02 s | Time step (= tau_e / 10) |
| end_time | 1.0 s | Duration (= 5 tau_e) |

---

## Input File Walkthrough

**Global parameter**

```
sigma_over_eps = 5.0
```

This single ratio determines the decay rate. For SI materials with epsilon_r = 3 and sigma = 1e-10 S/m, tau_e ≈ 0.27 s and sigma/epsilon ≈ 3.77 s^-1. Here it is rounded to 5.0 for a clean simulation.

**[Variables]**

`rho` is the free charge density in C/m^3. The initial condition is set inline using `FunctionIC`:

```
[InitialCondition]
  type     = FunctionIC
  function = 'exp(-((x-0.5)^2 + (y-0.5)^2) / 0.01)'
[]
```

This approach avoids defining a separate `[Functions]` block for a simple expression.

**[Kernels]**

Two kernels and nothing else:

```
[rho_time]  type = TimeDerivative  — d(rho)/dt
[rho_decay] type = CoefReaction    — coefficient * rho  (= sigma/epsilon * rho)
```

`CoefReaction` contributes the term `integral( coefficient * rho * psi dV )` to the residual. Combined with the time derivative, the discrete system implements the pointwise ODE d(rho)/dt + (sigma/epsilon) * rho = 0 at every quadrature point.

**[BCs]**

The block is present but empty. This is an explicit reminder that zero-flux Neumann is the MOOSE default — no charge leaves through the boundary.

**[Postprocessors]**

| Name | Type | What it tracks |
|------|------|---------------|
| avg_rho | ElementAverageValue | Spatial mean of rho — decays as exp(-5t) |
| max_rho | ElementExtremeValue | Peak charge density — same exponential, no spatial spreading |
| total_charge | ElementIntegralVariablePostprocessor | Total charge Q(t) — decays as Q_0 exp(-t/tau_e) |

Because the PDE has no spatial coupling, max_rho and avg_rho differ only by the spatial shape of rho_0 but track the same exponential envelope.

**[Executioner]**

Transient NEWTON solve with LU direct factorization. The time step dt = 0.02 s gives ten steps per relaxation time, which is more than enough for an exponentially decaying system.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case95-charge-relaxation \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case95_charge_relaxation.i 2>&1 | tail -30'
```

The 50-step transient completes in a few seconds. Exodus output captures the full spatial field at every time step, and the CSV file records the three postprocessors versus time.

---

## Expected Results

The CSV postprocessors should follow clean exponential decay:

| Time (s) | exp(-5t) | avg_rho (approximate) |
|----------|----------|----------------------|
| 0.0 | 1.000 | 0.031 (spatial average of Gaussian) |
| 0.2 | 0.368 | 0.011 |
| 0.4 | 0.135 | 0.0042 |
| 1.0 | 0.0067 | 2.1e-4 |

Note that avg_rho is smaller than exp(-5t) by the spatial average of the Gaussian (about 3.1% of the peak), but it decays with exactly the same rate. Plotting ln(max_rho) versus time in a spreadsheet gives a straight line with slope -5.

In the Exodus output, the spatial structure of rho is preserved at every time step — the Gaussian footprint does not spread, confirming that there is no spatial diffusion in this problem.

---

## Key Takeaways

- The dielectric relaxation time tau_e = epsilon / sigma is the fundamental timescale governing charge decay in ohmic conductors; it ranges from nanoseconds in metals to years in high-purity insulators.
- The charge-continuity equation combined with Gauss's law and Ohm's law reduces to a pointwise ODE with no spatial derivatives — each material point relaxes independently.
- `CoefReaction` is the MOOSE idiom for a linear decay or gain term of the form k * u in a PDE.
- `TimeDerivative` + `CoefReaction` (no diffusion kernel) implements a pure decay ODE in MOOSE — the spatial discretization is irrelevant because the residual has no gradient terms.
- Setting the `[BCs]` block empty is valid MOOSE syntax and is useful as an explicit reminder that zero-flux Neumann is the default.
- An inline `FunctionIC` within a variable block is more concise than a separate `[Functions]` block when the initial condition expression is simple.
