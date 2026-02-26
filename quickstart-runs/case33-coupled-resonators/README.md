# Case 33: Coupled Resonator Beating — Energy Exchange

## Overview

Two resonators (optical cavities, LC circuits, mechanical oscillators) coupled
by a small interaction κ exchange energy back and forth at the **beating
frequency** 2κ. This is the temporal analogue of spatial interference. If the
resonators also have loss rate γ, the energy envelope decays exponentially
while oscillating — producing a damped beat pattern.

This case models two coupled field amplitudes u and v using two scalar MOOSE
variables with diffusion, linear decay, and mutual coupling kernels. The
coupled-mode equations derive from Haus, *Electromagnetic Noise and Quantum
Optical Measurements* (Springer, 2000), Ch. 3. Because the equations are linear,
exact analytical solutions exist and provide direct verification.

---

## The Physics

### Governing Equations

The coupled-mode equations for two resonator amplitudes on a 2D domain:

```
∂u/∂t = D·∇²u − γ·u + κ·v
∂v/∂t = D·∇²v − γ·v + κ·u
```

Parameters:

| Symbol | Value | Meaning |
|--------|-------|---------|
| D      | 0.01  | Spatial diffusion coefficient |
| γ      | 0.5   | Energy decay rate (loss) |
| κ      | 3.0   | Inter-resonator coupling coefficient |

### Analytical Solution (Spatially Uniform Modes)

For uniform initial conditions u₀ = 1, v₀ = 0, the spatially homogeneous
solution is:

```
u(t) = exp(−γt) · cos(κt)
v(t) = exp(−γt) · sin(κt)
```

The energy oscillates between the two resonators with period T = π/κ ≈ 1.047 s
while the total energy (u² + v²) decays as exp(−2γt).

### Normal Modes

Adding and subtracting the two equations reveals the normal modes:

```
p = u + v  decays at rate γ − κ = −2.5  (grows if κ > γ)
q = u − v  decays at rate γ + κ =  3.5  (always decays)
```

The fast mode q dies quickly; the slow mode p persists (and may even grow if
κ > γ). The beating observed in u and v is the superposition of these two
normal modes at different rates.

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 20×20 QUAD4 | Unit square domain |
| `ADTimeDerivative` on u, v | ∂u/∂t and ∂v/∂t |
| `ADMatDiffusion` on u, v | D·∇²u and D·∇²v |
| `ADReaction` (rate=γ) on u, v | Linear decay −γ·u and −γ·v |
| `CoupledForce` (coeff=κ, v→u) | Coupling term +κ·v in u equation |
| `CoupledForce` (coeff=κ, u→v) | Coupling term +κ·u in v equation |
| `ConstantIC` u=1, v=0 | Initial conditions |
| `ADDirichletBC` u=v=0 on walls | Homogeneous Dirichlet boundaries |
| `ElementAverageValue` avg_u, avg_v | Track beat pattern over time |
| `Transient` executioner, dt=0.05 | Time integration to t=6s |

`CoupledForce` is the standard MOOSE kernel for injecting one variable as a
forcing term in another variable's residual. Its `coeff` parameter sets κ.
`ADReaction` with `rate = ${gamma}` implements the linear decay −γ·u.

---

## Expected Results

The CSV postprocessor output should show:

| Time | avg_u (expected) | avg_v (expected) | Physical state |
|------|-----------------|-----------------|----------------|
| 0.0  | ~1.0            | ~0.0            | All energy in u |
| 0.52 | ~0.0            | ~0.75           | Energy transferred to v |
| 1.05 | ~−0.57          | ~0.0            | Energy back in u (negative phase) |
| 2.0  | ~0.10           | ~0.27           | Beat amplitude decayed |
| 6.0  | ~0.0            | ~0.0            | Fully dissipated |

The beating period π/κ ≈ 1.047 s is visible in the oscillation of avg_u and avg_v.
Total energy avg_u² + avg_v² should decrease monotonically as exp(−2γt).
The spatial diffusion D = 0.01 is small enough that it only slightly smooths the
initial condition; the beating is dominated by the κ coupling term.

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case33-coupled-resonators \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case33_coupled_resonators.i 2>&1 | tail -20'
```

Output files produced:
- `case33_coupled_resonators_out.e` — Exodus file with u(x,y,t) and v(x,y,t)
- `case33_coupled_resonators_out.csv` — avg_u, avg_v vs. time (beating pattern)
