# Case 22: Charge Relaxation in an Ohmic Medium

## Overview

When free electric charge is placed inside an electrically conducting medium it does not
stay there — it migrates to the surfaces under the influence of its own electric field and
is neutralised by conduction current. This process, called **charge relaxation**, is one
of the fundamental time-scale problems of continuum electromechanics. The governing
physics is derived in Chapter 5 of Melcher's *Continuum Electromechanics* (MIT Press,
1981), specifically §5.9–5.10. The key result is that any initial free charge distribution
decays pointwise as a pure exponential in time, with the **relaxation time** τ_e = ε/σ,
where ε is the permittivity and σ is the electrical conductivity. For a good conductor
like copper (σ ≈ 6×10⁷ S/m, ε ≈ ε₀ ≈ 8.85×10⁻¹² F/m) the relaxation time is of order
10⁻¹⁹ s — instantaneous on any engineering time scale. For a poor conductor or a dielectric
liquid the time scale can be milliseconds to seconds, which is what this case models.

This case uses two coupled MOOSE variables to capture the full field problem. The free
charge density ρ_e is the primary evolution variable, governed by a simple first-order
ODE at each spatial point that is expressed through two kernels: `ADTimeDerivative` for
the rate term and `ADReaction` for the decay term. The electric potential φ is determined
at each instant by Poisson's equation, where the charge density acts as a source. Because
φ responds to ρ_e but does not feed back into ρ_e (in the linear, spatially uniform
conductivity limit), the physics decouple elegantly — ρ_e decays independently of φ, and
φ simply tracks the decaying charge.

Several MOOSE concepts are introduced for the first time in this case: using `ADReaction`
to model a linear volumetric sink, using `ADCoupledForce` to inject one variable as a
source term in another variable's equation, using `ADHeatConduction` as a generic
diffusion operator for a non-thermal problem (here Poisson's equation for the potential),
and using a HIT top-level variable to parameterise the decay rate so that it can be
changed in a single place.

---

## The Physics

### Governing Equations

The charge relaxation problem in a spatially uniform, linear, Ohmic medium consists of
two coupled equations:

```
∂ρ_e/∂t + (σ/ε)·ρ_e = 0          (charge relaxation)

-div(ε·grad(φ)) = ρ_e              (Poisson's equation)

E = -grad(φ)                        (electric field from potential)
```

The first equation is a first-order linear ODE in time at each spatial point. It has the
exact solution:

```
ρ_e(x, y, t) = ρ_e(x, y, 0) · exp(−t/τ_e)    where τ_e = ε/σ
```

The spatial distribution of the charge density is "frozen" — every point decays at the
same rate, preserving the shape of the initial distribution exactly. For a Gaussian blob
the blob simply shrinks in amplitude without changing its width.

The Poisson equation determines the electric potential from the instantaneous charge
distribution. Because ρ_e decays exponentially, so does the potential:

```
φ(x, y, t) = φ(x, y, 0) · exp(−t/τ_e)
```

where φ(x, y, 0) is the potential corresponding to the initial charge distribution via
Poisson's equation. Both fields collapse to zero at the same exponential rate.

### What Drives the Decay

In an Ohmic medium the conduction current density is J = σE. The current divergence
∇·J = σ∇·E = σρ_e/ε (by Gauss's law: ∇·εE = ρ_e). The conservation of charge law
∂ρ_e/∂t + ∇·J = 0 then gives directly ∂ρ_e/∂t = −(σ/ε)ρ_e. Charge is carried away
by the electric field it creates; the faster the medium conducts, the faster the field
drives the charge to the boundaries and neutralises it.

### Boundary Conditions

| Variable | Boundary       | Type      | Value | Physical meaning          |
|----------|----------------|-----------|-------|---------------------------|
| φ        | all four walls | Dirichlet | 0 V   | Grounded conducting walls |
| ρ_e      | all four walls | (natural) | zero flux | Charge cannot leave through walls |

The natural (zero-flux Neumann) boundary condition on ρ_e is the physically correct
choice: in the continuum model the charge decays in place by Ohmic conduction, it does
not flow through the boundary. The grounded walls set the reference potential.

### Parameters

| Symbol  | Name               | Value | Units     |
|---------|--------------------|-------|-----------|
| ε       | Permittivity       | 1.0   | F/m       |
| σ       | Conductivity       | 10.0  | S/m       |
| σ/ε     | Decay rate         | 10.0  | s⁻¹       |
| τ_e     | Relaxation time    | 0.1   | s         |
| L       | Domain side length | 1.0   | m         |

### Domain and Initial Condition

The domain is the unit square [0,1]×[0,1] discretised with a 30×30 uniform mesh of
quadrilateral (QUAD4) elements.

The initial free charge distribution is a Gaussian blob centred at (0.5, 0.5):

```
ρ_e(x, y, 0) = exp(−((x−0.5)² + (y−0.5)²) / 0.01)
```

The width parameter 0.01 gives a 1/e radius of 0.1 m. The peak value is 1. The initial
potential φ = 0 everywhere; the Poisson solve at the first timestep computes the actual
initial potential consistent with the charge distribution.

---

## Input File Walkthrough

### HIT Top-Level Variable

```
sigma_over_eps = 10.0
```

Defining the σ/ε ratio at the top of the file as a HIT variable means it is referenced
with `${sigma_over_eps}` inside the `ADReaction` kernel. To change the material (and thus
the relaxation time) only this one number needs to be edited.

### `[Mesh]`

A 30×30 uniform quadrilateral mesh over the unit square. The mesh is fine enough to
resolve the Gaussian blob (width ~0.1 m, captured by ~3 elements across the 1/e radius)
without being unnecessarily expensive.

### `[Variables]`

Two nodal (FIRST-order LAGRANGE) variables:

```
rho_e  --  free charge density [C/m²]
             IC: Gaussian blob exp(-((x-0.5)^2+(y-0.5)^2)/0.01)
             Governed by the relaxation ODE

phi    --  electric potential [V]
             IC: 0 everywhere
             Governed by Poisson's equation
```

The `FunctionIC` for `rho_e` evaluates an inline mathematical expression at every node.
No custom function block is needed because HIT function expressions can be written
directly inside the `InitialCondition` sub-block.

### `[Kernels]`

Four kernels define the two coupled PDEs:

**Charge relaxation equation (variable = rho_e):**

```
[rho_e_time]   ADTimeDerivative   -- ∂ρ_e/∂t
[rho_e_decay]  ADReaction         -- +(σ/ε)·ρ_e  (linear volumetric sink)
```

`ADReaction` adds ∫ rate·ρ_e·φ_i dV to the residual. With rate = σ/ε = 10, this is the
exact weak form of (σ/ε)·ρ_e. The `AD` prefix means automatic differentiation computes
the exact Jacobian contribution, enabling quadratic Newton convergence.

**Poisson equation (variable = phi):**

```
[phi_laplacian] ADHeatConduction   -- -div(ε·grad(φ))  (diffusion with ε as conductivity)
[phi_source]    ADCoupledForce     -- -ρ_e  (charge density as source)
```

`ADHeatConduction` solves −div(k·grad(u)) where k is a named material property. By
setting `thermal_conductivity = permittivity`, the same kernel that computes heat
diffusion is repurposed for Poisson's equation. The name "thermal_conductivity" is just
the kernel's internal property lookup — the actual material property is `permittivity`.

`ADCoupledForce` adds −∫ v·φ_i dV to the residual of the named variable, where v is
another MOOSE variable. This injects ρ_e as a source on the right-hand side of Poisson's
equation. (MOOSE residuals are written in the form R = 0, so the sign convention means
`ADCoupledForce` effectively adds ρ_e to the right-hand side.)

### `[BCs]`

One BC block covers all four boundaries for φ:

```
[phi_ground]  ADDirichletBC  phi = 0  on left, right, top, bottom
```

No BC block is needed for `rho_e`. The natural boundary condition (zero normal flux) is
automatically enforced when no BC is specified, which is exactly the physical condition.

### `[Materials]`

```
[permittivity]  ADGenericConstantMaterial  prop_names='permittivity'  prop_values='1.0'
```

A single material provides the permittivity used by `ADHeatConduction`. The `AD` prefix
(automatic differentiation) is required because the kernels are also `AD` objects and
they request material properties through the AD material system.

### `[Preconditioning]`

```
[SMP]  type = SMP  full = true
```

Single-Matrix Preconditioning with `full = true` builds a single global preconditioner
that includes off-diagonal Jacobian blocks from the rho_e–phi coupling. This is
important for robust convergence even though the coupling is one-way (ρ_e → φ), because
the Poisson equation has no time derivative and the system is not block-diagonal.

### `[Executioner]`

```
type       = Transient
solve_type = NEWTON
dt         = 0.01   (one-tenth of the relaxation time)
end_time   = 0.5    (five relaxation times)
```

The fixed timestep dt = 0.01 s gives 10 points per relaxation time, which is adequate
to resolve the exponential decay. Five relaxation times bring the charge to exp(−5) ≈
0.007 of its initial value — effectively zero for practical purposes.

### `[Postprocessors]`

```
avg_rho_e   ElementAverageValue  -- domain-mean charge density
max_rho_e   ElementExtremeValue  -- peak charge density (tracks the Gaussian maximum)
max_phi     ElementExtremeValue  -- peak electric potential
```

All three quantities should show the same exponential decay with time constant τ = 0.1 s.
The CSV output allows easy verification against the analytical solution.

---

## Running the Simulation

This case requires the `combined` application (or any MOOSE app built with the
`HEAT_TRANSFER` module, which provides `ADHeatConduction`). The base framework kernels
`ADTimeDerivative`, `ADReaction`, and `ADCoupledForce` are always available.

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case22-charge-relaxation \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case22_charge_relaxation.i 2>&1 | tail -30'
```

Expected console output: 50 Newton solves (one per timestep), each converging in 1–3
iterations. The nonlinear residual should drop well below the absolute tolerance of 1e-10.

Output files produced:
- `case22_charge_relaxation_out.e` — Exodus mesh with ρ_e and φ fields at every timestep
- `case22_charge_relaxation_out.csv` — avg_rho_e, max_rho_e, max_phi vs. time

---

## Expected Results

### Exponential Decay

The analytical solution is:

```
ρ_e(x, y, t)  = ρ_e(x, y, 0) · exp(−10·t)
max_rho_e(t)  = 1.0 · exp(−10·t)
avg_rho_e(t)  ∝ exp(−10·t)
max_phi(t)    ∝ exp(−10·t)
```

At selected times:

| t (s) | t/τ  | exp(−t/τ) | max_rho_e (approx) |
|--------|------|-----------|---------------------|
| 0.00   | 0    | 1.000     | 1.000               |
| 0.10   | 1    | 0.368     | 0.368               |
| 0.20   | 2    | 0.135     | 0.135               |
| 0.30   | 3    | 0.050     | 0.050               |
| 0.40   | 4    | 0.018     | 0.018               |
| 0.50   | 5    | 0.007     | 0.007               |

The `avg_rho_e` column in the CSV should match this table within numerical roundoff. A
log-linear plot of max_rho_e versus time should produce a perfect straight line with slope
−10 s⁻¹.

### Spatial Fields

The Gaussian shape of ρ_e is preserved throughout the simulation. At t = 0 the blob has
unit amplitude; by t = 0.5 s its amplitude has fallen to approximately 0.007. The spatial
width does not change because the decay rate (σ/ε) is uniform everywhere — every point
decays at the same rate.

The electric potential φ mirrors the charge distribution: a positive dome centred at
(0.5, 0.5) that collapses to zero as the charge dissipates. At t = 0 the potential
maximum is determined by Poisson's equation for the Gaussian source in a grounded box;
it then decays at exactly the same exponential rate as the charge.

### Verification Against Analytical Solution

To verify MOOSE's result, extract `max_rho_e` from the CSV and fit a line to
log(max_rho_e) vs. t. The slope should be −10 s⁻¹ (i.e., −σ/ε). The intercept should
reflect the initial peak of the Gaussian (1.0, but note that `ElementExtremeValue`
evaluates at quadrature points or nodes; the FunctionIC sets nodal values, so the initial
max should be very close to 1.0 for a fine mesh).

---

## Key Takeaways

- **Charge relaxation time τ = ε/σ** is the single parameter controlling how fast free
  charge disappears in an Ohmic medium. It ranges from ~10⁻¹⁹ s (metals) to seconds
  (resistive liquids). Melcher calls it the most fundamental time scale in
  electromechanics.

- **ADReaction** is the correct kernel for a linear volumetric decay term `rate·u`. It
  adds the term to both the residual and the Jacobian automatically via AD, enabling
  Newton convergence in a single iteration for this linear problem.

- **ADCoupledForce** injects one variable as a forcing term in another variable's
  equation. It is the building block for any one-way or two-way coupling between
  different physics fields.

- **ADHeatConduction repurposed as Laplacian**: MOOSE kernels are general mathematical
  operators, not just "heat" or "mechanics" objects. `ADHeatConduction` solves
  −div(k·grad(u)) for any scalar u and any material property k. Renaming the property
  from `thermal_conductivity` to `permittivity` changes the physics label without
  changing the mathematics.

- **Natural boundary condition for ρ_e**: when no BC is specified for a variable, MOOSE
  automatically applies zero normal flux (homogeneous Neumann). This is often the correct
  physical boundary condition and requires no explicit input block.

- **HIT top-level variables** (e.g., `sigma_over_eps = 10.0`) are the idiomatic way to
  expose key physical parameters at the top of an input file. They act as named constants
  that can be referenced with `${name}` anywhere in the file, making the input
  self-documenting and easy to modify.

- **The spatial distribution is frozen**: because the decay rate σ/ε is spatially
  uniform, the charge relaxation ODE decouples in space — every point evolves
  independently with the same time constant. The initial spatial pattern is preserved
  exactly, only its amplitude changes. This is a property of linear equations with
  constant coefficients.
