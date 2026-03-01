# Case 101: RC Transmission Line Transient — Voltage Diffusion on a Lossy Cable

## Overview

This case models the transient voltage response of an RC transmission line — a lossy cable where series resistance and shunt capacitance dominate over inductance. The governing equation is a pure diffusion equation for voltage rather than the full wave equation, and it describes signal propagation on submarine telegraph cables, VLSI interconnects at low frequencies, and biological axons (the Hodgkin-Huxley cable equation in its passive limit).

The telegrapher's equations for a distributed R-C line (inductance neglected) collapse to:

```
dv/dt = D * d²v/dx²,    D = 1 / (R * C)
```

where R [Ohm/m] and C [F/m] are the per-unit-length resistance and capacitance. A step voltage V₀ = 1 V is applied at x = 0; the far end (x = L = 1 m) is open-circuit (zero-flux Neumann). The semi-infinite analytic solution is:

```
v(x, t) = V0 * erfc(x / (2 * sqrt(D * t)))
```

and the characteristic diffusion delay to position x is t_delay = x² / (4D). This is MIT 6.641 Lecture 17 material, demonstrating that signal propagation on a resistive cable is diffusive — not wave-like — with a delay that grows quadratically with distance rather than linearly.

MOOSE solves this with a single variable `v`, the `ADTimeDerivative` kernel for the time term, and `ADMatDiffusion` for the spatial diffusion operator. The use of AD (automatic differentiation) variants provides exact Jacobians without symbolic derivation.

---

## The Physics

**Governing equation:**

```
dv/dt = D * d²v/dx²,    D = 1 / (R * C) = 0.1 m²/s
```

This is mathematically identical to the heat equation with temperature replaced by voltage, thermal diffusivity replaced by 1/(RC), and heat flux replaced by the displacement current density.

**Boundary conditions:**

| Boundary | Location | Condition | Value |
|----------|----------|-----------|-------|
| Source end | x = 0 | Dirichlet (ramped) | v rises from 0 to V₀ = 1 V over 0.01 s |
| Load end | x = 1 | Natural Neumann | dv/dx = 0 (open circuit, zero current) |

The small ramp at t = 0 (PiecewiseLinear from 0 to V₀ over 0.01 s) avoids a discontinuous initial condition that would otherwise require an extremely fine time step to resolve.

**Material properties:**

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| RC diffusivity | D | 0.1 | m²/s |
| Step voltage | V₀ | 1.0 | V |
| Cable length | L | 1.0 | m |

**Domain geometry:**

- Quasi-1D domain: [0, 1] × [0, 0.04], 100 × 2 QUAD4 elements
- Resolution: 100 elements along the cable axis, 0.01 m spacing
- Physical time range: 0 to 5 s (two diffusion time scales to x = 1)

**Analytic diffusion timescales:**

| Position | t_delay = x² / (4D) |
|----------|---------------------|
| x = 0.25 | 0.156 s |
| x = 0.50 | 0.625 s |
| x = 0.75 | 1.406 s |
| x = 1.00 | 2.500 s |

---

## Input File Walkthrough

### Mesh

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100
  ny   = 2
  xmin = 0
  xmax = 1.0
  ymin = 0
  ymax = 0.04
[]
```

The cable is modelled as a thin 2D strip. With only 2 elements in y, the geometry is essentially 1D — MOOSE uses standard 2D QUAD4 elements but the solution varies only in x. This avoids the corner-node complications of a 1D mesh while still using the standard 2D element library.

### Variables

```
[Variables]
  [v]
    initial_condition = 0
  []
[]
```

The voltage variable `v` starts at zero everywhere, representing an initially uncharged cable.

### Kernels

```
[Kernels]
  [v_time]
    type     = ADTimeDerivative
    variable = v
  []
  [v_diffusion]
    type        = ADMatDiffusion
    variable    = v
    diffusivity = diffusivity
  []
[]
```

`ADTimeDerivative` contributes the dv/dt term. `ADMatDiffusion` contributes -D nabla²v in weak form — in 1D this is exactly the right-hand side diffusion operator. Both use the AD (automatic differentiation) framework, so the Jacobian is computed exactly via forward-mode AD rather than finite differences.

### Boundary Conditions and Functions

```
[BCs]
  [source_step]
    type     = FunctionDirichletBC
    variable = v
    boundary = left
    function = step_voltage
  []
[]

[Functions]
  [step_voltage]
    type = PiecewiseLinear
    x    = '0      0.01    100'
    y    = '0      ${V0}   ${V0}'
  []
[]
```

`PiecewiseLinear` linearly interpolates between the x-y pairs. The voltage rises from 0 at t = 0 to V₀ = 1 at t = 0.01 s, then holds constant. The 100 s endpoint ensures V₀ is maintained for the full 5 s simulation. No BC is applied to the right boundary, so MOOSE automatically uses the natural (zero-flux) Neumann condition, representing an open-circuit load.

### Materials

```
[Materials]
  [rc_cable]
    type        = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '${D_rc}'
  []
[]
```

A single constant material property provides D = 0.1 m²/s uniformly across the domain. Using the AD variant (`ADGenericConstantMaterial`) is required because the kernel `ADMatDiffusion` expects AD material properties.

### Postprocessors

Five postprocessors are defined:

| Postprocessor | Location | Purpose |
|--------------|----------|---------|
| `v_quarter` | x = 0.25 | Tracks early voltage arrival |
| `v_mid` | x = 0.50 | Midpoint voltage history |
| `v_three_quarter` | x = 0.75 | Late voltage arrival |
| `v_load` | x = 1.00 | Open-circuit end |
| `avg_v` | domain-wide | Total charge fraction |

The staggered arrival times at these four positions directly demonstrate the quadratic delay law: the midpoint sees significant voltage at t ≈ 0.625 s while the far end does not reach half V₀ until t ≈ 2.5 s.

### Executioner

```
[Executioner]
  type     = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  dt       = 0.05
  end_time = 5.0
[]
```

One hundred uniform time steps of dt = 0.05 s. The diffusion stability criterion for explicit methods is dt < dx² / (2D) = 0.01² / 0.2 = 0.0005 s; since MOOSE uses an implicit method (backward Euler via `ADTimeDerivative`), the much larger dt = 0.05 s is stable. BoomerAMG is efficient for the symmetric positive-definite system produced by the Laplacian.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case101-transmission-line \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case101_transmission_line.i 2>&1 | tail -30'
```

---

## Expected Results

**Voltage front progression:**

The erfc wave front travels from left to right. The position x at which v = V₀/2 at time t is:

```
x_half = 2 * sqrt(D * t) * erfc_inv(0.5) ≈ 0.954 * sqrt(D * t)
```

So the "half-voltage front" reaches x = 0.25 m at t ≈ 0.069 s, x = 0.5 m at t ≈ 0.276 s, x = 0.75 m at t ≈ 0.621 s, and x = 1.0 m at t ≈ 1.103 s.

**Steady-state (t >> L²/D):**

At long times the voltage equalises to V₀ = 1 V everywhere along the cable. The average voltage `avg_v` approaches 1 V as the cable fully charges.

**Open-circuit end:**

Because dv/dx = 0 at x = 1, the voltage at the load end approaches V₀ from below, reaching ~0.5 V at t ≈ 1.1 s and ~0.9 V at t ≈ 2 s. The open end sees voltage last and charges slowest.

**CSV output columns:**

The postprocessor CSV file `case101_transmission_line_out.csv` contains time plus five voltage columns. Plotting `v_quarter`, `v_mid`, `v_three_quarter`, and `v_load` versus time on the same axes clearly shows the staggered, increasingly delayed arrival of the voltage front at successive positions — the signature of diffusive rather than wave-like propagation.

**Comparison to wave propagation:**

On an LC transmission line (Case 100 analogy in electromagnetics), a step voltage would arrive at x = L at time t = L / v_p with a sharp front. On the RC cable here, there is no sharp front: every point in the cable sees an instantaneous (but tiny) voltage response as soon as the step is applied, growing slowly as erfc. This is the key physical difference between inductive (wave) and resistive (diffusive) transmission lines.

---

## Key Takeaways

- **Diffusion vs. wave propagation**: The RC transmission line model produces voltage diffusion, not wave propagation. Signal delay scales as x², not x — a fundamental limitation of resistive interconnects at low frequencies.
- **ADTimeDerivative + ADMatDiffusion**: The AD kernel pair provides the complete diffusion equation with exact Jacobians via automatic differentiation, requiring no hand-coded residual or Jacobian terms.
- **PiecewiseLinear ramp BC**: A short linear ramp at t = 0 avoids the initial discontinuity that would arise from a perfect step, improving Newton convergence in the first time steps without affecting long-time behaviour.
- **Natural Neumann at open end**: Omitting a BC on the right boundary automatically enforces dv/dx = 0 — the natural condition for the variational statement of the diffusion equation, physically corresponding to zero current at the open-circuit load.
- **Quadratic delay law**: The delay time t ~ x²/(4D) means doubling the cable length quadruples the delay — a critical consideration for submarine cable design (the original motivation for Kelvin's analysis) and VLSI wire timing.
- **Connection to MIT 6.641 Lecture 17**: This is Zahn's RC cable diffusion model, part of a broader treatment of transmission lines that also covers the full TEM wave mode in LC lines and the transition between the two regimes.
