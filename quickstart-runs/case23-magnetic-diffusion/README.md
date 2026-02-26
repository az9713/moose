# Case 23: Magnetic Diffusion into a Conducting Slab

## Overview

When a magnetic field is suddenly applied to the surface of an electrical conductor,
it does not penetrate instantly. The conductor's finite resistivity slows the field's
advance: eddy currents induced near the surface shield the interior, and the field
can only diffuse inward over time. This phenomenon — **magnetic diffusion** — governs
the response time of transformers, electric motors, induction heaters, electromagnetic
brakes, and magnetic shielding enclosures.

This case models magnetic diffusion into a semi-infinite conducting slab using MOOSE.
The governing equation is:

```
∂B/∂t = D_m · ∂²B/∂x²       magnetic diffusion equation
```

where `D_m = 1/(μ₀ · σ)` is the magnetic diffusivity (m²/s), `μ₀` is the permeability
of free space (4π × 10⁻⁷ H/m), and `σ` is the electrical conductivity (S/m).

This equation is **mathematically identical to the heat equation** (Case 03), with the
substitutions:

| Heat equation | Magnetic diffusion |
|---------------|--------------------|
| Temperature T | Magnetic flux density B |
| Thermal diffusivity α = k/(ρ·cₚ) | Magnetic diffusivity D_m = 1/(μ₀·σ) |
| Heated surface (T=1) | Step-applied field surface (B=1) |
| Insulated far end (T=0) | Field-free far boundary (B=0) |

The MOOSE implementation is therefore a direct reuse of the transient diffusion pattern,
with only the variable name and the physical interpretation changed.

Reference: Melcher, J.R., *Continuum Electromechanics*, MIT Press (1981), Chapter 6,
Sections 6.2–6.3.

---

## The Physics

### The Governing Equation and Its Derivation

Starting from Faraday's law and Ampere's law (in the low-frequency, quasi-static limit
where displacement currents are negligible), combined with Ohm's law for a linear
conductor, Maxwell's equations reduce to the **magnetic diffusion equation**:

```
∂B/∂t = D_m · ∇²B

where   D_m = 1 / (μ₀ · σ)
```

This derivation is as follows. In the quasi-static limit, Faraday's law gives:

```
∇ × E = -∂B/∂t
```

Ohm's law relates the electric field to the current density:

```
J = σ · E
```

Ampere's law (no displacement current) closes the system:

```
∇ × B = μ₀ · J
```

Combining: take the curl of Ampere's law, substitute Faraday's law and Ohm's law,
and use the vector identity `∇ × (∇ × B) = ∇(∇·B) - ∇²B` with `∇·B = 0`:

```
∂B/∂t = (1/(μ₀σ)) · ∇²B = D_m · ∇²B
```

For a 1D slab (variation only in x):

```
∂B/∂t = D_m · ∂²B/∂x²
```

This is exactly the 1D heat equation with `D_m` playing the role of thermal diffusivity.

### Physical Interpretation

Before `t = 0`, the conductor is field-free (`B = 0` everywhere). At `t = 0`, an
external magnetic field of magnitude `B₀ = 1` (normalised) is suddenly applied at the
surface `x = 0`. The field tries to penetrate into the conductor.

The conductor resists this penetration. As the surface field is applied, it drives a
changing flux that, by Faraday's law, induces an electric field in the conductor. This
electric field drives eddy currents (by Ohm's law). The eddy currents, by Ampere's law,
create a magnetic field that opposes the applied field — shielding the interior.

Over time, the eddy currents dissipate energy as Joule heat. The shielding weakens,
and the applied field progressively penetrates deeper into the conductor. The result is
a diffusing front that moves as `δ(t) ≈ 2√(D_m · t)`.

### Analytical Solution

For a semi-infinite slab with a step-applied boundary condition, the exact analytical
solution is the complementary error function:

```
B(x, t) = erfc( x / (2 √(D_m · t)) )
```

where `erfc(z) = 1 - erf(z) = (2/√π) · ∫_z^∞ exp(-u²) du`.

Properties of this solution:
- At `x = 0`: `erfc(0) = 1` — the boundary condition is exactly satisfied.
- As `x → ∞`: `erfc(∞) = 0` — the far-field remains unaffected.
- At any fixed `x`, `B(x,t)` increases monotonically from 0 to 1 as `t → ∞`.
- The "penetration depth" where `B ≈ 0.16` (one e-folding) is `δ ≈ 2√(D_m·t)`.

### Penetration Depth at Key Times

With `D_m = 0.01`:

```
δ(t) = 2 √(D_m · t) = 2 √(0.01 · t) = 0.2 √t
```

| Time | Penetration depth δ | B at midpoint x=0.5 |
|------|---------------------|----------------------|
| t = 2 | δ ≈ 0.28 | erfc(0.5/0.28) ≈ erfc(1.77) ≈ 0.004 |
| t = 5 | δ ≈ 0.45 | erfc(0.5/0.45) ≈ erfc(1.12) ≈ 0.127 |
| t = 10 | δ ≈ 0.63 | erfc(0.5/0.63) ≈ erfc(0.79) ≈ 0.262 |
| t = 20 | δ ≈ 0.89 | erfc(0.5/0.89) ≈ erfc(0.56) ≈ 0.419 |

At `t = 20` the penetration depth is nearly equal to the slab thickness (1 m), so the
right-boundary condition `B = 0` begins to be felt. The simulation captures the full
transition from "thin skin" regime to "thick skin" regime.

### ASCII Domain Diagram

```
  y=0.04  +================================================+
          |                                                |
  (thin   |  B diffuses rightward →  →  →  →  →          |
   strip, |                                                |
  quasi   +================================================+
    1D)   |                                                |
  y=0     +================================================+
          x=0                                           x=1
          B=1                                           B=0
          (step applied)                          (far field)

Time evolution of B(x,t) profiles:
   B
   1 |*
     |**   t=20
     | **
     |  ***  t=10
     |   ***
     |    ****   t=5
     |     ****
     |      *****   t=2
     |       *****
   0 +--+----+----+----+----+-> x
     0  0.2  0.4  0.6  0.8  1.0

Each profile is an erfc curve; the front advances as δ ≈ 2√(D_m·t).
```

### The Skin Effect and D_m for Real Materials

In electrical engineering, the "skin depth" is defined as the depth at which B falls
to `1/e ≈ 37%` of its surface value at a sinusoidal steady state (frequency `f`):

```
δ_skin = √(2 · D_m / (2πf)) = √(1 / (π · f · μ₀ · σ))
```

This is the AC steady-state counterpart of the transient penetration depth. For copper
(`σ = 6×10⁷ S/m`) at 50 Hz: `δ_skin ≈ 9 mm`. At 1 MHz: `δ_skin ≈ 66 μm`.

The transient problem studied here is the DC (step-function) version: a sudden field
application, rather than a sinusoidally varying one.

Magnetic diffusivities for common conductors:

| Material | σ (S/m) | D_m = 1/(μ₀σ) (m²/s) |
|----------|---------|------------------------|
| Copper | 6.0 × 10⁷ | 13.3 × 10⁻³ |
| Aluminium | 3.5 × 10⁷ | 22.7 × 10⁻³ |
| Steel (low-C) | 6.0 × 10⁶ | 133 × 10⁻³ |
| Seawater | 4 | 1.99 × 10⁸ |
| Ionospheric plasma | ~10³ | ~8 × 10² |

Seawater and plasmas have much lower conductivity, so `D_m` is very large — fields
penetrate nearly instantaneously compared to engineering conductors.

---

## Input File Walkthrough

The input file is `case23_magnetic_diffusion.i`. It follows the Case 03 transient
diffusion pattern with two differences: the variable is named `B`, and the kernels use
the `AD` (automatic differentiation) variants for exact Jacobians.

### Header and Top-Level Variable

```
D_m = 0.01
```

A HIT top-level variable declares `D_m` once at the top of the file. The `[Materials]`
block references it as `'${D_m}'`. To change the magnetic diffusivity, edit only this
one line; the change propagates automatically to all blocks that use it.

### `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 50
  ny   = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.04
[]
```

The mesh is a 50×2 grid on a 1 m × 0.04 m rectangle. The physics varies only in `x`
(the diffusion direction), so the 2-element `y` dimension is the minimum needed to
form a valid 2D mesh. The thin aspect ratio (1 : 0.04) makes the domain quasi-1D.

50 elements along `x` gives a spatial resolution of `Δx = 0.02 m`. This is fine enough
to resolve the `erfc` profile even at early times when the front is steep.

### `[Variables]`

```
[Variables]
  [B]
    initial_condition = 0.0
  []
[]
```

`B` is the magnetic flux density normalised by the applied field (`B₀ = 1`). The initial
condition `B = 0` everywhere represents the field-free state of the conductor before
`t = 0`. MOOSE uses first-order Lagrange (linear) basis functions by default.

### `[Kernels]`

```
[Kernels]
  [B_time]
    type     = ADTimeDerivative
    variable = B
  []

  [B_diff]
    type        = ADMatDiffusion
    variable    = B
    diffusivity = diffusivity
  []
[]
```

Two kernels implement the two terms of the magnetic diffusion equation:

**`ADTimeDerivative`** — implements `∂B/∂t` in weak form:

```
R_i = ∫ φ_i · (∂B/∂t) dV
```

The `AD` prefix means MOOSE computes the Jacobian automatically via automatic
differentiation rather than a hand-coded analytic formula. For a linear PDE like this,
the result is identical to using `TimeDerivative`, but the AD variant is more robust for
extensions to nonlinear magnetic materials (μ depending on B).

**`ADMatDiffusion`** — implements `D_m · ∇²B` in weak form:

```
R_i = ∫ D_m · ∇φ_i · ∇B dV
```

The `diffusivity = diffusivity` line means: read the material property named
`'diffusivity'` from the material system. The property name on the right is defined in
the `[Materials]` block. This is the standard MOOSE pattern: kernels do not hardcode
coefficients — they look them up by name, so the material can later be made spatially
varying or temperature-dependent without touching the kernel.

### `[BCs]`

```
[BCs]
  [B_left]
    type     = DirichletBC
    variable = B
    boundary = left
    value    = 1.0
  []

  [B_right]
    type     = DirichletBC
    variable = B
    boundary = right
    value    = 0.0
  []
[]
```

Two Dirichlet boundary conditions:

- **Left boundary** (`x = 0`): `B = 1` — the step-applied external field. This is the
  driving condition; it is applied at `t = 0` and held fixed for all time.
- **Right boundary** (`x = 1`): `B = 0` — the far-field condition. The conductor is
  assumed to extend far enough that the field has not yet reached the right edge. The
  penetration depth at `t = 20` is `δ ≈ 0.89 m`, so the right BC does begin to
  influence the solution near the end of the simulation — this is expected.

The top and bottom boundaries have no BCs specified. Because the diffusion flux
`D_m · ∇B` is zero in the `y` direction (the solution does not vary in `y`), the
natural (zero-flux) condition that MOOSE applies by default is physically correct.

### `[Materials]`

```
[Materials]
  [mag_diff]
    type        = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '${D_m}'
  []
[]
```

`ADGenericConstantMaterial` is the automatic-differentiation version of
`GenericConstantMaterial`. It declares a single material property named `'diffusivity'`
with value `0.01` (substituted from the `D_m` top-level variable at parse time).

The property name `'diffusivity'` is the string that `ADMatDiffusion` looks up. The
name must match exactly. If you want to use a different name (e.g., `'mag_diffusivity'`),
you would write:

```
prop_names  = 'mag_diffusivity'
...
[B_diff]
  diffusivity = mag_diffusivity
[]
```

Using `ADGenericConstantMaterial` (rather than the non-AD version) is required when the
kernels use `AD` variants, so that the automatic differentiation chain is consistent
from material to kernel.

### `[Postprocessors]`

```
[Postprocessors]
  [avg_B]
    type     = ElementAverageValue
    variable = B
  []

  [max_B]
    type       = ElementExtremeValue
    variable   = B
    value_type = max
  []

  [min_B]
    type       = ElementExtremeValue
    variable   = B
    value_type = min
  []
[]
```

Three postprocessors track scalar summaries of the B field at each timestep:

- **`avg_B`**: spatial average of `B` over the slab. Starts at 0 and increases
  monotonically toward 0.5 (the average of a linear profile `B = 1 - x` at steady state).
- **`max_B`**: maximum B in the domain. After the first timestep this equals 1.0 (the
  left BC) and remains at 1.0 throughout, confirming the Dirichlet BC is enforced.
- **`min_B`**: minimum B in the domain. Decreases from the initial 0 as the field front
  advances; the minimum always occurs at the right boundary where `B → 0`.

These three numbers at each timestep are written to the CSV file and printed to the
console, providing a quick sanity check without opening ParaView.

### `[Executioner]`

```
[Executioner]
  type = Transient

  solve_type          = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type = ConstantDT
    dt   = 0.2
  []

  start_time = 0.0
  end_time   = 20.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
```

**`type = Transient`** — time-marching solve. MOOSE steps from `t = 0` to `t = 20`
in 100 steps of `dt = 0.2`.

**`solve_type = 'NEWTON'`** — because the AD variants of the kernels provide exact
Jacobians, the full Newton method (rather than the matrix-free PJFNK) is appropriate.
Newton converges in one iteration per timestep for this linear problem.

**`petsc_options`** — BoomerAMG algebraic multigrid preconditioner from the HYPRE
library. Efficient for the scalar diffusion operator arising from this problem.

**`dt = 0.2`** — the timestep. For the magnetic diffusion equation, the relevant time
scale is the diffusion time `τ = L² / D_m = 1 / 0.01 = 100 s`. The timestep `dt = 0.2`
is 1/500th of `τ`, providing good temporal resolution of the diffusing front.

The spatial Fourier number `Fo = D_m · dt / Δx²` at this step:

```
Fo = 0.01 × 0.2 / (0.02)² = 5.0
```

Backward Euler (the default MOOSE time integrator) is unconditionally stable for all
`Fo > 0`, so large `Fo` is acceptable. It does introduce temporal smoothing, which is
why the erfc front at early times is slightly rounded compared to the exact solution.

### `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- **`exodus = true`**: writes the full `B(x, y, t)` field at every timestep. The Exodus
  file contains 101 time frames (t = 0 through t = 20 at intervals of 0.2). In ParaView,
  animating through these frames shows the diffusion front advancing from left to right.
- **`csv = true`**: writes `avg_B`, `max_B`, and `min_B` vs time to a CSV file, one
  row per timestep.

---

## Running the Simulation

This case requires only the framework's `combined` application (no special physics
modules). It uses standard kernels (`ADTimeDerivative`, `ADMatDiffusion`) available
in the base MOOSE installation.

### Docker (Windows / cross-platform)

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case23-magnetic-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case23_magnetic_diffusion.i 2>&1 | tail -30'
```

### Native Linux / macOS

```bash
cd quickstart-runs/case23-magnetic-diffusion
combined-opt -i case23_magnetic_diffusion.i
```

### Parallel Execution

```bash
mpirun -n 4 combined-opt -i case23_magnetic_diffusion.i
```

The simulation runs 100 timesteps, each requiring one Newton iteration. Total wall time
is typically a few seconds on a single core.

Output files produced:
- `case23_magnetic_diffusion_out.e` — Exodus file with B(x,y,t) at all 101 time frames
- `case23_magnetic_diffusion_out.csv` — CSV with avg_B, max_B, min_B vs time

---

## Expected Results

### Postprocessor Time History

The CSV file will have this structure (approximate values):

```
time,avg_B,max_B,min_B
0,0,0,0
0.2,<small>,1.0,0
0.4,...,1.0,0
...
5.0,~0.108,1.0,0
10.0,~0.211,1.0,0
20.0,~0.384,1.0,0
```

Key observations:
- `max_B` locks to `1.0` after the first timestep and stays there (left BC enforced).
- `min_B` stays at `0.0` throughout (right BC enforced or far field still unperturbed).
- `avg_B` increases monotonically, measuring how much total flux has entered the slab.

### Analytical Comparison: avg_B(t)

The analytical average of the erfc profile over `x ∈ [0, 1]`:

```
avg_B(t) = ∫₀¹ erfc(x / (2√(D_m·t))) dx
```

For `D_m = 0.01`, this integral evaluates (numerically) to:

| Time | Analytical avg_B | MOOSE expected |
|------|-----------------|----------------|
| t = 2 | 0.045 | ~0.044 |
| t = 5 | 0.108 | ~0.107 |
| t = 10 | 0.211 | ~0.209 |
| t = 20 | 0.384 | ~0.381 |

The MOOSE result will underestimate slightly because backward Euler introduces temporal
smearing of the diffusion front (the exact erfc is more peaked). The discrepancy shrinks
with smaller `dt`.

### B Profile Shape at Selected Times

At each time, the B(x) profile follows the erfc curve. Extracting a line plot along
`y = 0.02` (the midline of the strip) from the Exodus file in ParaView shows:

```
t = 5 (δ ≈ 0.45):
  B
  1 |*
    |**
    | **
    |  ***
    |    *****
    |        ***
  0 +--+--+--+--+---> x
    0 .2 .4 .6 .8  1

t = 20 (δ ≈ 0.89):
  B
  1 |*
    |**
    | ***
    |   ****
    |      *****
    |           ***
  0 +--+--+--+--+---> x
    0 .2 .4 .6 .8  1
```

The front becomes less steep and more elongated as time increases. Near `t = 25` the
penetration depth would reach 1 m (the slab thickness), and the boundary condition
at `x = 1` would begin to appreciably alter the profile — this is visible at `t = 20`
as the right portion of the profile deviates slightly below the infinite-slab erfc.

---

## Analytical Verification

To verify the MOOSE solution against the exact erfc formula, extract a line plot from
the Exodus file at any timestep and compare:

```python
import numpy as np
from scipy.special import erfc

D_m = 0.01
t   = 5.0   # check at t=5

x  = np.linspace(0, 1, 200)
B_exact = erfc(x / (2 * np.sqrt(D_m * t)))

# Compare against values extracted from MOOSE CSV or Exodus
```

Alternatively, compare the postprocessor `avg_B` at `t = 5` against the analytical
integral:

```python
from scipy.integrate import quad
val, _ = quad(lambda x: erfc(x / (2 * np.sqrt(D_m * t))), 0, 1)
print(f"Analytical avg_B at t={t}: {val:.4f}")
# Expected: ~0.1081
```

---

## Key Takeaways

- **Mathematical analogy**: the magnetic diffusion equation is identical in form to
  the heat equation. Any MOOSE transient diffusion setup can model magnetic field
  penetration simply by renaming the variable and reinterpreting `D_m`.

- **Magnetic diffusivity** `D_m = 1/(μ₀σ)` plays the role of thermal diffusivity. High
  conductivity means small `D_m` — a good conductor shields its interior longer.

- **ADTimeDerivative and ADMatDiffusion**: the AD (automatic differentiation) kernel
  variants provide exact Jacobians without hand-coding. They are preferred for
  extensibility: if the conductivity later depends on temperature or B, the Jacobian
  remains exact with no extra work.

- **ADGenericConstantMaterial**: the AD-compatible material for constant properties.
  Must match the kernel variant (AD kernels require AD materials for consistency of
  the derivative chain).

- **HIT top-level variables** (`D_m = 0.01`): centralise parameters at the top of the
  file; reference them as `'${D_m}'` in sub-blocks. One edit changes everything.

- **Penetration depth** `δ ≈ 2√(D_m·t)`: the square-root-of-time scaling is universal
  for all diffusion problems. It quantifies how quickly information propagates.

- **erfc analytical solution**: the complementary error function is the exact solution
  for the semi-infinite slab with a step boundary condition. It is an important benchmark
  for validating any diffusion solver.

- **Skin effect in engineering**: at AC frequencies, the penetration depth becomes the
  "skin depth" `δ_skin = √(D_m/(πf))`. Magnetic diffusion is the time-domain equivalent.
  Understanding the transient problem builds intuition for transformer core losses,
  induction heating efficiency, and electromagnetic shielding performance.
