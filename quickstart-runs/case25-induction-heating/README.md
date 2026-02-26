# Case 25: Induction Heating — Magnetic Diffusion Generates Heat

## Overview

Induction heating is one of the signature applications of continuum electromechanics.
An alternating magnetic field is applied to the surface of a conducting slab. The
changing field drives eddy currents inside the conductor according to Faraday's law,
and those currents dissipate energy as Joule heat — but only in the thin surface layer
called the **electromagnetic skin depth**. This is why induction cooktops heat only
the metal pot, not the air or the glass-ceramic surface above the coil, and why
induction furnaces can selectively heat just the outer rim of a rotating workpiece.

The governing physics is developed in Chapter 6 of Melcher's *Continuum
Electromechanics* (MIT Press, 1981), specifically §6.7–6.8. This case builds
directly on:

- **Case 23** (Magnetic Diffusion): the single-field `∂B/∂t = D_m·∇²B` problem,
  which established the penetration depth concept.
- **Case 17** (Joule Heating): the coupled electro-thermal problem using `ADHeatConduction`
  and `ADHeatConductionTimeDerivative`.

The new element introduced here is the bridge between those two physics: extracting
the spatial gradient of one field (B) and using it to drive a source term in another
field (T). This pattern — compute gradient, square it, inject as a body force — recurs
throughout electromechanics and MHD.

---

## The Physics

### Governing Equations

The problem is formulated on a quasi-1D slab with the applied field oscillating at
the left surface (x = 0) and the right surface (x = 1) held at zero field:

```
∂B/∂t = D_m · ∂²B/∂x²                    (1)   magnetic diffusion

ρcp · ∂T/∂t = k · ∂²T/∂x² + Q_eddy      (2)   heat equation

Q_eddy = Q_coeff · (∂B/∂x)²              (3)   eddy-current heating
```

**Equation (1)** is the magnetic diffusion equation. It is mathematically identical
to the heat equation (Case 03), with the magnetic diffusivity D_m = 1/(μ₀σ) playing
the role of thermal diffusivity.

**Equation (2)** is the Fourier heat equation with a volumetric source term. The
source Q_eddy is the Joule dissipation per unit volume due to the eddy currents.

**Equation (3)** derives from Ampere's law in the low-frequency (magnetoquasistatic)
limit. The eddy current density is J = (1/μ₀)·∂B/∂x. The Joule power density is
Q = J²/σ = (1/μ₀²σ)·(∂B/∂x)². Collecting the constants into Q_coeff:

```
Q_coeff = 1/(μ₀²·σ) = D_m/μ₀
```

In normalised units this is set to 0.5 — a modest heating coefficient chosen so
that the temperature rise over 10 oscillation periods is physically plausible and
numerically well-conditioned.

### The Skin Depth

For a sinusoidally oscillating applied field at angular frequency ω = 2π/τ, the
steady-state B-field solution (Melcher §6.2) shows exponential penetration from the
surface with characteristic length:

```
δ = √(2·D_m / ω)
```

With D_m = 0.005 and τ = 0.5 s (ω = 4π ≈ 12.57 rad/s):

```
δ = √(2 × 0.005 / 4π) = √(0.01 / 4π) ≈ 0.028
```

Since the slab is 1 m thick and δ ≈ 0.028 m, the field is essentially zero for
x > 4δ ≈ 0.11 m. The 50-element mesh (element size 0.02 m) places approximately
1.4 elements within one skin depth — adequate for qualitative behaviour. A finer
mesh (100 elements) would give better spatial resolution of the heating peak.

The eddy-current heating Q ~ (∂B/∂x)² also peaks near the surface. The gradient of
an exponential e^(−x/δ)·sin(·) is proportional to e^(−x/δ)/δ, so Q ~ e^(−2x/δ)/δ².
The heating is therefore strongly concentrated at x ≲ δ. Interior regions (x ≫ δ)
see negligible heating and remain near the initial 300 K.

### Boundary Conditions

| Variable | Boundary | Type      | Value              | Physical meaning               |
|----------|----------|-----------|--------------------|-------------------------------|
| B        | left     | Dirichlet | sin(2πt/τ)         | Oscillating applied field      |
| B        | right    | Dirichlet | 0                  | Far-field: field decays to zero |
| T        | left     | Dirichlet | 300 K              | Left surface as a thermal sink |
| T        | right    | Dirichlet | 300 K              | Right surface as a thermal sink |

The top and bottom boundaries have no explicit BC, so the natural (zero-flux Neumann)
condition is applied automatically. This is correct for the quasi-1D geometry.

### Parameters

| Symbol    | Name                      | Value | Units        |
|-----------|---------------------------|-------|--------------|
| D_m       | Magnetic diffusivity      | 0.005 | m²/s         |
| k         | Thermal conductivity      | 1.0   | W/(m·K)      |
| ρcp       | Volumetric heat capacity  | 1.0   | J/(m³·K)     |
| Q_coeff   | Eddy-heating coefficient  | 0.5   | (normalised) |
| τ         | Oscillation period        | 0.5   | s            |
| δ         | Skin depth                | ~0.028| m            |
| L         | Slab thickness            | 1.0   | m            |

---

## Input File Walkthrough

### HIT Top-Level Variables

```
D_m     = 0.005
k_th    = 1.0
rho_cp  = 1.0
Q_coeff = 0.5
tau     = 0.5
```

All five key physical parameters are declared at the top of the file and referenced
with `${name}` throughout. Changing D_m automatically updates both the magnetic
diffusion material property and the skin depth physics. Changing tau updates both
the applied-field BC and the comment in the header.

### `[Mesh]`

A 50×2 uniform quadrilateral mesh over the quasi-1D domain [0,1]×[0,0.04]. The
aspect ratio (25:1 per element) is large but acceptable for a quasi-1D problem where
the solution is y-independent. The 50 elements in x give an element size of 0.02 m,
placing about 1.4 elements per skin depth.

### `[Variables]`

Two nodal (FIRST-order LAGRANGE) variables:

```
B   --  magnetic flux density [T, normalised]
         IC: 0 everywhere (field-free conductor at t=0)

T   --  temperature [K]
         IC: 300 K everywhere (uniform reference state)
```

### `[AuxVariables]` and `[AuxKernels]`

This is the key new pattern for this case. The eddy-current heating source Q_eddy
depends on (∂B/∂x)², which is a nonlinear function of the gradient of a primary
variable. MOOSE does not allow kernel residuals to reference gradients of other
variables directly through the standard coupling mechanism. The workaround uses two
AuxVariable / AuxKernel pairs:

**Step 1 — Extract the gradient:**

```
[dBdx]           order = CONSTANT, family = MONOMIAL
[calc_dBdx]      type = VariableGradientComponent
                 gradient_variable = B
                 component = x
                 execute_on = 'TIMESTEP_END'
```

`VariableGradientComponent` computes the x-component of grad(B) at element centres
and stores it in the CONSTANT MONOMIAL AuxVariable `dBdx`. The `CONSTANT MONOMIAL`
space is used because gradients of LAGRANGE variables are piecewise-constant on each
element. `execute_on = 'TIMESTEP_END'` means dBdx reflects the converged B field at
the end of each timestep, ready to be used as a source at the next step.

**Step 2 — Compute the heating:**

```
[eddy_heat]      order = CONSTANT, family = MONOMIAL
[calc_eddy_heat] type = ParsedAux
                 coupled_variables = 'dBdx'
                 expression = '${Q_coeff} * dBdx * dBdx'
                 execute_on = 'TIMESTEP_END'
```

`ParsedAux` evaluates an algebraic expression involving AuxVariables and stores the
result in another AuxVariable. This avoids the need to write a custom C++ material
or AuxKernel for a simple algebraic operation.

### `[Kernels]`

Five kernels define the two coupled PDEs:

**Magnetic diffusion (variable = B):**

```
[B_time]   ADTimeDerivative    -- ∂B/∂t
[B_diff]   ADMatDiffusion      -- D_m·∇²B  (reads 'mag_diffusivity' material property)
```

**Heat equation (variable = T):**

```
[T_time]   ADHeatConductionTimeDerivative  -- ρcp·∂T/∂t
[T_diff]   ADHeatConduction               -- k·∇²T
[T_source] CoupledForce                   -- +Q_eddy (from eddy_heat AuxVariable)
```

`CoupledForce` adds +∫ v·ψ_i dV to the residual for T, where v = eddy_heat. The
positive sign is correct: heating adds energy to the system (positive source on the
right-hand side of the PDE). Note that `CoupledForce` (without the AD prefix) is
used here because eddy_heat is an AuxVariable, not a primary variable, and non-AD
coupling to AuxVariables is standard practice.

### Lagged Coupling Pattern

The B–T coupling is **one-way and lagged by one timestep**:

```
timestep n:   solve B_n (using T-independent equations)
              solve T_n (using eddy_heat from step n-1 as source)
TIMESTEP_END: update dBdx and eddy_heat from B_n
timestep n+1: solve B_{n+1}, T_{n+1} (now using eddy_heat from step n)
```

This is not a fully implicit coupled solve. The lag error is O(dt) — for dt = 0.01 s
and slowly varying Q_eddy, this is negligible. For problems where the heating changes
rapidly the timestep should be refined or a fully implicit scheme should be used.

### `[BCs]`

The oscillating applied field uses `ADFunctionDirichletBC` with an inline function
expression. The `${tau}` substitution inserts the period at parse time:

```
[B_left]  ADFunctionDirichletBC  function = 'sin(2*pi*t/${tau})'
```

MOOSE's built-in parser recognises `pi` as the mathematical constant. The factor
`2*pi/0.5 = 4π` sets the angular frequency ω.

### `[Materials]`

Two AD-compatible material blocks:

```
[mag_diff]   ADGenericConstantMaterial   mag_diffusivity = D_m
[thermal]    ADGenericConstantMaterial   thermal_conductivity, specific_heat, density
```

`ADHeatConductionTimeDerivative` requires both `density` and `specific_heat` as
separate properties — it forms ρ·cp internally. The product ρcp = density ×
specific_heat = 1.0 × 1.0 = 1.0, matching the governing equation.

### `[Executioner]`

```
type       = Transient
solve_type = NEWTON
dt         = 0.01     (50 steps per oscillation period)
end_time   = 5.0      (10 full oscillation periods)
```

With 50 timesteps per period the oscillating boundary condition is well-resolved
temporally. Ten periods allow the transient heating to reach a quasi-periodic
steady state where the temperature profile oscillates about a time-averaged shape.

---

## Running the Simulation

This case requires the `combined` application, which provides `ADHeatConduction`,
`ADHeatConductionTimeDerivative`, `VariableGradientComponent`, and `ParsedAux`.
The magnetic diffusion kernels (`ADTimeDerivative`, `ADMatDiffusion`, `CoupledForce`)
are part of the base MOOSE framework.

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case25-induction-heating \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case25_induction_heating.i 2>&1 | tail -30'
```

Expected console output: 500 Newton solves (one per timestep), each converging in
1–4 iterations. The B equation (linear in B) converges in 1–2 iterations. The T
equation (linear in T for fixed eddy_heat) also converges in 1–2 iterations. The
combined nonlinear residual should drop well below 1e-8 at every step.

Output files produced:

- `case25_induction_heating_out.e` — Exodus mesh with B, T, dBdx, eddy_heat fields
  at every timestep, suitable for animation in ParaView.
- `case25_induction_heating_out.csv` — max_T, avg_T, max_B, max_eddy_heat vs. time,
  suitable for plotting with Python or a spreadsheet.

---

## Expected Results

### Magnetic Field

The B field oscillates sinusoidally at x = 0 with unit amplitude. It decays
exponentially with x at the skin depth δ ≈ 0.028 m, so by x = 0.1 m the amplitude
is below exp(−0.1/0.028) ≈ 0.03 — less than 3% of the surface value. The interior
(x > 0.1 m) is effectively magnetically shielded.

### Eddy-Current Heating

The heating power density Q_eddy = 0.5·(∂B/∂x)² peaks just inside the left surface.
Its peak value oscillates with the applied field: maximum when |∂B/∂x| is largest
(near the zero-crossing of sin(2πt/τ), where the time derivative is maximum) and
zero twice per period. The time-averaged heating profile decays as e^(−2x/δ).

### Temperature Field

At early times (t < τ) the temperature is nearly uniform at 300 K. As oscillations
continue, heat accumulates in the skin-depth layer. The steady periodic state shows:

- Near-surface (x ≲ 0.05 m): elevated temperature, rising with time until conduction
  to the left boundary (held at 300 K) balances the heating.
- Interior (x > 0.1 m): temperature remains near 300 K throughout.

Approximate max_T values at selected times:

| t (s) | Oscillations | max_T (approx) |
|--------|-------------|----------------|
| 0.5    | 1           | ~302 K         |
| 1.0    | 2           | ~304 K         |
| 2.0    | 4           | ~307 K         |
| 5.0    | 10          | ~314 K         |

The exact values depend on Q_coeff. With Q_coeff = 0.5 the temperature rise is
modest — roughly 1 K per two oscillation periods — reflecting the moderate heating
coefficient in normalised units.

### CSV Postprocessors

The `max_eddy_heat` column will oscillate roughly as sin²(2πt/τ), peaking at 0, τ/2,
τ, etc. (when |∂B/∂x| is maximum at the surface). The `max_T` and `avg_T` columns
should increase monotonically with a superimposed small oscillation at frequency 2ω
(twice the field frequency, because heating goes as the square of B).

---

## Key Takeaways

- **VariableGradientComponent** is the standard AuxKernel for extracting one spatial
  component of a primary variable's gradient into a CONSTANT MONOMIAL AuxVariable.
  Gradients of LAGRANGE variables are piecewise-constant on elements, making CONSTANT
  MONOMIAL the correct discretisation space.

- **ParsedAux** evaluates any algebraic expression involving coupled AuxVariables at
  element centres without requiring custom C++ code. It is the go-to tool for derived
  quantities such as squared gradients, products of fields, or any quantity that can
  be written as a mathematical formula.

- **Lagged one-way coupling via AuxVariables** is a practical pattern for multiphysics
  problems where the coupling is not stiff: compute a derived quantity at TIMESTEP_END
  from converged primary fields, then use it as a forcing term in the next timestep.
  The lag error is first-order in dt and is acceptable when dt is small relative to
  the physical time scale of the coupling.

- **CoupledForce** adds a variable (here an AuxVariable) as a volumetric source to
  another variable's equation. It is the correct kernel when the source is represented
  as a MOOSE variable rather than a material property or analytic function.

- **Skin depth concentration**: the eddy-current heating is not uniform but exponentially
  concentrated near the surface within one skin depth δ = √(2D_m/ω). This is the
  physical basis of induction hardening, surface heat treatment, and electromagnetic
  shielding. Increasing ω (higher frequency) or decreasing D_m (better conductor)
  makes the heating more surface-localised.

- **Magnetic vs. thermal time scales**: the magnetic diffusion time scale τ_m ~ L²/D_m
  (here ~200 s for L=1) is much longer than the oscillation period τ = 0.5 s, which
  is why the quasi-static oscillatory skin-depth picture holds. The thermal diffusion
  time scale τ_th ~ L²·ρcp/k = 1 s means the temperature responds on the same order
  as the field oscillations, giving the transient temperature build-up observed here.

- **Modular physics composition**: this case demonstrates MOOSE's strength for
  multi-physics problems. The magnetic diffusion equation (Case 23) and the heat
  equation (Case 17) each use their own well-tested kernels and materials. The only
  new code required is the AuxKernel pair that bridges them — no custom C++ kernel
  was needed.
