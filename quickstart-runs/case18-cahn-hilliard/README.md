# Case 18: Cahn-Hilliard Spinodal Decomposition

## Overview

This case models **spinodal decomposition** — the spontaneous unmixing of a binary
mixture into two distinct phases. Start with a nearly-uniform mixture (concentration
c = 0.5 everywhere, plus tiny random noise) and watch it separate: regions rich in
component A (c near 1) and regions rich in component B (c near 0) grow and coarsen
over time until large, clearly-separated domains fill the domain.

This is the governing physics behind:
- Metallic alloy decomposition (e.g., Fe-Cr steel aging)
- Polymer blend morphology during processing
- Oil-water separation in emulsions
- Pattern formation in biological systems

The simulation uses the **Cahn-Hilliard equation**, a fourth-order PDE that conserves
total solute (the average concentration stays at 0.5 for all time) while driving the
system toward lower free energy. MOOSE solves it in **split form** — two coupled
second-order equations — which avoids the need for special high-continuity elements.

---

## The Physics

### Free Energy and Why Phases Separate

Every binary mixture has a **free energy** F(c). For a system that undergoes spinodal
decomposition, F has a characteristic **double-well shape** with two minima — one near
c = 0 and one near c = 1. This case uses:

    F(c) = c^2 * (1 - c)^2

This polynomial has minima at c = 0 and c = 1 and a local maximum (energy barrier) at
c = 0.5. A uniform mixture at c = 0.5 sits at the top of the barrier: it is unstable.
Any small fluctuation in concentration lowers the local free energy and grows
spontaneously. This is **spinodal decomposition** — there is no nucleation barrier;
perturbations of every wavelength are thermodynamically downhill from the start.

### The Cahn-Hilliard Equation

The full Cahn-Hilliard equation is:

    dc/dt = div( M * grad( dF/dc - kappa * laplacian(c) ) )

where:

- `c`       — concentration (ranges from 0 to 1)
- `t`       — time
- `M`       — atomic mobility (controls kinetics, set to 1.0)
- `F(c)`    — bulk free energy density (double-well)
- `dF/dc`   — bulk chemical potential (thermodynamic driving force)
- `kappa`   — gradient energy coefficient (set to 1.0; penalises sharp interfaces)
- `laplacian(c)` — the Laplacian of c

The term `dF/dc - kappa * laplacian(c)` is the **generalized chemical potential** w.
The Laplacian term adds an energetic cost to spatial gradients of concentration, which
regularises the interface and sets the equilibrium interface width proportional to
sqrt(kappa / |d2F/dc2|).

### Why the Split Form?

The Cahn-Hilliard equation is fourth-order in space (it contains a biharmonic term
`laplacian(laplacian(c))`). Standard C0 Lagrange finite elements only guarantee
continuity of the function, not its gradient, so the biharmonic term cannot be
assembled directly. The **split form** rewrites the problem as two coupled second-order
equations:

    dc/dt = div( M * grad(w) )          [Equation 1 — on variable w]
    w     = dF/dc - kappa * laplacian(c) [Equation 2 — on variable c]

Now each equation is second-order and can be assembled with standard linear (FIRST
LAGRANGE) elements. The cost is one additional degree of freedom per node (w), but
this is much cheaper than using Hermite or other C1 elements.

### Conservation of Mass

The Cahn-Hilliard equation is a **conservative** PDE. In the absence of sources and
with periodic or no-flux boundary conditions, the total amount of each component is
preserved:

    d/dt [ integral(c) dV ] = 0

This is a direct consequence of the divergence form of Equation 1: the flux `M*grad(w)`
integrates to zero over a periodic domain. You can verify this by watching the
`avg_c` postprocessor, which should remain at 0.5 for the entire simulation.

### Coarsening Dynamics

After the initial spinodal decomposition (small domains form quickly), the system
enters a **coarsening** regime. Smaller domains have more surface area per unit volume
and therefore higher interfacial energy (from the kappa term). They are absorbed by
larger domains, reducing total interface length and lowering the total free energy.
This coarsening follows the Lifshitz-Slyozov power law:

    L(t) ~ t^(1/3)

where L is the characteristic domain size. You will observe this as the plot evolves:
the pattern becomes coarser and coarser, but the total number of distinct regions
decreases monotonically.

---

## Input File Walkthrough

### `[Mesh]` Block

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 40
  ny   = 40
  xmin = 0
  xmax = 25
  ymin = 0
  ymax = 25
[]
```

A 40x40 element mesh on a 25x25 domain gives element size h = 0.625. The interface
width is set by kappa and the free energy curvature; with kappa=1 and the double-well
F, the interface half-width is approximately sqrt(2*kappa / |d2F/dc2|_max) ~ 1.4
length units. The mesh resolves this with roughly two elements per interface half-width,
which is adequate for visualisation-quality results.

### `[Variables]` Block

Two variables are solved simultaneously:

**`c`** (concentration): the primary order parameter. Initialized with `RandomIC` to
c = 0.5 +/- 0.06. This small noise seeds the instability; without it, the perfectly
uniform c=0.5 state would have zero driving force (by symmetry of the residual) and
nothing would happen.

**`w`** (chemical potential): the split auxiliary variable. Starts at zero (default
InitialCondition); it is immediately updated by the solver on the first Newton
iteration.

### `[Kernels]` Block

Four kernels implement the two coupled equations:

**`c_dot`** (`CoupledTimeDerivative` on w): Adds the `(dc/dt, phi_w)` term to
Equation 1. Note the unusual structure: `CoupledTimeDerivative` operates on variable
`w` but involves the time derivative of coupled variable `c`. This correctly assembles
the mass-matrix term in the w equation.

**`w_res`** (`SplitCHWRes` on w): Adds `(M * grad(w), grad(phi_w))`. This is the
diffusion of chemical potential — the flux term that moves concentration around.

**`c_res`** (`SplitCHParsed` on c): Adds the chemical potential equation. It reads
`dF/dc` from the `DerivativeParsedMaterial` and assembles `(w - dF/dc + kappa*laplacian(c), phi_c) = 0`.
The Jacobian is also assembled automatically using the second derivative `d2F/dc2`.

### `[Materials]` Block

**`free_energy`** (`DerivativeParsedMaterial`): Parses the expression `c^2*(1-c)^2`
and automatically differentiates it to produce dF/dc and d2F/dc2 as named material
properties. `SplitCHParsed` requests these derivatives by name at every quadrature
point. Setting `derivative_order = 2` instructs the material to precompute up to the
second derivative.

**`const`** (`GenericConstantMaterial`): Provides `kappa_c = 1.0` (gradient energy
penalty) and `M = 1.0` (mobility) as uniform constant material properties.

### `[Executioner]` Block

**`solve_type = NEWTON`**: Full Newton with an analytically assembled Jacobian. The
DerivativeParsedMaterial provides the Jacobian contributions automatically.

**`scheme = bdf2`**: Second-order backward differentiation formula time integration.
BDF2 is more accurate than backward Euler for transient problems and is recommended
for Cahn-Hilliard simulations.

**`IterationAdaptiveDT`**: Starts at dt=0.5, grows by 1.5x when Newton converges in
fewer than 8 iterations, and cuts by 0.5x when it needs more. Spinodal kinetics are
fastest early; the time step can grow substantially during late-stage coarsening.

---

## Running the Simulation

```bash
cd quickstart-runs/case18-cahn-hilliard

# Run with the combined application
/path/to/combined-opt -i case18_cahn_hilliard.i
```

The combined application (`modules/combined`) includes the `PhaseField` module, which
provides `SplitCHParsed`, `SplitCHWRes`, `CoupledTimeDerivative`, and
`DerivativeParsedMaterial`. These kernels and materials are not available in the basic
test application.

Estimated wall time on a modern laptop: 2-5 minutes for the full 100-unit run.

Progress can be monitored via the CSV output, which writes `avg_c`, `bulk_energy`,
`min_c`, and `max_c` at every time step.

---

## Expected Results

### Early Time (t < 5)

The small random perturbations in c are amplified. The linear stability analysis of
the Cahn-Hilliard equation predicts that wavelengths in the range 0 < k < k_c will
grow, where k_c = sqrt(|d2F/dc2| / kappa). For this double-well with kappa=1 and
d2F/dc2 = -0.5 at c=0.5, the most unstable wavelength is lambda* ~ 2*pi/k* ~ 9 length
units. You will see a fine-scale microstructure forming with a characteristic spacing
consistent with this prediction.

### Intermediate Time (5 < t < 50)

The domains are well-defined but small. The bulk energy `bulk_energy` has dropped
substantially from its initial value (since more material sits near c=0 and c=1 rather
than at c=0.5). The `min_c` postprocessor approaches 0 and `max_c` approaches 1.
The `avg_c` stays at 0.5 throughout, confirming mass conservation.

### Late Time (t > 50)

Coarsening dominates. Smaller domains shrink and disappear; larger domains grow.
The number of distinct phase regions decreases, and the interface network becomes
coarser. The bulk energy approaches a plateau as the remaining interfaces are few
and slowly moving.

### Postprocessor Trends

| Quantity    | Expected Trend                                        |
|-------------|-------------------------------------------------------|
| `avg_c`     | Constant at 0.5 (mass conservation)                  |
| `bulk_energy` | Monotonically decreasing toward a plateau           |
| `min_c`     | Decreases from ~0.44 toward 0 as phases separate     |
| `max_c`     | Increases from ~0.56 toward 1 as phases separate     |
| `dt_size`   | Grows over time as kinetics slow during coarsening   |

---

## Key Takeaways

| Concept | Where in the Input File |
|---------|------------------------|
| Split Cahn-Hilliard formulation | `[Kernels]` — two equations on `c` and `w` |
| Double-well free energy | `[Materials]` — `DerivativeParsedMaterial` |
| Automatic differentiation of F | `derivative_order = 2` in `free_energy` |
| Periodic boundary conditions | `[BCs]` — `auto_direction = 'x y'` |
| Random initial condition | `[Variables]` — `RandomIC` |
| BDF2 time integration | `[Executioner]` — `scheme = bdf2` |
| Adaptive time stepping | `[Executioner]` — `IterationAdaptiveDT` |
| Mass conservation diagnostic | `[Postprocessors]` — `avg_c` |
| Free energy monitoring | `[Postprocessors]` — `bulk_energy` |
