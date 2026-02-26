# Case 42: Sod Shock Tube — 1D Riemann Problem

## Overview

The Sod shock tube is the canonical test problem for compressible flow solvers.
A thin membrane at x = 0.5 separates a high-pressure gas on the left from a
low-pressure gas on the right. When the membrane ruptures at t = 0, three
distinct waves emerge from the diaphragm location simultaneously:

1. A **leftward rarefaction fan** — a smooth, isentropic expansion that
   accelerates the left-side gas toward the interface.
2. A **rightward contact discontinuity** — a surface separating gas that
   originated from the two sides; density and velocity jump across it but
   pressure and velocity are continuous.
3. A **rightward shock wave** — a sharp discontinuity that compresses and
   accelerates the undisturbed right-side gas, governed by the
   Rankine-Hugoniot jump conditions.

This problem was introduced by Gary Sod in 1978 as a standard benchmark for
numerical methods applied to the Euler equations. An exact analytical solution
exists, making it ideal for validating both wave propagation speed and
post-wave state values.

Reference: Rieutord, *Fluid Dynamics* (Springer, 2015), Ch. 5, Sec. 5.5.

---

## Governing Equations

The compressible Euler equations in conservative form govern an inviscid ideal gas:

```
d/dt [rho    ]   d/dx [rho*u              ]   [0]
d/dt [rho*u  ] + d/dx [rho*u^2 + p        ] = [0]
d/dt [rho*E  ]   d/dx [(rho*E + p)*u      ]   [0]
```

where:
- rho = density
- u = velocity
- p = pressure
- E = specific total energy (internal + kinetic)

For an ideal gas with ratio of specific heats gamma = 1.4:

```
p = (gamma - 1) * rho * e       (equation of state)
E = e + u^2/2                   (total specific energy)
e = p / ((gamma - 1) * rho)     (specific internal energy)
```

The conservative variables evolved by MOOSE are rho, rho*u, and rho*E.

### Rankine-Hugoniot Jump Conditions

At the shock wave, mass, momentum, and energy fluxes must be continuous
across the discontinuity. If the shock moves at speed s, then:

```
[rho*(u - s)]  = 0
[rho*u*(u - s) + p]  = 0
[(rho*E + p)*(u - s)] = 0
```

These conditions uniquely determine the post-shock state given the pre-shock
state and the shock Mach number.

---

## Initial and Boundary Conditions

### Initial Conditions (Sod problem, standard values)

| Region | rho   | p   | u | rho*E |
|--------|-------|-----|---|-------|
| Left (x < 0.5)  | 1.000 | 1.0 | 0 | 2.500 |
| Right (x > 0.5) | 0.125 | 0.1 | 0 | 0.250 |

The initial total energy per unit volume is:

```
rho*E = p / (gamma - 1)    (at rest, kinetic energy = 0)

Left:   rho*E = 1.0 / 0.4 = 2.5
Right:  rho*E = 0.1 / 0.4 = 0.25    (= rho * E_specific = 0.125 * 2.0)
```

### Boundary Conditions

Reflective (implicit) boundary conditions at x = 0 and x = 1 using
`CNSFVHLLCMassImplicitBC`, `CNSFVHLLCMomentumImplicitBC`, and
`CNSFVHLLCFluidEnergyImplicitBC`. These enforce zero normal velocity at the
domain walls while allowing the interior waves to develop naturally. The domain
is wide enough that neither the rarefaction nor the shock reaches the walls
before t = 0.2.

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 1D, 200 cells, [0, 1] | Uniform finite volume mesh |
| `IdealGasFluidProperties` (gamma=1.4) | Equation of state for air |
| `CNSFVMassHLLC` | HLLC Riemann flux for continuity equation |
| `CNSFVMomentumHLLC` | HLLC Riemann flux for momentum equation |
| `CNSFVFluidEnergyHLLC` | HLLC Riemann flux for energy equation |
| `CNSFVHLLCMassImplicitBC` | Reflective wall BC for density |
| `CNSFVHLLCMomentumImplicitBC` | Reflective wall BC for momentum |
| `CNSFVHLLCFluidEnergyImplicitBC` | Reflective wall BC for energy |
| `ConservedVarValuesMaterial` | Reconstructs pressure/velocity from conservative vars |
| `FunctionIC` piecewise constant | Left/right initial states separated at x=0.5 |
| `ExplicitSSPRungeKutta` order=2 | Strong Stability Preserving RK2 time integrator |
| `ElementExtremeValue` max_rho, min_rho | Tracks density range (rarefaction vs shock) |
| `ElementIntegralVariablePostprocessor` | Mass and energy conservation checks |

### HLLC Flux Scheme

The Harten-Lax-van Leer-Contact (HLLC) approximate Riemann solver provides
upwinding that respects the three-wave structure of the Euler equations. It
introduces a contact wave estimate alongside the left and right acoustic wave
speeds, recovering the exact solution for isolated contact discontinuities.
This makes it significantly more accurate than the original HLL scheme for
problems where the contact wave carries a density jump — exactly the situation
in the Sod problem.

### Time Step and CFL Condition

With 200 cells on [0, 1] the cell width is dx = 0.005. The maximum wave speed
at t = 0 is c_left = sqrt(gamma * p / rho) = sqrt(1.4) ≈ 1.183. The CFL
condition requires dt < dx / c_max, giving dt < 0.005 / 1.183 ≈ 0.0042.
The chosen dt = 5e-4 provides a CFL number of approximately 0.12, well within
the stability region of SSP-RK2.

---

## Expected Results at t = 0.2

| Feature | Location | Value |
|---------|----------|-------|
| Left domain (undisturbed) | x < 0.26 | rho = 1.0, p = 1.0 |
| Rarefaction fan (left edge) | x ≈ 0.26 | Begins smooth expansion |
| Rarefaction fan (right edge) | x ≈ 0.49 | Ends at contact state |
| Contact discontinuity | x ≈ 0.69 | rho jumps: 0.426 left, 0.265 right |
| Post-shock state | 0.69 < x < 0.85 | rho ≈ 0.265, p ≈ 0.303 |
| Shock wave front | x ≈ 0.85 | Sharp density/pressure jump |
| Right domain (undisturbed) | x > 0.85 | rho = 0.125, p = 0.1 |

### Conserved Quantities

- `total_mass` should remain constant at 0.125 * 1.0 + 0.125 * 0.5 * 2 = ...
  Numerically: integral of rho over [0, 1] = initial value ≈ 0.5625
- `total_rho_E` should remain constant at the initial value ≈ 0.5 * 2.5 + 0.5 * 0.25 = 1.375

### Postprocessor Values at t = 0.2

| Postprocessor | Expected value |
|---------------|---------------|
| max_rho | ≈ 0.426 (post-contact density) |
| min_rho | ≈ 0.125 (undisturbed right state) |
| max_rho_u | > 0 (gas moving right through contact/shock) |
| total_mass | ≈ 0.5625 (conserved) |
| total_rho_E | ≈ 1.375 (conserved) |

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case42-sod-shock-tube \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case42_sod_shock_tube.i 2>&1 | tail -20'
```

Output files produced:
- `case42_sod_shock_tube_out.e` — Exodus file with rho, rho_u, rho_E at every timestep
- `case42_sod_shock_tube_out.csv` — Postprocessor time series (mass, energy, extremes)

The run takes approximately 400 time steps (end_time=0.2, dt=5e-4) and
completes in under a minute with combined-opt.

---

## Key Takeaways

1. **Three-wave structure**: The Euler equations have three characteristic wave
   families. Any Riemann problem for a perfect gas produces exactly these three
   waves (rarefaction, contact, shock), whose speeds and post-wave states are
   determined by the initial jump and the equation of state.

2. **Conservation is exact**: Mass and total energy are globally conserved by
   the finite volume scheme to machine precision; only boundary flux terms can
   change them, and the reflective walls contribute zero net flux here.

3. **Contact vs shock**: The contact discontinuity is a linearly degenerate
   wave — it propagates without change and carries no pressure jump. The shock
   is a genuinely nonlinear wave that irreversibly converts kinetic energy to
   internal energy via compression, increasing entropy.

4. **HLLC accuracy**: The HLLC solver exactly captures the contact wave speed,
   producing a sharper contact discontinuity than HLL while remaining stable
   across the shock. With 200 cells the shock is captured in about 3-4 cells.

5. **Godunov-type methods**: This case illustrates the core idea of modern
   compressible flow numerics: solve a local Riemann problem at each cell face
   to compute the inter-cell flux. HLLC is one of the most widely used
   approximate Riemann solvers in production CFD codes.
