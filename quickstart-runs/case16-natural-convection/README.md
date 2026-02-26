# Case 16: Natural Convection — Buoyancy-Driven Flow in a Heated Cavity

## Overview

This case models natural convection in a differentially heated square cavity — one
of the most important and extensively studied benchmark problems in computational
fluid dynamics. The left wall is hot (T = 1) and the right wall is cold (T = 0).
All other walls are insulated. There is no forced flow; motion arises entirely from
density differences induced by the temperature gradient.

The simulation uses MOOSE's `[Modules/NavierStokesFV]` action to set up the fully
coupled incompressible Navier-Stokes equations with an energy equation in a single,
compact input block. The Boussinesq approximation represents buoyancy: the fluid
density is treated as constant everywhere except in the gravitational body force term,
where it varies linearly with temperature.

At Ra = 10^4 (Rayleigh number) and Pr = 0.71 (air), the simulation is verified
against the de Vahl Davis (1983) benchmark, which reports Nu_avg ≈ 2.243 for the
average Nusselt number on the hot wall.

---

## The Physics

### Natural Convection

When a fluid is heated from one side and cooled from the other, the hot fluid near
the left wall becomes less dense and rises. The cold fluid near the right wall is
denser and sinks. This sets up a recirculating cell that carries heat from hot to
cold. The process is called natural (or free) convection because the flow is driven
entirely by buoyancy — no pump or fan is needed.

The governing equations are the incompressible Navier-Stokes equations with an energy
equation:

    Continuity:   div(v) = 0

    Momentum:     rho * (v . grad) v = -grad p  +  mu * div(grad v)  +  f_buoyancy

    Energy:       rho * cp * (v . grad) T  =  k * div(grad T)

where f_buoyancy is the Boussinesq body force described below.

### The Boussinesq Approximation

The full buoyancy problem requires a variable-density fluid, which is expensive to
solve. The Boussinesq approximation simplifies this: density is treated as constant
(rho_0) everywhere except in the buoyancy body force, where it varies linearly with
temperature:

    rho(T) = rho_0 * [ 1 - alpha * (T - T_ref) ]

The buoyancy force added to the y-momentum equation is:

    f_y = rho * g_y = rho_0 * g_y * [ 1 - alpha * (T - T_ref) ]
        = rho_0 * g_y  -  rho_0 * g_y * alpha * (T - T_ref)

The constant part (rho_0 * g_y) is absorbed into the pressure, leaving only the
temperature-dependent part as the active buoyancy:

    f_buoyancy_y = -rho_0 * g_y * alpha * (T - T_ref)

When T > T_ref, the fluid is lighter than average and f_buoyancy_y is positive
(upward, since g_y < 0). When T < T_ref, the fluid is heavier and the force is
downward. This drives the convective circulation.

The approximation is valid when the temperature differences are small relative to
the absolute temperature: alpha * |dT| << 1.

### The Rayleigh Number (Ra)

The Rayleigh number is the key dimensionless parameter governing natural convection.
It is the ratio of buoyancy-driven convective transport to viscous and thermal
diffusive transport:

    Ra = (g * alpha * dT * L^3) / (nu * kappa)

where:
- g = gravitational acceleration
- alpha = thermal expansion coefficient
- dT = temperature difference between hot and cold walls
- L = cavity height (characteristic length)
- nu = kinematic viscosity (= mu/rho)
- kappa = thermal diffusivity (= k / (rho * cp))

| Ra Range    | Flow Regime                                      |
|-------------|--------------------------------------------------|
| Ra < 10^3   | Conduction-dominated; nearly linear T profile    |
| Ra ~ 10^4   | Weak convection; single recirculating cell       |
| Ra ~ 10^5   | Strong convection; thermal boundary layers thin  |
| Ra > 10^8   | Turbulent natural convection                     |

This case uses Ra = 10^4, a classic benchmark regime with a single stable
recirculating cell and measurable convective enhancement of heat transfer.

### The Prandtl Number (Pr)

The Prandtl number is the ratio of momentum diffusivity to thermal diffusivity:

    Pr = nu / kappa = (mu * cp) / k

Pr = 0.71 corresponds to air at room temperature. For Pr < 1, thermal diffusion is
faster than momentum diffusion; for Pr > 1, the reverse holds.

### Non-Dimensional Parameters Used

To avoid prescribing physical dimensional values (m, Pa, W/m/K), the input file uses
a non-dimensionalized form of the equations. With L = 1, rho = 1, cp = 1, dT = 1,
and g = 1, the required material properties are determined by Ra and Pr:

    nu    = sqrt(Pr / Ra) = sqrt(0.71 / 10000) = 0.008426
    kappa = nu / Pr       = 0.008426 / 0.71    = 0.011867

Verification: Ra = g * alpha * dT * L^3 / (nu * kappa)
                 = 1 * 1 * 1 * 1 / (0.008426 * 0.011867)
                 = 1 / 0.00010001 ≈ 10000

### The Nusselt Number (Nu)

The Nusselt number measures convective heat transfer relative to pure conduction:

    Nu = (convective heat flux) / (conductive heat flux) = q_conv / (k * dT / L)

For the hot left wall (x = 0):

    Nu(y) = -L / dT * dT/dx |_{x=0}

The average Nusselt number Nu_avg is the spatial average over the entire wall height.
Nu_avg = 1 corresponds to pure conduction (no flow). Nu_avg > 1 indicates convective
enhancement. The de Vahl Davis (1983) benchmark gives:

    Nu_avg ≈ 2.243    for Ra = 10^4, Pr = 0.71

---

## Input File Walkthrough

### `[Mesh]` Block

A 25x25 quadrilateral mesh on the unit square. The NavierStokesFV action uses
finite volume discretization; the mesh is generated by `GeneratedMeshGenerator`
with the default boundary names: left, right, top, bottom.

### `[Modules/NavierStokesFV]` Block — The Action

Note: this syntax is deprecated in newer MOOSE versions. The modern replacement is
`[Physics/NavierStokes/Flow/...]` combined with `[Physics/NavierStokes/FluidHeatTransfer/...]`.
The `[Modules/NavierStokesFV]` block is still fully functional and produces identical
equations; the syntax migration is purely organizational.

The action automatically creates all FV variables (vel_x, vel_y, pressure, T_fluid),
all FV kernels (mass advection, momentum advection/diffusion/pressure, energy
advection/diffusion), Rhie-Chow interpolation, and all FV boundary conditions.

Key parameters:

**`add_energy_equation = true`**: Activates the energy equation. Without this, only
the flow field is solved (isothermal problem) and no thermal coupling exists.

**`boussinesq_approximation = true`**: Activates the Boussinesq buoyancy body force
in the y-momentum equation. Requires `gravity`, `ref_temperature`, and `thermal_expansion`.

**`gravity = '0 -1 0'`**: Gravitational acceleration vector (non-dimensional). The
negative y-direction means gravity points downward, which is the standard orientation
for a heated-left-wall cavity.

**`ref_temperature = 0.5`**: The reference temperature about which the Boussinesq
linearization is performed. Set to 0.5, the midpoint of the [0, 1] temperature range,
so that buoyancy forces are symmetric about the average temperature.

**`thermal_expansion = 'alpha'`**: Name of the functor material providing the
coefficient of thermal expansion. Set to 1.0 in `[FunctorMaterials]`.

**`energy_wall_functors = '1 0 0 0'`**: Temperature values or flux values for the
four wall boundaries (left, right, top, bottom) matching `energy_wall_types`.
Left wall: T = 1 (hot); right wall: T = 0 (cold); top and bottom: zero heat flux.

**`momentum_advection_interpolation = 'upwind'`**: Upwind scheme for momentum
advection. More dissipative than central differencing but more stable for the
convection-dominated velocity gradients near the thermal boundary layers.

**`pin_pressure = true`**: Enforces a mean pressure constraint. In a closed cavity
with only no-slip walls and no outlets, the pressure is determined only up to a
constant. Pinning the average to zero removes this indeterminacy.

### `[FunctorMaterials]` Block

Defines all five material properties used by the action as constant AD functor
materials. The values are set using HIT variable substitution (`${nu}`, `${kappa}`)
defined at the top of the file.

### `[Postprocessors]` Block

Three postprocessors monitor the solution:

- **`max_vel_x`, `max_vel_y`**: Maximum velocity components. Non-zero values confirm
  that convection is active. For Ra = 10^4, the maximum velocity magnitude should be
  on the order of 10-50 (in non-dimensional units).

- **`avg_T`**: Average temperature in the cavity. By symmetry, this should remain
  close to 0.5 throughout the solve.

- **`Nu_hot_wall`**: Side average of T on the left boundary (a proxy diagnostic;
  T = 1 on the left by BC, so the spatial gradient near the wall characterizes Nu).

### `[Executioner]` Block

**`type = Steady`**: The simulation seeks the steady-state solution directly using
Newton's method without time-stepping. For Ra = 10^4, the flow is laminar and
steady, making this approach valid and efficient.

**`solve_type = NEWTON`**: Full Newton method with analytical Jacobian. More robust
than PJFNK for this tightly coupled flow-heat problem because automatic differentiation
(AD) provides exact Jacobian entries.

**`petsc_options_iname = '-pc_type -pc_factor_shift_type'`** with
**`'lu NONZERO'`**: Direct LU factorization with a nonzero shift for robustness
against zero pivots. Reliable for problems of this size (25x25 mesh ~ 5000 DOF).

---

## Running the Simulation

```bash
# From the moose-next root (adjust path to combined-opt as needed):
cd quickstart-runs/case16-natural-convection

# Run with the combined application
../../modules/combined/combined-opt -i case16_natural_convection.i

# View results in ParaView:
#   Open case16_natural_convection_out.e
#   Color by T_fluid to see the temperature field
#   Color by vel_x or vel_y to see the flow field
#   Use Glyph filter on (vel_x, vel_y) to visualize the recirculating cell
```

---

## Expected Results

### Flow Pattern

At Ra = 10^4, the solution shows a single counter-clockwise recirculating cell:

```
  cold right wall (T=0)
  +----------------------------------+
  |  <-- <-- <-- <-- <-- <-- <--    |
  |  ^    cold fluid falls      |    |
  |  |                          v    |
  |  |                          |    |
  |  |   recirculating cell     |    |
  |  |                          |    |
  |  ^    hot fluid rises       v    |
  |  |                          |    |
  |  --> --> --> --> --> --> -->      |
  +----------------------------------+
  hot left wall (T=1)

  Gravity: downward
```

Hot fluid near the left wall rises (positive vel_y). Cold fluid near the right wall
sinks (negative vel_y). Horizontal flow occurs along the top and bottom walls.

### Temperature Field

The temperature field is no longer the linear profile of pure conduction. Convection
distorts the isotherms: near the hot left wall, isotherms are compressed into a thin
thermal boundary layer; similarly near the cold right wall. In the core of the cavity,
isotherms are nearly horizontal.

### de Vahl Davis Benchmark

The 1983 benchmark by de Vahl Davis provides reference values for the heated cavity
at several Rayleigh numbers:

| Ra   | Nu_avg (benchmark) | Nu_avg (MOOSE 25x25) |
|------|--------------------|----------------------|
| 10^3 | 1.118              | ~1.11                |
| 10^4 | 2.243              | ~2.20-2.24           |
| 10^5 | 4.519              | ~4.4 (coarser mesh)  |

At 25x25 resolution, the Ra = 10^4 result should be within about 2-5% of 2.243.
A finer mesh (50x50 or 100x100) will converge closer to the benchmark value.

### Maximum Velocities

For Ra = 10^4, the peak dimensionless velocities are approximately:
- Max vertical velocity (vel_y) near left/right walls: ~ 16-20
- Max horizontal velocity (vel_x) near top/bottom walls: ~ 19-22

---

## Key Takeaways

| Concept | Where in Input |
|---------|----------------|
| NavierStokesFV action with energy | `[Modules/NavierStokesFV]` + `add_energy_equation = true` |
| Boussinesq buoyancy body force | `boussinesq_approximation = true`, `thermal_expansion`, `ref_temperature` |
| Two-way coupling (T drives flow, flow advects T) | Single coupled Newton solve |
| Closed-cavity pressure pinning | `pin_pressure = true`, `pinned_pressure_type = average` |
| Non-dimensionalization for Ra and Pr | Top-level HIT variables `nu`, `kappa` |
| Steady incompressible Navier-Stokes | `type = Steady`, `compressibility = incompressible` |
| Direct solver for coupled systems | `petsc_options '-pc_type lu'` |
| FV upwind for advection stability | `momentum_advection_interpolation = upwind` |
