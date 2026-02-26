# Case 37: Rayleigh-Benard Convection Onset

## Overview

Rayleigh-Benard convection is the canonical instability of a fluid layer heated
from below and cooled from above. When the bottom plate is hot and the top plate is
cold, the fluid near the bottom becomes buoyant and wants to rise while the denser
cold fluid above wants to sink. For small temperature differences, viscosity and
thermal diffusion suppress any motion and heat is carried by conduction alone — the
temperature profile is a linear gradient from hot to cold. Above a sharply defined
critical Rayleigh number, Ra_c = 1708, this conductive state becomes linearly
unstable and the fluid spontaneously organizes into a pattern of convective rolls.

The transition at Ra_c = 1708 is one of the most precisely determined instability
thresholds in all of fluid mechanics and was established analytically by Rayleigh
(1916) for free-slip plates and confirmed numerically and experimentally for no-slip
(rigid) plates. This case simulates Ra = 2000, just 17% above the critical value,
so the convective rolls are weak but clearly established. The Nusselt number Nu > 1
at steady state confirms that convective heat transport augments pure conduction.
The aspect ratio 2:1 domain (width = 2, height = 1) accommodates exactly one pair
of counter-rotating rolls at the most unstable horizontal wavenumber k_c ≈ 3.117.

This case deliberately contrasts with Case 16 (natural convection in a differentially
heated cavity with a hot left wall, Ra = 10000, steady solver). Case 37 uses a hot
bottom plate, operates near the onset threshold, and runs as a transient so the roll
formation process is captured in time.

---

## The Physics

### Governing Equations

The Boussinesq incompressible Navier-Stokes equations coupled to an energy equation:

```
Continuity:   div(u) = 0

Momentum:     rho * (du/dt + u.grad(u)) = -grad(p) + mu * lap(u) + f_buoy

Energy:       dT/dt + u.grad(T) = kappa * lap(T)
```

where the Boussinesq buoyancy body force is:

```
f_buoy = -rho * alpha * (T - T_ref) * g
```

With gravity pointing in the -y direction (g = (0, -1, 0)), a fluid parcel hotter
than T_ref experiences an upward (positive y) buoyancy force; a colder parcel is
pushed downward. This is the mechanism that drives convective overturning once the
temperature difference is large enough to overcome viscous and diffusive damping.

### The Rayleigh Number and Critical Threshold

The Rayleigh number Ra characterises the ratio of buoyancy-driven convective transport
to the product of viscous and thermal diffusive transport:

```
Ra = g * alpha * dT * H^3 / (nu * kappa)
```

where H is the layer height (characteristic length), dT is the temperature difference
between bottom and top plates, nu is the kinematic viscosity, and kappa is the thermal
diffusivity. When Ra is small, diffusion wins and the fluid remains quiescent. Above
Ra_c = 1708 (rigid no-slip plates, first established by Jeffreys 1928 and confirmed
by Pellew & Southwell 1940), the conduction solution is linearly unstable to infinitesimal
perturbations. Rolls grow exponentially until nonlinear saturation sets the amplitude.

The non-dimensional parameters for this case:

| Symbol | Value  | Formula                          |
|--------|--------|----------------------------------|
| Ra     | 2000   | target (above Ra_c = 1708)       |
| Pr     | 0.71   | air at room temperature          |
| nu     | 0.01884 | sqrt(Pr/Ra) = sqrt(0.71/2000)   |
| kappa  | 0.02653 | nu/Pr = 0.01884/0.71            |

Verification: Ra = 1/(nu * kappa) = 1/(0.01884 * 0.02653) ≈ 2000.

### The Prandtl Number

The Prandtl number Pr = nu/kappa is the ratio of momentum diffusivity to thermal
diffusivity. For air Pr = 0.71, meaning thermal diffusion is slightly faster than
viscous diffusion. For water Pr ≈ 7; for liquid metals Pr << 1.

### The Most Unstable Mode

Linear stability analysis of the Rayleigh-Benard problem shows that the most
dangerous perturbation has a horizontal wavenumber k_c = pi/sqrt(2) ≈ 2.221 and
a critical Ra_c = 27*pi^4/4 ≈ 1708 (for no-slip plates the exact value is 1707.76).
The corresponding roll wavelength is lambda_c = 2*pi/k_c ≈ 2*H. For the unit-height
domain used here, lambda_c ≈ 2, so the 2:1 domain (width = 2) accommodates exactly
one wavelength, which is why a single pair of rolls forms rather than a disordered
pattern.

The sinusoidal temperature perturbation in the initial condition,
T = 1 - y + 0.01*sin(pi*x)*sin(pi*y), projects directly onto this most unstable
mode and triggers organised roll formation without waiting for numerical noise to
amplify.

### The Nusselt Number

The Nusselt number Nu measures convective heat transport enhancement over pure
conduction:

```
Nu = (total heat flux) / (conductive heat flux) = q_total / (k * dT / H)
```

For the conduction solution Nu = 1 exactly. Above onset, Nu grows as:

```
Nu - 1 ≈ C * (Ra - Ra_c) / Ra_c    (weakly nonlinear theory near onset)
```

so at Ra = 2000, with Ra_c = 1708, we expect Nu - 1 ≈ 0.17 * C, giving Nu
modestly above 1. The postprocessors T_hot_wall and T_cold_wall provide a proxy:
if heat is being redistributed by convection, the surface-averaged temperatures
deviate from the pure-conduction values (1 and 0 respectively) due to the flow
sweeping hot and cold fluid to different regions of each plate.

---

## Input File Walkthrough

### Top-level Variables

```
nu    = 0.01884   # kinematic viscosity = dynamic viscosity (rho=1)
kappa = 0.02653   # thermal diffusivity = conductivity (rho*cp=1)
```

These are computed from Ra = 2000 and Pr = 0.71 using nu = sqrt(Pr/Ra) and
kappa = nu/Pr. Setting rho = 1 and cp = 1 means nu = mu and kappa = k, so the
same numerical values serve double duty as both kinematic and dynamic quantities.

### `[Mesh]` Block

A 40x20 quadrilateral mesh on the 2x1 rectangle. The 2:1 aspect ratio matches the
critical roll wavelength so a single counter-rotating pair fits naturally. GeneratedMesh
assigns boundary names: bottom (y = 0, hot plate), top (y = 1, cold plate), left and
right (x = 0 and x = 2, side walls).

### `[Modules/NavierStokesFV]` Block

The NavierStokesFV action sets up the complete coupled incompressible flow plus energy
system from a single compact block. Key parameters for this case:

**`add_energy_equation = true`**: Activates the energy (temperature) equation. Without
this the flow is isothermal and no buoyancy is possible.

**`boussinesq_approximation = true`**: Adds the Boussinesq body force
-rho*alpha*(T - T_ref)*g to the momentum equation. Requires `gravity`,
`ref_temperature`, and `thermal_expansion`.

**`gravity = '0 -1 0'`**: Gravity points in the -y direction. With the hot plate at
the bottom (y = 0) and cold at the top (y = 1), this is the destabilising orientation
that produces Rayleigh-Benard convection. If gravity were reversed, the hot-bottom
configuration would be stabilising.

**`ref_temperature = 0.5`**: The Boussinesq linearisation reference. Set to the mean
temperature (midpoint of [0,1]) so that buoyancy forces are antisymmetric about the
centre of the layer.

**`wall_boundaries = 'left right top bottom'`** with
**`momentum_wall_types = 'noslip noslip noslip noslip'`**: All four walls are rigid
no-slip boundaries. This is the classic Rayleigh-Benard setup with rigid plates and
insulated side walls, yielding the threshold Ra_c = 1708.

**`energy_wall_types = 'heatflux heatflux fixed-temperature fixed-temperature'`**
**`energy_wall_functors = '0 0 0 1'`**: The functor order matches wall_boundaries
(left, right, top, bottom):
- left: zero heat flux (insulated)
- right: zero heat flux (insulated)
- top: fixed T = 0 (cold plate)
- bottom: fixed T = 1 (hot plate)

Note that fixed-temperature comes before heatflux in the type string but the last
two entries (top = 0, bottom = 1) are the temperatures. The ordering strictly follows
the `wall_boundaries` list.

**`pin_pressure = true`**: In a closed cavity with only no-slip walls and no pressure
inlets/outlets, the absolute pressure level is undefined (only gradients matter).
Pinning the average pressure to zero removes this null space and allows the linear
solver to converge.

### `[FunctorMaterials]` Block

Five constant AD functor material properties used by the action: rho = 1.0, mu = nu,
k = kappa, cp = 1.0, alpha = 1.0. With rho = cp = 1 and alpha = 1 the Boussinesq
buoyancy magnitude is fully controlled by the gravity vector and temperature difference.

### `[Functions]` and `[FVICs]` Blocks

The temperature initial condition is:

```
T(x, y, 0) = 1 - y + 0.01 * sin(pi*x) * sin(pi*y)
```

The term (1 - y) is the linear conduction profile. The small sinusoidal perturbation
of amplitude 0.01 breaks the symmetry and seeds the most unstable mode. Without this
perturbation the conduction state is an exact solution even above Ra_c and the solver
might not depart from it (or would take much longer to grow from numerical noise).

`FVFunctionIC` applies a `ParsedFunction` to a finite-volume variable (T_fluid).

### `[Postprocessors]` Block

Five diagnostics track the solution:

- **`max_vel_x`, `max_vel_y`**: Maximum velocity components. In the conduction state
  both are zero (or the initial 1e-15). Once rolls form, these grow to measurable
  values. For Ra just above onset the velocities are small but clearly nonzero.

- **`avg_T`**: Average temperature in the domain. Should remain near 0.5 throughout
  by the antisymmetry of the heating configuration.

- **`T_hot_wall`**: Surface-averaged T on the bottom boundary. The BC fixes T = 1
  pointwise, so this postprocessor confirms the BC is being applied and serves as
  a reference.

- **`T_cold_wall`**: Surface-averaged T on the top boundary. Similarly fixed at T = 0.

A proper Nusselt number requires integrating the normal temperature gradient over the
hot wall. These postprocessors provide the simpler proxy; a more precise Nu calculation
would use `FVPostprocessors` or side-integrated flux objects from the navier_stokes
module.

### `[Executioner]` Block

**`type = Transient`**: Unlike Case 16 (Steady), this case resolves the time evolution
so the roll-formation transient is visible. The conduction state at t = 0 transitions
to the convective roll pattern over a few thermal diffusion times.

**`dt = 0.5`** with **`IterationAdaptiveDT`**: The adaptive time-stepper starts at
dt = 0.5 and adjusts based on Newton iteration count. During the fast growth phase
(exponential amplification of the instability) the nonlinear solver may require more
iterations and the time step is cut back; once the rolls reach steady amplitude, fewer
iterations are needed and the step grows.

**`end_time = 20.0`**: Approximately 20 thermal diffusion times (kappa * t / H^2 with
kappa ≈ 0.027 and H = 1 gives diffusion time ~ 37 time units). At end_time = 20,
the rolls should be well into their nonlinear saturation phase and approaching steady
state.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case37-rayleigh-benard \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case37_rayleigh_benard.i 2>&1 | tail -30'
```

Output files produced:
- `case37_rayleigh_benard_out.e` — Exodus file with vel_x, vel_y, pressure, T_fluid
  at each output timestep
- `case37_rayleigh_benard_out.csv` — postprocessor time series (max_vel_x,
  max_vel_y, avg_T, T_hot_wall, T_cold_wall)

To visualise in ParaView:
- Open the Exodus file and animate
- Color by T_fluid to watch the conduction profile deform into roll-distorted isotherms
- Color by vel_y to see the ascending (positive) and descending (negative) plumes
- Apply Glyph filter to (vel_x, vel_y) to display the roll circulation vectors

---

## Expected Results

### Conduction State (t = 0)

At t = 0 the temperature is the linear profile T = 1 - y plus a tiny perturbation.
Velocities are effectively zero. The postprocessor CSV shows max_vel_x and max_vel_y
at 1e-15 (the initial condition floor).

### Roll Formation (t = 1 to 10)

The perturbation grows exponentially. Because Ra = 2000 is only 17% above Ra_c = 1708,
the growth rate sigma = kappa/H^2 * (Ra/Ra_c - 1) * pi^2 ≈ 0.027 * 0.17 * 10 ≈ 0.046
per unit time, giving an e-folding time of about 22 time units. With a 1% initial
perturbation, the rolls reach order-unity amplitude around t = 5 to 10. The adaptive
time-stepper may cut back dt during this growth phase.

### Convective Steady State (t > 15)

The rolls saturate through nonlinear effects and reach a quasi-steady pattern:

```
   cold plate (T = 0)
   +--------------------+
   |  down   |   up     |
   |    <----|---->     |
   |  cold   |  hot     |
   |    <----|---->     |
   |  down   |   up     |
   +--------------------+
   hot plate (T = 1)

   Gravity: downward
```

One pair of counter-rotating rolls fills the 2x1 domain. The left roll descends
(negative vel_y, cold fluid sinking) and the right roll ascends (positive vel_y,
hot fluid rising), or vice versa depending on which symmetry the small perturbation
selects.

### Quantitative Indicators

| Quantity       | Conduction state | Convective state |
|----------------|-----------------|-----------------|
| max_vel_y      | ~0              | 0.1 to 0.5      |
| avg_T          | ~0.5            | ~0.5 (symmetric)|
| Nu (inferred)  | 1.0             | > 1.0           |

The exact convective velocity depends on how far Ra is above Ra_c. Weakly nonlinear
theory (Malkus & Veronis 1958) predicts the roll amplitude scales as sqrt(Ra - Ra_c),
so at Ra = 2000 the rolls are about sqrt(292/1708) ≈ 0.41 of the saturation amplitude
reached at, say, Ra = 3416.

### Comparison with Case 16

| Feature          | Case 16                   | Case 37                        |
|------------------|---------------------------|--------------------------------|
| Heating geometry | Hot left wall             | Hot bottom plate               |
| Ra               | 10000                     | 2000 (just above Ra_c = 1708)  |
| Executioner      | Steady                    | Transient                      |
| Roll topology    | Single side-wall cell     | Horizontal convective rolls    |
| Physical context | Cavity natural convection | Rayleigh-Benard onset          |

---

## Key Takeaways

- The critical Rayleigh number Ra_c = 1708 for rigid no-slip plates is a precise
  bifurcation threshold: below it, the conductive state is the unique stable
  solution; above it, convective rolls are the attractor.

- The Boussinesq approximation reduces the variable-density problem to a constant-
  density flow with a temperature-dependent body force. It is valid when
  alpha * dT << 1 — satisfied here since alpha = 1 and dT = 1, which is marginal;
  the non-dimensional formulation ensures internal consistency.

- Seeding the most unstable mode (sin(pi*x)*sin(pi*y) perturbation) dramatically
  accelerates roll formation compared to waiting for numerical round-off to grow.
  This is essential for transient simulations near onset where growth rates are slow.

- The 2:1 aspect ratio is chosen to match the critical wavelength lambda_c = 2H.
  A 4:1 domain would admit two pairs of rolls; a 1:1 domain would suppress the most
  unstable mode and shift the effective Ra_c upward.

- Adaptive time-stepping (IterationAdaptiveDT) automatically accelerates during
  the slow initial phase and the quasi-steady phase, while refining through the
  rapid roll-formation transient. This is more efficient than a fixed small dt.

- The `pin_pressure = true` option is essential for closed-cavity problems. Without
  it the linear system for pressure is singular (the pressure is determined only up
  to an additive constant in an incompressible closed domain).

- Running at Ra = 2000 rather than Ra = 10000 (Case 16) keeps the convective
  velocities small, making the weakly nonlinear analytical predictions quantitatively
  testable and the nonlinear solver convergence more reliable.

---

## Reference

Rieutord, M., *Fluid Dynamics: An Introduction* (Springer, 2015), Chapter 7,
Section 7.5: Rayleigh-Benard Convection.

Chandrasekhar, S., *Hydrodynamic and Hydromagnetic Stability* (Oxford, 1961),
Chapter 2: The Stability of a Layer of Fluid Heated Below.

Malkus, W. V. R. & Veronis, G. (1958). Finite amplitude cellular convection.
*Journal of Fluid Mechanics*, 4, 225-260.
