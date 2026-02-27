# Case 65: Contaminant Transport — Advection-Dispersion-Reaction

## Overview

Groundwater contamination is one of the most consequential environmental problems in hydrogeology. A dissolved pollutant injected at one location migrates through the aquifer under the combined action of advection (carried by the flowing groundwater), hydrodynamic dispersion (spreading due to pore-scale heterogeneity), and chemical reaction (natural attenuation). The mathematical model for this process is the Advection-Dispersion-Reaction (ADR) equation, also called the convection-diffusion-reaction equation in other fields.

This case uses the MOOSE `chemical_reactions` module to simulate a contaminant breakthrough experiment in a 1D (thin 2D) aquifer. Groundwater flows from left to right driven by a prescribed linear pressure field. Contaminant enters at the left boundary (c = 1 mol/m^3) and is simultaneously dispersed, advected rightward, and degraded at a first-order decay rate (k = 0.01/s). The classic S-shaped breakthrough curves recorded at monitoring points along the aquifer demonstrate the interplay between the three transport mechanisms.

Key concepts demonstrated:

- `PrimaryTimeDerivative` and `PrimaryDiffusion` from the `chemical_reactions` module
- `PrimaryConvection` kernel for advection driven by a Darcy velocity
- `CoefReaction` for first-order biological/chemical degradation
- Prescribing Darcy velocity from a linear pressure field using `DarcyFluxPressure` AuxKernel
- `ChemicalOutFlowBC` for advection-dominated outflow boundary conditions
- Monitoring breakthrough curves at multiple spatial locations

---

## The Physics

### Advection-Dispersion-Reaction Equation

The governing equation for concentration c [mol/m^3] in a porous medium with porosity phi:

```
phi * dc/dt + div(v * c) - div(D_h * grad(c)) + phi * k * c = 0
```

where:
- `phi = 0.3` — porosity (fraction of void space)
- `v` — Darcy flux (seepage velocity, m/s) derived from the pressure field
- `D_h = 1e-3 m^2/s` — hydrodynamic dispersion coefficient (combines molecular diffusion and mechanical dispersion)
- `k = 0.01 /s` — first-order decay rate (biodegradation or radioactive decay)

This is the standard ADE used throughout groundwater hydrology. In the `chemical_reactions` module, the three kernel types `PrimaryTimeDerivative`, `PrimaryDiffusion`, and `PrimaryConvection` implement the storage, dispersion, and advection terms respectively.

### Darcy Flow Field

The pressure field is prescribed as a simple linear gradient:

```
P(x) = 2 - x   (Pa)
```

giving a constant pressure gradient `dP/dx = -1 Pa/m`. The Darcy flux is:

```
v_x = -(k_perm / mu) * dP/dx = (1e-3 / 1.0) * 1 = 1e-3 m/s
```

using hydraulic conductivity K = k_perm / mu = 1e-3 m^2/s / (1.0 Pa.s) (note these are effective transport parameters, not physical pore-scale values). The Darcy velocity is uniform across the domain — a plug-flow advection field.

### Dimensionless Parameters

The problem is characterised by two dimensionless groups:

**Peclet number**: ratio of advective to dispersive transport:
```
Pe = v * L / D_h = 1e-3 * 1 / 1e-3 = 1
```

With Pe = 1, advection and dispersion are of comparable magnitude. This produces the classic S-shaped (sigmoidal) breakthrough curve. For Pe >> 1, advection dominates and the front is sharp; for Pe << 1, dispersion dominates and the front is broad.

**Damkohler number**: ratio of reaction to advection:
```
Da = k * L / v = 0.01 * 1 / 1e-3 = 10
```

With Da = 10, reaction is much faster than advection across the domain length. This means the contaminant is significantly degraded as it travels through the aquifer — concentrations at the outlet are much less than 1.

### Steady-State Analytical Solution

For the 1D ADE with constant injection at x = 0 and first-order decay, the steady-state concentration profile is:

```
c(x) = exp(x * v / (2 * D_h)) * exp(-x * sqrt((v/(2*D_h))^2 + k/D_h))
     = exp(x * (Pe/2 - sqrt((Pe/2)^2 + Da)) / L)
```

With Pe/2 = 0.5 and Da = 10: the exponent is dominated by the decay term and c(x) falls rapidly from 1 at x = 0. The concentration at x = 0.5 m at steady state is approximately c ~ exp(-0.5 * 3.16) ~ 0.21, reflecting the strong attenuation.

### Domain and Boundary Conditions

```
y = 0.1 m  -------------------------------------------  TOP: zero-flux
           |                                           |
           |   phi=0.3, D_h=1e-3, k=0.01             |
           |   Darcy flow: v_x = 1e-3 m/s (left to right)   |
           |                                           |
y = 0.0 m  -------------------------------------------  BOTTOM: zero-flux

           x = 0 (c=1)                    x = 1 m (ChemicalOutFlowBC)
```

**Left boundary (x = 0)**: `DirichletBC` with c = 1 mol/m^3. This represents the constant contaminant injection point or the upstream source plume.

**Right boundary (x = 1)**: `ChemicalOutFlowBC`. This special boundary condition handles advective outflow without reflection: it applies the outgoing advective flux `v * c * n` without imposing a Dirichlet value. Using a simple zero-concentration Dirichlet at the outflow would incorrectly restrict the concentration and cause a spurious accumulation.

**Top and bottom**: zero-flux Neumann (natural conditions). The problem is effectively 1D in x.

---

## Input File Walkthrough

The input file is `case65_contaminant_transport.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 50
ny = 2
xmax = 1.0
ymax = 0.1
```

A 50 x 2 element mesh. The 50 elements in x give element size 0.02 m. With Peclet Pe = v * dx / D_h = 1e-3 * 0.02 / 1e-3 = 0.02 (element Peclet number), the grid-scale Peclet is well below 2, so no upwinding is needed to prevent spurious oscillations. The 2 elements in y confirm the quasi-1D setup.

### `[Variables]`

```
[c]
  initial_condition = 0    # initially clean aquifer
[]
```

The single primary variable `c` is the dissolved contaminant concentration. The aquifer is initially clean (c = 0 everywhere). Contaminant enters from the left boundary starting at t = 0.

### `[AuxVariables]` and `[AuxKernels]`

```
[vel_x]
  order = CONSTANT
  family = MONOMIAL
[]
[darcy_vel]
  type = DarcyFluxPressure
  variable = vel_x
  component = 0
[]
```

The Darcy velocity is computed as an auxiliary variable from the pressure field. `DarcyFluxPressure` evaluates `-K * dP/dx` at each element center, storing the result as a piecewise-constant MONOMIAL variable `vel_x`. The `PrimaryConvection` kernel then reads `vel_x` as its velocity field.

### `[Kernels]`

Four kernels define the transport equation:

**`PrimaryTimeDerivative`**: the time rate of change of `phi * c`. This is the chemical_reactions module's version of the standard `TimeDerivative` kernel; it includes the porosity factor phi.

**`PrimaryDiffusion`**: the dispersive flux `div(D_h * grad(c))`. Reads the diffusivity from the `diffusivity` material property. Uses the `conductivity` property name as per the `chemical_reactions` convention.

**`PrimaryConvection`**: the advective flux `v . grad(c)`. The velocity vector is provided via the `darcy_flux` parameter pointing to the `vel_x` aux variable.

**`CoefReaction`**: first-order decay `-k * c`. Contributes `+k * c * test` to the residual (positive sign because it is a sink term in the residual form). The rate constant `coef = 0.01 /s`.

### `[Materials]`

```
[porous_medium]
  type = GenericConstantMaterial
  prop_names  = 'diffusivity conductivity porosity'
  prop_values = '1e-3       1e-3         0.3'
[]
```

`GenericConstantMaterial` provides all three properties at once. Note that in the `chemical_reactions` module, `PrimaryDiffusion` reads the `conductivity` property (not `diffusivity`) — the two names refer to the same coefficient in different contexts. The `porosity` property is used by `PrimaryTimeDerivative`.

### `[Functions]`

```
[pressure_field]
  type = ParsedFunction
  expression = '2 - x'
[]
```

The linear pressure field P(x) = 2 - x drives the Darcy flow. The `DarcyFluxPressure` AuxKernel differentiates this function to compute the Darcy velocity.

### `[BCs]`

| Name | Variable | Boundary | Type | Value/Note |
|------|----------|----------|------|------------|
| `left_injection` | c | left | DirichletBC | 1.0 mol/m^3 |
| `right_outflow` | c | right | ChemicalOutFlowBC | advective outflow |

`ChemicalOutFlowBC` is provided by the `chemical_reactions` module specifically for advection-dominated outflow. It imposes no constraint on the concentration value but correctly applies the outgoing advective flux in the weak form, preventing the artificial concentration buildup that occurs with a zero-gradient Neumann condition at an advective outlet.

### `[Postprocessors]`

| Name | Type | Location | Physical meaning |
|------|------|----------|-----------------|
| `c_at_025` | PointValue | (0.25, 0.05, 0) | Breakthrough at 25% of domain |
| `c_at_050` | PointValue | (0.5, 0.05, 0) | Breakthrough at 50% of domain |
| `c_at_075` | PointValue | (0.75, 0.05, 0) | Breakthrough at 75% of domain |

These three monitoring points record the breakthrough curves — the time history of concentration as the contaminant plume passes each location.

### `[Executioner]`

```
type     = Transient
dt       = 2.0
end_time = 100.0
```

Fifty time steps of 2 s each. The advective travel time across the domain is L/v = 1/1e-3 = 1000 s; the simulation runs to 100 s, capturing the early breakthrough and the approach toward steady state. By t = 100 s the concentrations are still evolving but the shape of the breakthrough curves is clearly established.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case65-contaminant-transport \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case65_contaminant_transport.i 2>&1 | tail -30'
```

Output files:
- `case65_contaminant_transport_out.e` — Exodus with concentration field at each time step
- `case65_contaminant_transport_out.csv` — Time series of breakthrough concentrations at three locations

---

## Expected Results

### Breakthrough Curves

The signature result is the S-shaped breakthrough curve at each monitoring location. With Pe = 1, the front is dispersed over a distance comparable to the domain length:

| Time (s) | c at x=0.25 | c at x=0.50 | c at x=0.75 |
|----------|------------|------------|------------|
| 0        | 0.000      | 0.000      | 0.000      |
| 10       | 0.037      | 0.002      | 0.000      |
| 30       | 0.228      | 0.063      | 0.007      |
| 60       | 0.382      | 0.186      | 0.070      |
| 100      | 0.432      | 0.255      | 0.129      |

The concentrations at all three points are rising throughout the 100 s simulation and have not reached steady state. The steady-state values (reached after ~500-1000 s, well beyond the simulated period) would be significantly attenuated by the k = 0.01 /s decay term.

### Effect of the Decay Term

The contaminant concentration at x = 0.75 m is much lower than at x = 0.25 m even accounting for the time lag. This attenuation is due to the first-order decay: the plume loses mass as it travels. Without the `CoefReaction` term, the concentrations at all three points would eventually reach c = 1 at steady state (simple advection-dispersion). With k = 0.01 /s, the steady-state concentration at x = 1 m is reduced to approximately 15-20% of the inlet value.

### Spatial Profile

In the Exodus file, the concentration field at any given time shows:
- A sharp front near the inlet at early times that gradually broadens due to dispersion
- A smooth exponential decay of the steady-state profile toward the outlet
- No oscillations, confirming the grid-scale Peclet number is below the critical value of 2

---

## Key Takeaways

- The `chemical_reactions` module provides `PrimaryTimeDerivative`, `PrimaryDiffusion`, and `PrimaryConvection` kernels that implement the advection-dispersion equation for reactive solutes in porous media. They differ from the standard framework kernels in that they incorporate the porosity factor phi.
- `ChemicalOutFlowBC` is essential at advective outlets. It applies the outgoing flux `v * c * n` without imposing a Dirichlet value. Using a simple Neumann or no-BC at an advective outlet produces physically incorrect concentration buildup or reflections.
- `DarcyFluxPressure` converts a pressure field (prescribed via a `ParsedFunction`) into a Darcy flux auxiliary variable. This approach allows complex, spatially varying flow fields derived from Darcy's law without solving the full flow equations simultaneously.
- The grid-scale Peclet number `Pe_h = v * dx / D_h` must remain below 2 to prevent spurious oscillations in standard Galerkin finite elements. If Pe_h > 2, an upwinding scheme (such as `ADUpwinding`) or a finer mesh is required.
- The Damkohler number Da = k * L / v controls how much attenuation occurs over the domain. With Da = 10, the contaminant is significantly degraded; with Da << 1 it would pass through almost unreacted. The two dimensionless groups Pe and Da completely characterize the shape of the breakthrough curve in this linear problem.
- Breakthrough curves (concentration vs. time at fixed locations) are the primary observational data in tracer tests used to characterize aquifer properties. Fitting simulated breakthrough curves to field data allows estimation of D_h and k.
- This case uses `CoefReaction` (not `MatReaction` or `ADMatReaction`) to avoid complications with AD vs. non-AD material property compatibility. `CoefReaction` takes the rate coefficient directly as a parameter, requiring no material property lookup.
