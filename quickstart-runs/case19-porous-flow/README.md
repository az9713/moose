# Case 19: Darcy Flow with Heat Transport in Porous Media

## Overview

This case couples fluid flow and heat transport inside a porous rock matrix —
the governing physics of geothermal reservoirs, underground CO2 storage, and
petroleum extraction. A pressure gradient drives single-phase liquid through a
2-D rectangular domain (3 m x 2 m), and hot fluid injected from the left boundary
creates a thermal plume that slowly advects rightward. Watching the plume migrate
is a fundamental verification of coupled thermo-hydro (TH) porous-media models.

The simulation uses MOOSE's `PorousFlow` module, specifically the
`PorousFlowBasicTHM` Action, which assembles the entire coupled system (mass
balance, energy balance, material hierarchy, Darcy flow) from a compact input
block rather than requiring every kernel and material to be listed individually.

Key features demonstrated:
- `PorousFlowBasicTHM` Action for one-line TH setup
- `SimpleFluidProperties` equation of state for water-like fluid
- Pressure-driven Darcy flow with pressure inlet and outlet BCs
- Hot-fluid injection BC with thermal retardation by rock matrix
- Transient simulation over 5000 s observing thermal plume migration

---

## The Physics

### Darcy's Law

Flow in a porous medium obeys Darcy's law:

```
q = -(k / mu) * grad(P)
```

where:
- `q`   -- Darcy flux (volumetric flow per unit area) [m/s]
- `k`   -- permeability tensor [m^2]  (here k = 1e-12 m^2, ~1 millidarcy)
- `mu`  -- dynamic viscosity [Pa.s]   (here 0.001 Pa.s, liquid water at 20 C)
- `P`   -- pore pressure [Pa]

The mass balance (continuity equation) for the fluid is:

```
d(phi * rho_f) / dt  +  div(rho_f * q)  =  0
```

where `phi = 0.3` is porosity. For an incompressible fluid this reduces to
`div(q) = 0` — the Laplace equation in pressure with Darcy velocity as the flux.

### Heat Transport in Porous Media

The energy balance for the coupled fluid-rock system is:

```
d(E_total) / dt  +  div(rho_f * h_f * q)  -  div(lambda * grad(T))  =  0
```

where:
- `E_total = phi * rho_f * u_f + (1-phi) * rho_s * c_s * T`  -- total volumetric
  energy (fluid internal energy + rock sensible heat)
- `rho_f * h_f * q` -- advective energy flux (enthalpy carried by Darcy flow)
- `lambda`  -- effective thermal conductivity of the saturated medium [W/(m K)]
- `h_f`    -- specific enthalpy of the fluid [J/kg]

The first term is heat storage in both fluid and rock. The second term is advection
of heat by the Darcy velocity. The third term is conductive heat diffusion.

### Thermal Retardation

The thermal front moves slower than the fluid itself. The retardation factor R is:

```
R = [ phi * rho_f * c_f + (1 - phi) * rho_s * c_s ] / (phi * rho_f * c_f)
```

For this case:
- phi = 0.3,  rho_f = 1000 kg/m^3,  c_f = 4186 J/(kg K)
- rho_s = 2600 kg/m^3,  c_s = 800 J/(kg K)

```
R = [ 0.3 * 1000 * 4186 + 0.7 * 2600 * 800 ] / (0.3 * 1000 * 4186)
  = [ 1255800 + 1456000 ] / 1255800
  = 2711800 / 1255800  ~  2.16
```

The thermal front travels at approximately 1/R = 46% of the fluid Darcy velocity.

### Domain and Boundary Conditions

```
y=2   (adiabatic, no-flux Neumann -- automatic)
      +-----------------------------------------------+
      |                                               |
P=1.1 |  --> Darcy flow -->  T plume advects right   | P=1.0
T=350 |  (pressure gradient drives fluid rightward)  | T=300
(hot) |                                               | (cool)
      |  Initial state: P=1.0 MPa, T=300 K           |
      |  Mesh: 30 x 20 (600 elements)                |
      +-----------------------------------------------+
y=0   (adiabatic, no-flux Neumann -- automatic)
     x=0                                            x=3
```

Darcy velocity (estimated):
```
q  =  k/mu * dP/L  =  1e-12 / 0.001 * 0.1e6 / 3  ~  3.3e-8 m/s
```

At this Darcy velocity, fluid traversal time is ~91000 s. The thermal front,
retarded by the rock, traverses in ~91000 * R ~ 200000 s.  Over the 5000 s
simulation, the front advances roughly `3.3e-8 / R * 5000 ~ 0.076 m` — a
visible but partial advance, clearly showing the thermal plume beginning to
migrate without yet reaching the outlet.

---

## Input File Walkthrough

### `[GlobalParams]` Block

```
[GlobalParams]
  PorousFlowDictator = dictator
[]
```

The `PorousFlowDictator` UserObject is the central bookkeeper for the
`PorousFlow` module. It tracks the number of fluid phases, fluid components,
and PDE variables so that kernels and materials can coordinate. The name
`dictator` is a convention; all PorousFlow objects look for it via this
`GlobalParams` entry so it does not need to be declared in each object.

### `[Variables]` Block

```
[Variables]
  [porepressure]
    initial_condition = 1e6
  []
  [temperature]
    initial_condition = 300
    scaling = 1e-6
  []
[]
```

Two primary variables are declared:
- `porepressure` starts at 1 MPa everywhere (hydrostatic equilibrium before
  the inlet BC is applied).
- `temperature` starts at 300 K everywhere (ambient, before hot injection).

The `scaling = 1e-6` for temperature is important: pressure values are O(1e6)
while temperature values are O(300), a factor of ~3300 difference in magnitude.
Without scaling, the nonlinear solver would see a poorly-conditioned Jacobian.
The scaling multiplies the temperature equation residual by 1e-6 internally,
bringing both equations to a similar numerical magnitude.

### `[PorousFlowBasicTHM]` Action Block

```
[PorousFlowBasicTHM]
  porepressure  = porepressure
  coupling_type = ThermoHydro
  gravity       = '0 0 0'
  fp            = simple_fluid
  multiply_by_density = true
[]
```

This Action automatically generates:
- A `PorousFlowDictator` UserObject named `dictator`
- The fluid mass balance kernel with Darcy flux
- The fluid energy balance kernel with advection and conduction
- Time derivatives for both porepressure and temperature
- All internal PorousFlow material objects linking fluid properties

Setting `coupling_type = ThermoHydro` activates both the pressure (H) and
temperature (T) equations. Using `coupling_type = Hydro` would solve pressure
alone (isothermal flow) without temperature. Setting `multiply_by_density = true`
ensures the mass balance is in terms of mass flux (rho*q) rather than volume
flux, which gives conservative mass accounting.

### `[FluidProperties]` Block

```
[FluidProperties]
  [simple_fluid]
    type                 = SimpleFluidProperties
    density0             = 1000
    viscosity            = 0.001
    thermal_expansion    = 0
    cp                   = 4186
    thermal_conductivity = 0.6
  []
[]
```

`SimpleFluidProperties` is an analytical equation of state for a slightly
compressible liquid with a constant viscosity and specific heat. With
`thermal_expansion = 0`, density does not change with temperature, simplifying
the flow problem to focus on thermal advection without buoyancy-driven flow.
The properties match liquid water at approximately 20 C.

### `[Materials]` Block

Four PorousFlow materials are required:

**Porosity** (`PorousFlowPorosityConst`): constant phi = 0.3 everywhere.

**Permeability** (`PorousFlowPermeabilityConst`): the 9-component permeability
tensor (3x3, row-major) is `'1e-12 0 0  0 1e-12 0  0 0 1e-12'` — isotropic
with k = 1e-12 m^2 (approximately 1 millidarcy, characteristic of sandstone).

**Rock thermal energy** (`PorousFlowMatrixInternalEnergy`): contributes the
`(1-phi) * rho_s * c_s * T` term to the energy balance (heat stored in rock).
Without this material, the thermal retardation factor R would equal 1 (no rock
thermal mass) and the thermal front would travel at the fluid Darcy velocity.

**Thermal conductivity** (`PorousFlowThermalConductivityIdeal`): the effective
bulk conductivity of the saturated medium, given as a 3x3 tensor. The value
2.0 W/(m K) represents a mixture of water (0.6 W/(m K)) and rock (3-4 W/(m K)).

### `[BCs]` Block

Pressure BCs establish the steady-state flow field. The 0.1 MPa pressure
difference over 3 m drives a Darcy flux of ~3.3e-8 m/s rightward.

Temperature BCs model continuous injection of hot fluid from the left
(350 K) with the right boundary held at ambient (300 K). The 50 K temperature
contrast is large enough to produce a clearly visible thermal plume in ParaView.

Top and bottom walls have the default zero-flux Neumann condition (no heat
or mass passes through them), consistent with an insulated, impermeable boundary.

### `[Executioner]` Block

```
[Executioner]
  type       = Transient
  solve_type = NEWTON
  dt         = 100
  end_time   = 5000
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
```

`solve_type = NEWTON` uses the full (exact) Jacobian provided analytically by
PorousFlow kernels. This converges faster than PJFNK for the strongly-coupled
TH system but requires that all Jacobian contributions be correct. PorousFlow
provides these automatically via automatic differentiation.

`dt = 100` s gives 50 time steps to `end_time = 5000` s. This step size keeps
the CFL-like condition for heat advection below 1 per cell, ensuring accuracy.

---

## Running the Simulation

The `PorousFlow` module is part of MOOSE's `combined` application:

```bash
# From the moose-next root:
cd modules/combined && make -j8

# Then run the case:
cd /path/to/quickstart-runs/case19-porous-flow
/path/to/modules/combined/combined-opt -i case19_porous_flow.i
```

If MUMPS is not available on your system, replace the solver options in the
`[Executioner]` block with:

```
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre    boomeramg'
```

or for a simpler option:

```
petsc_options_iname = '-pc_type'
petsc_options_value = 'lu'
```

---

## Expected Results

### Pressure Field

From the first time step onward, pressure reaches a near-steady linear gradient
from 1.1 MPa (left) to 1.0 MPa (right). Because the fluid is only slightly
compressible and the pressure BCs are fixed, the pressure distribution barely
changes from step to step. This is a verification check: if you plot
`porepressure` in ParaView at any time, you should see a smooth linear gradient
in x and uniformity in y (no y-dependence because gravity is zero).

### Temperature Field

The temperature evolution is the main observable:

```
t =    0 s:  T = 300 K everywhere (initial state)
t =  500 s:  Thin warm band at left wall (plume just beginning to enter)
t = 2000 s:  Plume extends ~0.04 m into domain; steep thermal front
t = 5000 s:  Plume extends ~0.08-0.10 m; front still well left of center
```

The plume is elongated in x and centered in y (because the Darcy flow is
uniform and horizontal). The hot front is relatively sharp (conduction is
moderate; advection dominates at this Peclet number).

### Postprocessor CSV

The `avg_T` postprocessor rises slowly from 300 K as the hot plume occupies
a growing fraction of the domain. By t = 5000 s, the domain-average temperature
is only modestly above 300 K because the plume has only penetrated a small
fraction of the domain.

`max_T` stays near 350 K from early times once the inlet BC is active,
confirming that the hot inlet boundary condition is maintained throughout.

---

## Key Takeaways

| Concept | Description |
|---------|-------------|
| `PorousFlowBasicTHM` Action | Single block sets up full coupled TH system |
| Darcy flow | Pressure-driven flow in porous media (not inertial) |
| Thermal retardation | Thermal front slower than fluid due to rock heat storage |
| Variable scaling | Essential when two variables differ by orders of magnitude |
| `PorousFlowDictator` | Central registry for PorousFlow variable bookkeeping |
| `SimpleFluidProperties` | Analytical EOS for water-like single-phase fluid |
| Permeability tensor | 9-component (3x3) material property for anisotropic media |
| `multiply_by_density` | Conserves fluid mass (not volume) in the balance equation |
| Direct solver (LU/MUMPS) | Robust for small-to-medium coupled problems |

---

## Experiments to Try

**Experiment 1: Increase the permeability**
Change `permeability = '1e-10 0 0  0 1e-10 0  0 0 1e-10'` (100x higher).
The Darcy velocity increases 100-fold, so the thermal front advances much
faster. Reduce `dt` and `end_time` accordingly to see the front reach the
right boundary.

**Experiment 2: Enable gravity**
Change `gravity = '0 -9.81 0'` to include gravitational body force on the
fluid. For the current horizontal pressure gradient this adds a small
buoyancy component. For a vertical domain or a large density contrast
between hot and cold fluid, gravity-driven (thermosiphon) flow would emerge.

**Experiment 3: Isothermal flow only**
Change `coupling_type = ThermoHydro` to `coupling_type = Hydro` and remove
the `[temperature]` variable and temperature BCs. This solves only the
pressure-flow problem. The pressure should immediately reach the linear
steady-state gradient — a fast, simple verification of the Darcy solver.

**Experiment 4: Anisotropic permeability**
Replace the isotropic tensor with an anisotropic one:
`permeability = '1e-12 0 0  0 1e-14 0  0 0 1e-14'`
Now kx >> ky. Flow is still predominantly horizontal, but the reduced
vertical permeability causes the pressure field to develop a slight
y-curvature near the boundaries. The thermal plume remains symmetric
in y because the flow velocity field does too.

**Experiment 5: Run longer to see full plume breakthrough**
Change `end_time = 200000` and `dt = 2000`. The thermal front should
reach the outlet at roughly t ~ 91000 * R ~ 200000 s. At breakthrough,
`avg_T` will start rising much more sharply and `max_T` at the outlet
will jump from 300 K toward 350 K. This is the classic thermal
breakthrough curve used to characterize geothermal or tracer tests.
