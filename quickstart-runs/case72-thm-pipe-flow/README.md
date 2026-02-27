# Case 72: THM Pipe Flow — 1D Single-Phase Compressible

## Overview

The thermal_hydraulics module (THM) in MOOSE provides a component-based domain-specific language (DSL) for modeling thermal-hydraulic systems: networks of pipes, heat exchangers, pumps, valves, and boundary condition components connected to form a system. Rather than writing the Variables, Kernels, and BCs blocks manually, the user declares physical components and the module automatically assembles the governing equations, discretization, and coupling terms. This approach is inspired by the Modelica language philosophy and is particularly well-suited for nuclear thermal-hydraulic safety analysis.

This case models 1D compressible single-phase flow of an ideal gas through a pipe with significant friction. The fluid properties are those of a diatomic ideal gas (gamma = 1.4, like air or nitrogen). An inlet supplies a fixed mass flow rate of 2 kg/s at 500 K; an outlet holds the exit pressure at 2e5 Pa. The system starts far from steady state (T = 300 K, p = 1e5 Pa throughout), and the simulation captures the transient as the flow accelerates, pressure adjusts to meet the boundary conditions, and the temperature front propagates from inlet to outlet under the influence of advection and frictional dissipation.

Key concepts demonstrated:

- `[FluidProperties]` block with `IdealGasFluidProperties` for the equation of state
- `[Closures]` block with `Closures1PhaseSimple` for friction and heat transfer correlations
- `[Components]` block as the THM Component DSL: `FlowChannel1Phase`, `InletMassFlowRateTemperature1Phase`, `Outlet1Phase`
- `scaling_factor_1phase` for conditioning the compressible flow linear system
- `PointValue` postprocessors to track spatial temperature and pressure profiles along the pipe

---

## The Physics

### 1D Compressible Euler Equations with Friction

The THM module solves the area-averaged 1D Euler equations for mass, momentum, and energy:

```
d(A*rho)/dt + d(A*rho*vel)/dx = 0

d(A*rho*vel)/dt + d(A*(rho*vel^2 + p))/dx = -f * rho*vel^2 * A / (2*D_h)

d(A*rho*E)/dt + d(A*vel*(rho*E + p))/dx = -f * rho*vel^3 * A / (2*D_h)
```

where:
- A = 1 m^2 — pipe cross-sectional area
- rho — fluid density
- vel — fluid velocity
- p — pressure
- E = e + vel^2/2 — specific total energy (e = specific internal energy)
- f = 10 — Darcy friction factor (high value for visible demonstration effect)
- D_h = A^(1/2) = 1 m — hydraulic diameter

The three equations are solved simultaneously for the conservative variables (rho*A, rho*vel*A, rho*E*A). Primitive variables (T, p, vel) are recovered from these via the equation of state.

### Ideal Gas Equation of State

The `IdealGasFluidProperties` object computes thermodynamic properties from first principles:

```
p = rho * R_specific * T           [equation of state]
e = cv * T                          [internal energy]
cv = R_specific / (gamma - 1)       [specific heat at constant volume]
```

For gamma = 1.4 (diatomic ideal gas): the ratio cv/R = 1/(gamma-1) = 2.5 and cp/R = gamma/(gamma-1) = 3.5. The default MOOSE `IdealGasFluidProperties` uses air-like parameters with molecular weight M = 28.97 g/mol, giving R_specific = R/M = 287 J/(kg*K), cv = 717.5 J/(kg*K), and cp = 1004.5 J/(kg*K).

### Acoustic vs. Convective Timescales

Two very different timescales govern the transient:

1. **Acoustic timescale**: L/c where c = sqrt(gamma * R * T) ~ sqrt(1.4 * 287 * 300) ~ 347 m/s. For L = 1 m: t_acoustic ~ 3 ms.
2. **Convective timescale**: L/vel where vel ~ m_dot/(rho*A) ~ 2/(1.2*1) ~ 1.7 m/s. For L = 1 m: t_convective ~ 0.6 s.

Pressure adjusts on the acoustic timescale (3 ms << dt = 25 ms), so the pressure field reaches its new equilibrium almost instantaneously. Temperature propagates on the convective timescale, so the temperature front is still traversing the pipe at t = 0.5 s.

### Domain

- Pipe: L = 1 m, A = 1 m^2, oriented along the x-axis, no gravity
- 50 finite volume elements (dx = 0.02 m each)
- Time: t in [0, 0.5 s], dt = 0.025 s (20 time steps)

---

## Input File Walkthrough

The input file is `case72_thm_pipe_flow.i`.

### `[FluidProperties]`

```
[fp]
  type = IdealGasFluidProperties
  gamma = 1.4
[]
```

Declares a fluid properties object named `fp`. The THM components that reference `fp = fp` will use this equation of state for all thermodynamic calculations including density, internal energy, enthalpy, speed of sound, and transport properties. No Mesh, Variables, Kernels, or BCs blocks are written — the Component DSL generates all of these automatically based on the `[Components]` block.

### `[Closures]`

```
[simple_closures]
  type = Closures1PhaseSimple
[]
```

`Closures1PhaseSimple` provides the simplest possible closure relations: a constant (user-specified) friction factor and no wall heat transfer (adiabatic pipe). More advanced closures such as `Closures1PhaseTHM` provide the Churchill friction correlation (covering laminar, transitional, and turbulent regimes automatically) and Dittus-Boelter heat transfer. The `simple` closure is appropriate for prototype problems and cases where friction is explicitly prescribed.

### `[Components]` — The THM DSL

The three components form the complete system:

```
[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe:in'
  m_dot = 2
  T = 500
[]
```

`InletMassFlowRateTemperature1Phase` prescribes the inlet mass flow rate (m_dot = 2 kg/s) and total temperature (T_in = 500 K). The `input = 'pipe:in'` parameter specifies which end of which component this inlet connects to. The boundary condition is implemented as a characteristic-based inflow condition consistent with the hyperbolic character of the Euler equations.

```
[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  gravity_vector = '0 0 0'
  length = 1.0
  n_elems = 50
  A = 1.0
  initial_T = 300
  initial_p = 1e5
  initial_vel = 1
  f = 10.0
  closures = simple_closures
  fp = fp
  scaling_factor_1phase = '1 1 1e-5'
[]
```

`FlowChannel1Phase` is the core component: a 1D finite volume pipe. `position` and `orientation` define the pipe geometry in 3D space. `gravity_vector = '0 0 0'` disables gravitational body force. The three initial conditions (initial_T, initial_p, initial_vel) define the starting state. `f = 10.0` is the Darcy friction factor passed to the closures. `closures = simple_closures` and `fp = fp` reference the declared closures and fluid properties objects.

```
[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = 2e5
[]
```

`Outlet1Phase` prescribes the outlet back pressure (p_out = 2e5 Pa). The boundary condition permits supersonic outflow (all characteristics exit) or subsonic outflow (one characteristic enters, carrying the pressure information).

### `scaling_factor_1phase`

```
scaling_factor_1phase = '1 1 1e-5'
```

This is a critical setting for numerical conditioning. The three entries scale the three conservative variables: mass density (rho*A), momentum density (rho*vel*A), and total energy density (rho*E*A). The energy variable has units of J/m and values ~ 2e5 J/m (much larger than mass density ~ 1.4 kg/m and momentum ~ 2.8 kg/(m*s)). Without scaling, the Jacobian is ill-conditioned by a factor of ~1e5 and Newton iterations converge slowly or diverge. The factor 1e-5 on the energy equation reduces it to the same order of magnitude as the mass and momentum equations.

### `[Postprocessors]`

```
[T_near_inlet]   type = PointValue  variable = T  point = '0.01 0 0'
[T_mid]          type = PointValue  variable = T  point = '0.5 0 0'
[T_near_outlet]  type = PointValue  variable = T  point = '0.99 0 0'
[p_near_inlet]   type = PointValue  variable = p  point = '0.01 0 0'
[p_near_outlet]  type = PointValue  variable = p  point = '0.99 0 0'
```

`PointValue` postprocessors sample the THM-generated variables T and p at specific spatial locations. The THM module creates these variables automatically from the component definitions — the user never declares `[Variables]` for T or p. Sampling near (but not exactly at) the inlet/outlet boundaries (x = 0.01 and x = 0.99) avoids boundary layer artifacts at the component interfaces.

### `[Executioner]`

```
type = Transient
solve_type = NEWTON
end_time = 0.5
dt = 0.025
petsc_options_iname = '-pc_type'
petsc_options_value = 'lu'
[Quadrature]
  type = GAUSS
  order = SECOND
[]
```

NEWTON solver (as opposed to PJFNK) uses the exact analytic Jacobian provided by the THM module. Direct LU factorization handles the well-conditioned (after scaling) 3-variable-per-element system efficiently for 50 elements (150 unknowns total). The `GAUSS SECOND` quadrature is appropriate for the piecewise-constant (first-order finite volume) elements used by the THM discretization.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case72-thm-pipe-flow \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case72_thm_pipe_flow.i 2>&1 | tail -30'
```

Output files:
- `case72_thm_pipe_flow_out.e` — Exodus with T, p, vel, rho fields along the pipe over time
- `case72_thm_pipe_flow_out.csv` — Time series of temperatures and pressures at three axial locations

---

## Expected Results

### Transient Evolution

| Time (s) | T_near_inlet (K) | T_mid (K) | T_near_outlet (K) | p_near_inlet (Pa) | p_near_outlet (Pa) |
|----------|------------------|-----------|-------------------|-------------------|--------------------|
| 0.000    | 300              | 300       | 300               | 100,000           | 100,000            |
| 0.025    | 461              | 387       | 394               | 199,224           | 199,769            |
| 0.050    | 486              | 388       | 394               | 200,532           | 200,020            |
| 0.075    | 495              | 387       | 393               | 200,026           | 200,001            |
| 0.100    | 498              | 387       | 393               | 200,014           | 200,000            |
| 0.200    | 500              | 391       | 391               | 200,015           | 200,000            |
| 0.300    | 500              | 415       | 389               | 200,015           | 200,000            |
| 0.400    | 500              | 453       | 389               | 200,014           | 200,001            |
| 0.500    | 500              | 483       | 395               | 200,014           | 200,000            |

### Pressure Adjustment (Nearly Instantaneous)

The pressure jumps from 1e5 Pa to approximately 2e5 Pa within the very first time step (t = 0.025 s). Acoustic pressure waves propagate through the gas at c ~ 347 m/s, so the 1 m pipe is traversed in ~3 ms — far shorter than the 25 ms time step. After the first step, p_near_outlet tracks the outlet boundary condition at exactly 2e5 Pa. The small elevation at p_near_inlet (~200,015 Pa vs. 200,000 Pa at the outlet) is the friction pressure drop: delta_p = f * rho * vel^2 * L / (2 * D_h) ~ 10 * 1.4 * 1.4^2 / 2 ~ 14 Pa.

### Temperature Front Propagation

The temperature response is much slower than the pressure response. T_near_inlet rises from 300 K toward the 500 K inlet boundary condition within the first 0.1 s (limited by the time for the inlet BC to propagate one element inward). T_mid then rises progressively as the hot fluid convects from inlet toward outlet. By t = 0.5 s, T_mid has reached ~483 K — the thermal front has nearly traversed the full pipe length. The convective timescale L/vel ~ 1/(2/1.4) ~ 0.7 s is consistent with this behavior.

---

## Key Takeaways

- The `[Components]` block is the defining feature of the THM module. Instead of writing Variable, Kernel, BC, and InterfaceKernel blocks, the user declares physical components (pipes, inlets, outlets, heat exchangers) and the module assembles the complete equation system automatically. This reduces input file complexity by an order of magnitude for large networks with many pipe segments and junctions.
- `scaling_factor_1phase = '1 1 1e-5'` is essential for the compressible flow equations. The three conservative variables — mass density, momentum density, and total energy density — have vastly different magnitudes for typical gas flow conditions. Without this scaling, the Jacobian condition number is approximately 1e5 and Newton iteration will converge poorly or diverge. Always set the energy scaling factor to bring all three equations to the same order of magnitude.
- Pressure adjusts on the acoustic timescale (L/c ~ milliseconds) while temperature adjusts on the convective timescale (L/vel ~ seconds). These timescales differ by three orders of magnitude. The implicit time integration used by THM handles both stably with a single time step size, at the cost of solving a larger coupled system per step compared to explicit methods.
- `IdealGasFluidProperties` is a simple, general equation of state for compressible gas flow. For liquid-dominated systems (nuclear reactor coolant loops, steam generators), `StiffenedGasFluidProperties` or tabulated EOS objects are available in MOOSE's fluid properties library. The THM components are agnostic to the specific EOS — swapping `fp` objects changes the physics without modifying the component or executioner setup.
- `Closures1PhaseSimple` uses the user-prescribed friction factor f directly as the Darcy friction factor. For engineering applications, `Closures1PhaseTHM` provides the Churchill correlation for friction (covering laminar, transitional, and turbulent regimes automatically) and the Dittus-Boelter correlation for wall heat transfer. The closures system separates physics modeling (correlations) from the governing equation structure.
- `InletMassFlowRateTemperature1Phase` is the preferred inlet BC for forced-flow systems where a pump or compressor sets the mass flow rate. `InletStagnationPressureTemperature1Phase` is used when the inlet stagnation conditions are prescribed instead. The choice of inlet type determines which characteristics are specified at the boundary and must be consistent with whether the flow is subsonic or supersonic.
- THM scales naturally from this single-pipe prototype to full system models: multiple pipe segments connected by junctions, heat exchangers coupling primary and secondary loops, pumps and turbines, pressurizers, and steam generators. The same component-based input structure applies regardless of system complexity, making THM the primary tool for nuclear thermal-hydraulic safety analysis within the MOOSE ecosystem.
