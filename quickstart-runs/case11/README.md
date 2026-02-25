# Case 11: Adaptive Time Stepping with IterationAdaptiveDT

## Overview

This case demonstrates adaptive time stepping, where the simulation automatically
adjusts how large each time step is based on how hard the nonlinear solver had to
work during the previous step. If the solver converged quickly (few Newton iterations),
the problem is "easy" and we can safely take a larger step next time. If the solver
struggled (many Newton iterations), the problem is "hard" and we should take a smaller
step to stay accurate and avoid divergence.

The physical problem is a 2D transient heat equation with a step-change heat source
applied at t=0. This creates a strong initial transient (the temperature field rises
quickly at first) followed by an approach to a nearly steady state (temperature changes
slowly). A fixed time step sized for the stiff initial phase would waste thousands of
unnecessary tiny steps during the smooth final phase. Adaptive time stepping fixes this
automatically.

This case also introduces `MatDiffusion` (diffusion with a material-defined coefficient),
`BodyForce` as a heat source, and diagnostic `Postprocessors` for monitoring the
simulation as it runs.

---

## The Physics

### Physical Problem in Plain English

Imagine a square metal plate (1m x 1m) with all four edges held at zero temperature
(like a heat sink clamping the plate). At time t=0, a volumetric heat source switches
on inside the plate — think of resistive heating or a neutron flux turning on suddenly.
We want to track how the temperature field evolves from zero (initial condition) to its
final steady state where heat generation equals heat conducted to the walls.

The initial moment after the source turns on is physically dramatic: the center of the
plate heats up rapidly because the heat has not yet had time to diffuse to the cold
walls. The temperature gradient evolves quickly. Later, the solution approaches steady
state and changes slowly. Adaptive time stepping exploits this: use small dt early,
large dt late.

### Governing Equation

The transient heat equation:

    rho * cp * dT/dt = div( k * grad(T) ) + Q

Symbol explanations:

- T           — temperature field [K or degC]
- t           — time [s]
- rho         — density [kg/m^3]
- cp          — specific heat capacity [J/(kg K)]
- k           — thermal conductivity [W/(m K)]
- Q           — volumetric heat source [W/m^3]
- div         — divergence operator
- grad(T)     — temperature gradient vector

In this case, `rho * cp = 1` (implicitly, since `TimeDerivative` adds dT/dt without a
coefficient) and `k = 1.0` (from the material), so the equation simplifies to:

    dT/dt = div( grad(T) ) + 100

The source Q = 100 W/m^3 is large relative to the domain and conductivity, which is
why the initial transient is fast and nonlinear solver work is high at early times.

### Boundary and Initial Conditions

| Boundary      | Type      | Value | Physical Meaning                         |
|---------------|-----------|-------|------------------------------------------|
| left          | Dirichlet | 0     | Edge held at zero temperature (heat sink)|
| right         | Dirichlet | 0     | Edge held at zero temperature (heat sink)|
| top           | Dirichlet | 0     | Edge held at zero temperature (heat sink)|
| bottom        | Dirichlet | 0     | Edge held at zero temperature (heat sink)|
| Initial cond. | (default) | 0     | Plate starts uniformly at T = 0          |

All four walls act as perfect heat sinks. MOOSE initializes all nodal T values to zero
by default (the `[Variables]` block default initial condition is zero).

### ASCII Diagram of the Domain

```
y=1    T=0 (Dirichlet, heat sink)
       +------------------------------------+
       |                                    |
       |   dT/dt = div(k grad T) + Q        |
  T=0  |                                    |  T=0
  (D)  |   k = 1.0,  Q = 100 W/m^3         |  (D)
       |                                    |
       |   Plate starts at T = 0 everywhere |
       |   Source switches on at t = 0      |
       |                                    |
       |   Mesh: 30 x 30 elements           |
       +------------------------------------+
y=0    T=0 (Dirichlet, heat sink)
      x=0                                 x=1

Evolution:
  t=0+  (small dt): Temperature rises fast near center.
                    Newton solver works hard: many iterations.
                    --> dt stays small

  t>0.3 (large dt): Temperature nearly steady.
                    Newton solver converges easily: few iterations.
                    --> dt grows to large values
```

### Why Fixed dt Is Wasteful

Suppose you use a fixed time step sized for accuracy at early times, say dt = 0.001.
To reach t = 1.0, you need 1000 steps. But after t = 0.3, the solution barely changes
per step — you could safely use dt = 0.05 or even dt = 0.1 without losing accuracy.

Using adaptive time stepping, dt starts at 0.001 and can grow to dt = 0.064 or more
by the end. Instead of 1000 steps, you might need only ~30-50. The savings in solve
time are proportional to the number of steps avoided.

### How Iteration Count Guides dt

MOOSE's `IterationAdaptiveDT` monitors how many Newton (nonlinear) iterations the
solver used during each time step. The logic:

```
  if (actual_iterations < optimal_iterations - iteration_window):
      # Solver converged very easily -- step too small, grow it
      dt_next = dt_current * growth_factor

  elif (actual_iterations > optimal_iterations + iteration_window):
      # Solver struggled -- step too large, shrink it
      dt_next = dt_current * cutback_factor

  else:
      # Within the acceptable window -- step is about right, keep it
      dt_next = dt_current
```

This is a feedback-control mechanism: Newton iteration count is the "sensor", and dt
is the "actuator". The target is to converge in approximately `optimal_iterations`
Newton iterations every step.

---

## Input File Walkthrough

### `[Mesh]` Block

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]
```

A fixed 30x30 quadrilateral mesh on the unit square. Unlike Case 10, this mesh does
not adapt spatially — all refinement here is in the time dimension. The 30x30 grid
gives 900 elements and 961 nodes, sufficient to resolve the smooth temperature field
that develops over time.

Note the difference from Case 10's syntax: here we use `type = GeneratedMesh` directly
in the `[Mesh]` block rather than the `[gen]` sub-block syntax. Both are valid; this
older syntax is still widely used.

### `[Variables]` Block

```
[Variables]
  [T]
  []
[]
```

A single nodal variable T (temperature). Default: FIRST-order LAGRANGE basis functions.
The initial condition defaults to zero, matching the physical setup (plate starts cold).

### `[Kernels]` Block

```
[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []
  [diffusion]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []
  [source]
    type     = BodyForce
    variable = T
    value    = 100.0
  []
[]
```

Three kernels work together to build the weak form of the heat equation:

**`TimeDerivative`**: Contributes `integral( dT/dt * v ) dOmega` to the residual. This
is the time rate of change term. MOOSE uses backward Euler time integration by default
for the Transient executioner.

**`MatDiffusion`**: Contributes `-integral( k * grad(T) . grad(v) ) dOmega`. This is
the heat conduction term. Unlike the `Diffusion` kernel (which uses a hardcoded
coefficient of 1), `MatDiffusion` reads the conductivity from a material property
named `k`. The material property `k = 1.0` is defined in the `[Materials]` block.

**`BodyForce`**: Contributes `integral( Q * v ) dOmega` where Q is the `value`
parameter. With `value = 100.0`, this adds a uniform 100 W/m^3 heat source to every
element in the domain. There is no time dependence in this source — it is on from
t=0 onward (a step function), creating the strong initial transient.

### `[BCs]` Block

```
[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]
```

A single `DirichletBC` object applied to all four boundaries simultaneously. The
`boundary` parameter accepts a space-separated list of boundary names, so one object
handles all four walls. This is more concise than writing four separate BC objects.

All walls are held at T = 0 for all time. Combined with the initial condition of
T = 0 everywhere, the plate starts at steady state (T=0 everywhere with zero source),
then the source switches on and drives the system to a new steady state.

### `[Materials]` Block

```
[Materials]
  [props]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '1.0'
  []
[]
```

`GenericConstantMaterial` creates material properties with fixed scalar values.
Here, it defines the thermal conductivity `k = 1.0`. The `MatDiffusion` kernel
reads this property by name at every quadrature point during assembly.

Why use a material property instead of hardcoding the diffusivity in the kernel?
Because materials can be spatially varying (different values in different subdomains),
temperature-dependent (nonlinear), or derived from complex expressions. Starting with
`GenericConstantMaterial` is the simplest case, and switching to a more complex
material later requires no change to the kernel.

### `[Postprocessors]` Block

```
[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []
  [dt_pp]
    type = TimestepSize
  []
[]
```

Three diagnostic quantities computed at every time step:

**`ElementExtremeValue`** with `value_type = max`: Finds the maximum value of T across
all quadrature points in all elements. This tracks the hottest point in the domain
(which will be near the center). Useful for checking against safety limits.

**`ElementAverageValue`**: Computes the volume-weighted average of T:
`avg_T = (1/|Omega|) * integral(T) dOmega`. For a unit square, |Omega| = 1, so this
is just the integral of T. Tracks how the overall energy content grows toward steady
state.

**`TimestepSize`**: Reports the current dt value. Critically useful for verifying that
the adaptive time stepper is actually growing dt as expected. Watching this column in
the CSV output is the best way to confirm the adaptive algorithm is working.

### `[Executioner]` Block

```
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 5
    growth_factor = 2.0
    cutback_factor = 0.5
    iteration_window = 2
  []

  start_time = 0.0
  end_time   = 1.0
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
```

The `Transient` executioner drives the time-stepping loop. Key parameters:

**`start_time = 0.0`**, **`end_time = 1.0`**: The simulation runs for 1 second of
physical time. With the given parameters, this is long enough for the solution to
closely approach steady state.

**`nl_rel_tol = 1e-8`**: Newton iteration converges when the nonlinear residual is
reduced by a factor of 10^8 from its initial value. This tight tolerance means Newton
will work harder when the problem is stiff (early times) and less hard when it is
smooth (later times) — exactly the signal that drives the adaptive time stepper.

#### `[TimeStepper]` Sub-block — the Heart of This Case

**`type = IterationAdaptiveDT`**: Selects the iteration-count-driven adaptive time
stepper.

**`dt = 0.001`**: The initial time step size at t=0. This is small enough to accurately
capture the rapid initial transient when the source first switches on.

**`optimal_iterations = 5`**: Target Newton iteration count per time step. If the
solver converges in exactly 5 Newton iterations, the time step is considered "right-
sized" and is not changed.

**`growth_factor = 2.0`**: When the solver converges in fewer than
`optimal_iterations - iteration_window = 5 - 2 = 3` iterations, dt is multiplied by
2.0. This aggressive doubling allows dt to grow rapidly once the solution smooths out.

**`cutback_factor = 0.5`**: When the solver takes more than
`optimal_iterations + iteration_window = 5 + 2 = 7` iterations, dt is multiplied by
0.5. This halving prevents divergence when the step is too large.

**`iteration_window = 2`**: An acceptable band of +/- 2 iterations around the target.
If actual iterations are in the range [3, 7] = [5-2, 5+2], dt is not changed. This
prevents constant oscillation of dt when the solver count fluctuates by 1-2 iterations.

#### How the Algorithm Plays Out

```
t=0.000, dt=0.001:  Source just switched on. Large gradients.
                    Newton needs ~8 iterations.
                    8 > optimal+window (7) --> cutback: dt *= 0.5 = 0.0005

t=0.001, dt=0.0005: Source still strong. Newton needs ~7 iterations.
                    7 == optimal+window --> no change: dt = 0.0005

t=0.002, dt=0.0005: ...Solution growing smoothly...
                    Newton needs ~4 iterations.
                    4 < optimal-window (3)? No (4 >= 3). --> no change.

t=...              ...Eventually Newton drops to ~2 iterations...
                    2 < 3 --> grow: dt *= 2.0

t=0.010, dt=0.001:  Newton still fast --> grow again: dt = 0.002
t=0.012, dt=0.002:  Newton still fast --> grow: dt = 0.004
...
t=0.100, dt=0.032:  Solution nearly steady. Newton converges in 1 step.
                    --> dt grows to 0.064, then 0.128, etc.
```

The final time steps near t=1.0 may have dt = 0.256 or larger — more than 100x
larger than the initial dt = 0.001. The entire simulation completes in tens of steps
rather than thousands.

### `[Outputs]` Block

```
[Outputs]
  exodus = true
  csv    = true
[]
```

Two output formats:

**`exodus = true`**: Writes the spatial solution (the T field at every time step) to an
Exodus file. You can animate this in ParaView to see the temperature field evolving.

**`csv = true`**: Writes the postprocessor values (max_T, avg_T, dt_pp) to a CSV file
with one row per time step. This is the best way to verify adaptive time stepping.

---

## What Happens When You Run This

Run the case with:

```bash
./moose_test-opt -i case11_adaptive_dt.i
```

The console output shows each time step on one line. Watch for the dt column and the
Newton iteration count. You will see something like:

```
Time Step 1, time = 0.001, dt = 0.001
  0 Nonlinear |R| = 2.8e+01
  1 Nonlinear |R| = 1.2e+00
  ...
  6 Nonlinear |R| = 3.1e-10  (converged)
  (6 iterations > optimal+window=7? No. Borderline.)

Time Step 2, time = 0.0015, dt = 0.0005
  ...converged in 5 iterations...
  (5 is optimal, no change)

...many steps later...

Time Step 25, time = 0.35, dt = 0.064
  0 Nonlinear |R| = 1.4e-01
  1 Nonlinear |R| = 8.2e-09  (converged in 1 iteration!)
  --> dt grows: 0.064 * 2.0 = 0.128

Time Step 26, time = 0.48, dt = 0.128
  ...
```

The postprocessor table printed to console at each step shows:

```
+----------------+-------+-------+-------+
| time           | max_T | avg_T | dt_pp |
+----------------+-------+-------+-------+
| 0.001          | 0.098 | 0.019 | 0.001 |
| 0.0015         | 0.147 | 0.028 | 0.0005|
| ...            | ...   | ...   | ...   |
| 1.000          | 6.21  | 2.48  | 0.256 |
+----------------+-------+-------+-------+
```

You can see max_T growing from 0 toward its steady-state value, avg_T rising similarly,
and dt_pp growing from 0.001 to much larger values as the simulation proceeds.

---

## Output Files

| File | Description |
|------|-------------|
| `case11_adaptive_dt_out.e` | Exodus file: T field at every time step |
| `case11_adaptive_dt_out.csv` | CSV: time, max_T, avg_T, dt_pp at each step |

### The CSV File

The CSV file (`case11_adaptive_dt_out.csv`) has this structure:

```
time,max_temp,avg_temp,dt_pp
0,0,0,0
0.001,0.09812...,0.01922...,0.001
0.0015,0.14718...,0.02883...,0.0005
...
```

Columns:
- `time`: Physical simulation time at the end of this step
- `max_T`: Maximum temperature over the whole domain at this time
- `avg_T`: Domain-average temperature at this time
- `dt_pp`: The time step size used for this step

The `dt_pp` column is the key diagnostic. Plot it vs. time (or even just scan it by
eye in the CSV) to confirm that dt grows over time. In a well-functioning adaptive
run, dt should grow by roughly 2 orders of magnitude from start to finish.

### Visualizing in ParaView

1. Open `case11_adaptive_dt_out.e` in ParaView.
2. Click Apply, color by T.
3. Use the "Play" button to animate. The early time steps will show rapid change
   (many frames); later steps will show slow, nearly imperceptible change even though
   each frame covers a much larger time interval.
4. Note: the time slider steps through each time step. With adaptive stepping, the
   time between frames is not constant.

---

## Interpreting the Results

### Physical Description

The temperature field evolves in three phases:

1. **Early transient (t = 0 to ~0.05)**: Heat source drives temperatures up rapidly.
   The pattern is symmetric about the center (x=0.5, y=0.5) because all BCs are
   symmetric. The center rises fastest; wall temperatures are fixed at zero.

2. **Mid transient (t = 0.05 to ~0.3)**: Temperature profile spreads. The center
   reaches its peak rate of rise and starts slowing. dt is now growing rapidly.

3. **Approach to steady state (t = 0.3 to 1.0)**: The solution changes very slowly.
   Large time steps are safe. The solver converges in 1-2 Newton iterations per step.

### Steady-State Reference

At steady state (dT/dt = 0), the equation becomes:

    div( grad(T) ) + 100 = 0     with T = 0 on all walls

The exact solution is a Fourier series. The maximum temperature at the center (0.5, 0.5)
at steady state for k=1, Q=100, unit square domain is:

    T_center_ss = Q * 16 / (k * pi^4) * sum_{odd m, n} 1/(m*n*(m^2+n^2)*pi^2/4)

Approximately: T_max_ss ~ 6.0-6.5

By t=1.0, the simulation should show max_T approaching this value.

### Correctness Verification

Three checks:
1. The temperature field should be symmetric: T(x, y) = T(1-x, y) = T(x, 1-y)
2. max_T should approach ~6.2 by t=1.0
3. The dt_pp column in the CSV should show monotonically growing values (with
   possible brief cutbacks if the solver momentarily struggles)

---

## Key Concepts Learned

| Concept | What It Is | Where in Input File |
|---------|-----------|---------------------|
| Transient executioner | Time-dependent simulation loop | `type = Transient` |
| IterationAdaptiveDT | Adjusts dt based on Newton iteration count | `[TimeStepper]` block |
| optimal_iterations | Target Newton iteration count | `optimal_iterations = 5` |
| growth_factor | dt multiplier when solver is fast | `growth_factor = 2.0` |
| cutback_factor | dt multiplier when solver is slow | `cutback_factor = 0.5` |
| iteration_window | Tolerance band around optimal | `iteration_window = 2` |
| MatDiffusion | Diffusion kernel with material-defined coefficient | `[Kernels]` |
| BodyForce | Constant volumetric source term | `[Kernels]` |
| GenericConstantMaterial | Defines material properties as constants | `[Materials]` |
| TimestepSize postprocessor | Reports current dt to CSV/console | `[Postprocessors]` |

---

## Experiments to Try

**Experiment 1: Change the initial dt**
In the `[TimeStepper]` block, change `dt = 0.001` to `dt = 0.01`. The early time steps
will be 10x larger. Does the solver still converge? Does it take more Newton iterations
at early times? How does this change the final total number of time steps?

**Experiment 2: Adjust growth_factor**
Change `growth_factor = 2.0` to `growth_factor = 1.5`. The time step will grow more
slowly, requiring more steps to reach t=1.0. Is the accuracy better? Then try
`growth_factor = 4.0` — does the solver ever fail to converge because dt grew too fast?

**Experiment 3: Narrow the iteration_window**
Change `iteration_window = 2` to `iteration_window = 0`. Now the window is exact:
only exactly `optimal_iterations = 5` Newton iterations avoids a change. Does dt
oscillate up and down more? Is the simulation less efficient?

**Experiment 4: Change the source value**
Change `value = 100.0` (BodyForce) to `value = 10.0`. The transient is gentler. Does
the initial dt cutback happen less? Does the adaptive stepper grow dt faster because
the Newton solver converges more easily from the start?

**Experiment 5: Extend end_time**
Change `end_time = 1.0` to `end_time = 10.0`. By t=10 the solution is truly at
steady state. Observe how dt grows to very large values in the late phase. How many
total time steps does the simulation take to cover 10 seconds? Compare to the ~1000
steps a fixed dt=0.001 run would need.
