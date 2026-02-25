# Case 03: Transient Heat Equation with Source Term

## Overview

This case introduces three major new ideas simultaneously, all of which appear in
nearly every real engineering simulation:

1. **Transient (time-dependent) solving**: the solution evolves over time, not just
   at a single steady state. MOOSE marches forward in time using a time-stepping scheme.

2. **A volumetric heat source**: an internal forcing term `Q` adds energy to the domain
   at every point continuously. Without it, the solution would be trivial. With it, the
   solution evolves toward a non-trivial steady state.

3. **Material properties**: thermal conductivity `k` is now defined through the
   `[Materials]` system, making it easy to vary `k` spatially, temporally, or even
   nonlinearly in more advanced problems.

Additionally, this case introduces two new MOOSE object types:

- **Postprocessors**: scalar quantities computed from the solution at each timestep —
  here, the spatial average temperature and the peak temperature.
- **`MatDiffusion` kernel**: a generalization of `Diffusion` that reads a material
  property `k` instead of assuming unit diffusivity.

What you will learn:
- How to set up a time-dependent solve with the `Transient` executioner
- What a `TimeDerivative` kernel is and how it couples the time axis to the PDE
- How the `[Materials]` block works and why material properties are separate from kernels
- How to define and use `[Postprocessors]` for scalar time-history outputs
- How to read and interpret a CSV time-history file
- What "approach to steady state" looks like numerically: the average temperature rises
  quickly at first, then asymptotically approaches a constant limit
- How convergence tolerances (`nl_rel_tol`, `nl_abs_tol`) control solve quality

The exact long-time solution (steady state) for this problem is known analytically and
can be used to verify the simulation at late times.

---

## The Physics

### The Physical Problem in Plain English

Imagine the same flat square plate from Case 02: one meter on each side, all four walls
held at exactly 0 degrees (cold walls). But now, instead of being passive, the plate
contains a distributed electrical heater that deposits 1 watt per cubic meter uniformly
throughout the interior.

At time t=0 the plate is everywhere at 0 degrees (cold start). Then we switch on the
heater. Heat begins to accumulate in the plate, raising the temperature. But as the
plate warms, the temperature difference between the plate interior and the cold walls
grows, which drives heat out through the walls faster. Eventually, heat generated
inside exactly balances heat lost through the walls, and the temperature stops rising.
This equilibrium is the steady state.

This problem models many real engineering situations:
- A nuclear fuel pellet generating heat from fission, cooled at its surface
- An electronic chip dissipating power, cooled by a heat sink at its edges
- A chemical reaction vessel with exothermic reaction, cooled at the walls
- A geothermal body with radiogenic heat generation, cooled at its surface

The transient behavior — how fast the temperature rises and how it distributes spatially
as it approaches steady state — is the information this case provides.

### The Governing Equation

The transient heat equation (also called the diffusion equation or parabolic PDE):

```
         dT
rho*cp * -- = div(k * grad T) + Q    on [0,1]x[0,1], t > 0
         dt
```

Written out fully in two dimensions:

```
         dT     d     dT         d     dT
rho*cp * -- = --(k * --) + --(k * --) + Q(x, y, t)
         dt     dx    dx         dy    dy
```

Every symbol defined:
- `T(x, y, t)` — temperature, a function of space and time. This is the unknown we solve.
- `rho` — material density (kg/m³). How heavy the material is per unit volume.
- `cp` — specific heat capacity (J/(kg·K)). How much energy is needed to raise 1 kg by 1 K.
- `rho * cp` — volumetric heat capacity (J/(m³·K)). Energy storage per unit volume per degree.
- `dT/dt` — rate of temperature change with respect to time. How fast the plate is warming.
- `k` — thermal conductivity (W/(m·K)). How easily heat flows through the material.
- `div(k * grad T)` — the diffusion term: describes heat spreading from hot to cold regions.
- `grad T` — temperature gradient vector: points in direction of steepest temperature rise.
- `k * grad T` — heat flux vector (Fourier's law): heat flows opposite to the gradient.
- `div(...)` — divergence of flux: net heat arriving at a point minus heat leaving it.
- `Q` — volumetric heat source (W/m³). Energy added per unit volume per unit time.

**Simplified parameters for this case:**

This case sets `rho = cp = k = 1` and `Q = 1` (all in consistent dimensionless units):

```
dT/dt = div(grad T) + 1    on [0,1]x[0,1], t > 0
```

These choices are made to keep the input file simple while preserving all the important
physics. In a real problem you would set these to the actual material properties of your
substance.

**Why `rho*cp*dT/dt` and not just `dT/dt`?**

The product `rho*cp` is the thermal mass of the material — how much energy it can store
per degree of temperature rise per unit volume. A large `rho*cp` means the material heats
up slowly (like water or concrete). A small `rho*cp` means it heats up quickly (like air).
With `rho = cp = 1`, we get unit thermal mass, so temperature changes on the same time
scale as the heat input.

### Boundary Conditions

All four walls are held at zero temperature (Dirichlet conditions):

```
T = 0    on all four boundaries:
         x=0 (left), x=1 (right), y=0 (bottom), y=1 (top)
```

Physical meaning: the walls are perfect conductors connected to a large heat reservoir
at 0 degrees. No matter how much heat accumulates in the interior, the walls stay cold.
This is an idealization (a simplification of a problem where the walls have a very large
cooling capacity).

### Initial Condition

```
T(x, y, 0) = 0    everywhere in the domain
```

The plate starts cold. At t=0 we simultaneously switch on the heater and begin the
simulation.

**What is "transient" vs "steady state"?**

- **Steady state**: the solution does not change with time (`dT/dt = 0` everywhere).
  The heat equation becomes the Poisson equation: `div(grad T) + Q = 0`.
- **Transient**: the solution changes with time. The `dT/dt` term in the PDE is nonzero.
  We need to track how the solution evolves from `t=0` to some final time.

### Steady-State Solution

As `t → infinity`, the transient solution approaches a steady state where
`dT/dt = 0` everywhere. The steady-state equation is:

```
-div(grad T) = Q = 1    on [0,1]x[0,1]
T = 0 on all walls
```

This is the Poisson equation (see Case 04 for a related problem). On the unit square with
`Q=1` and zero Dirichlet BCs, the exact steady-state solution is an infinite series:

```
T_ss(x, y) = sum_{m=1,3,5,...} sum_{n=1,3,5,...}
             16 / (pi^4 * m * n * (m^2 + n^2))
             * sin(m*pi*x) * sin(n*pi*y)
```

The first term (m=1, n=1) dominates: the steady-state temperature is roughly bell-shaped,
peaking at the center (x=0.5, y=0.5). The peak steady-state temperature can be computed
approximately as `T_max ≈ 0.0737` for Q=1.

The spatial average of the steady-state temperature is `T_avg = Q * L^2 / (8*pi^2) ≈ 0.0127`
for Q=1 on the unit square (this follows from integrating the series term by term).

### Time Scale

How long does it take for the plate to reach steady state? The characteristic thermal
time scale is:

```
tau = rho * cp * L^2 / (pi^2 * k) ≈ L^2 / (pi^2) ≈ 0.1 seconds
```

where L=1 is the domain size. After about `3*tau ≈ 0.3 seconds`, the solution is close
to steady state. The simulation runs to `t = 0.5 seconds`, which is about `5*tau`,
capturing essentially the full transient evolution plus approach to steady state.

### ASCII Domain Diagram

```
  y=1  T=0 (top wall, cold)
       +------------------------------------------+
       |          .   .   .   .   .               |
       |       .    warm interior    .             |
  T=0  |     .   Q=1 everywhere        .         |  T=0
(left) |    .   Heat rises from center  .        | (right)
       |   .      Peak T at (0.5, 0.5)   .       |
       |    .                           .         |
       |       .    warm interior    .             |
       |          .   .   .   .   .               |
       +------------------------------------------+
  y=0  T=0 (bottom wall, cold)

Time evolution of average temperature:
    T_avg
    ^
0.013|                              _______________
     |                          ___/
     |                      ___/
     |                  ___/
     |              ___/
     |          ___/
     |      ___/
     |   __/
     |__/
  0  +-------------------------------------------> time
    t=0                                         t=0.5

The plate heats up quickly at first, then asymptotically approaches steady state.
```

---

## Input File Walkthrough

The input file is `case03_heat_transient.i`. There are three new blocks compared to
Case 02: `[Materials]`, `[Postprocessors]`, and a modified `[Executioner]`.

### Header Comments

```
# ============================================================
# Case 3: Transient Heat Equation
# rho*cp * dT/dt = div(k*grad T) + Q
# rho=cp=k=1, Q=1, T=0 on all walls, T_0=0
# ============================================================
```

The comment documents the full PDE with all parameters, then lists the specific values
chosen for this case. The initial condition `T_0=0` is explicit. This header gives a
complete specification of the mathematical problem independent of the MOOSE syntax.

---

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]
```

Identical to Case 02 but with default domain bounds: `xmin=0, xmax=1, ymin=0, ymax=1`.
When `xmin/xmax/ymin/ymax` are omitted from a `GeneratedMesh`, MOOSE uses `[0,1]^dim`
by default. The unit square mesh is exactly the same: 400 elements, 441 nodes.

---

### Block: `[Variables]`

```
[Variables]
  [T]
  []
[]
```

The variable is now named `T` (for temperature) instead of `u`. The name is arbitrary
but `T` is clearer than `u` for a heat problem. All defaults still apply: first-order
Lagrange elements, zero initial condition.

The zero initial condition `T(x,y,0) = 0` everywhere is exactly what we want — the
plate starts at 0 degrees (the same temperature as the walls). If you needed a non-zero
initial temperature distribution, you would add:

```
[Variables]
  [T]
    initial_condition = 100    # uniform 100-degree initial temperature
  []
[]
```

Or use a `[ICs]` block with a `FunctionIC` for a spatially varying initial condition.

---

### Block: `[Kernels]`

```
[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []

  [heat_conduction]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []

  [heat_source]
    type     = BodyForce
    variable = T
    value    = 1.0
  []
[]
```

This case has three kernels, each contributing one term of the PDE. Together they
implement `dT/dt = div(grad T) + 1`.

**`[time_deriv]` — TimeDerivative kernel**

```
[time_deriv]
  type     = TimeDerivative
  variable = T
[]
```

The `TimeDerivative` kernel implements the time derivative term in the weak form:

```
R_i = integral( phi_i * dT/dt  dV )
```

This is how MOOSE couples time to the PDE. Internally, `dT/dt` is approximated using
the **backward Euler** time integration scheme (MOOSE's default for implicit transient):

```
dT/dt ≈ (T^{n+1} - T^n) / dt
```

where `T^{n+1}` is the unknown solution at the new time step and `T^n` is the known
solution at the previous time step. Backward Euler is first-order accurate in time
and unconditionally stable (it does not "blow up" for any time step size).

The kernel uses `_u_dot[_qp]` (MOOSE's notation for the time derivative at quadrature
point `_qp`) and `_du_dot_du[_qp]` for the Jacobian contribution. These are provided
automatically by MOOSE's time integration infrastructure.

**Why is this kernel separate from the spatial kernels?**

Separation of concerns: the spatial operator (`MatDiffusion`) knows nothing about
time, and the temporal operator (`TimeDerivative`) knows nothing about space. You
can mix and match: add `TimeDerivative` to any steady-state problem to make it
transient. Remove it to recover the steady-state equation. This modularity is a core
MOOSE design principle.

**`[heat_conduction]` — MatDiffusion kernel**

```
[heat_conduction]
  type        = MatDiffusion
  variable    = T
  diffusivity = k
[]
```

`MatDiffusion` implements the heat conduction term in weak form:

```
R_i = integral( k * grad(phi_i) . grad(T)  dV )
```

The key difference from `Diffusion` (Case 01/02): instead of assuming unit diffusivity,
`MatDiffusion` reads the conductivity from a **material property** named `k`. The name
`k` in `diffusivity = k` is just a string that MOOSE looks up in the material system
at each quadrature point.

Why use a material property instead of a hardcoded coefficient?
- **Flexibility**: `k` can later be a function of temperature, space, or time without
  changing the kernel
- **Separation of concerns**: the PDE operator (kernel) is separate from material
  behavior (material object)
- **Composability**: multiple kernels can read the same material property; changing the
  material object changes all kernels that use it simultaneously

The `k` property is defined in the `[Materials]` block, described below.

At each quadrature point, MOOSE evaluates `_diffusivity[_qp]` (the material property
value) and uses it in the kernel computation. This lookup is efficient — the material
property is computed once per quadrature point and cached.

**`[heat_source]` — BodyForce kernel**

```
[heat_source]
  type     = BodyForce
  variable = T
  value    = 1.0
[]
```

`BodyForce` implements a volumetric source term in weak form:

```
R_i = integral( phi_i * Q  dV )
```

`value = 1.0` sets `Q = 1.0` uniformly over the entire domain at all times.
This is the internal heat generation rate (1 W/m³ in dimensional terms).

The `BodyForce` kernel contributes to the right-hand side of the linear system:
it adds energy at every interior point every timestep. Without it, with T=0 initial
conditions and T=0 Dirichlet BCs, the solution would stay at T=0 forever.

Alternatives:
- `function = my_function` — use a function instead of `value` for spatially or
  temporally varying source
- A custom kernel for a source that depends on the solution `T` (e.g., a nonlinear
  reaction term)

---

### Block: `[BCs]`

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

A single `DirichletBC` sub-block applies to all four walls simultaneously by listing all
four boundary names as a space-separated string in single quotes:
`boundary = 'left right top bottom'`.

This is more concise than writing four separate sub-blocks. MOOSE applies T=0 at all
nodes on all four named boundaries — a total of (21+21+19+19) = 80 boundary nodes
(corners are shared, so we subtract the 4 corners from the double-counting).

**Physical meaning**: the walls act as perfect heat sinks at 0 degrees. All the heat
generated by `Q=1` inside the domain must eventually flow out through these walls.
At steady state, the heat generation rate exactly balances the heat loss rate.

---

### Block: `[Materials]`

```
[Materials]
  [thermal_props]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '1.0'
  []
[]
```

The `[Materials]` block is new in Case 03. It defines material properties that kernels
and other MOOSE objects can read at quadrature points.

**`type = GenericConstantMaterial`**

The simplest material type: declares one or more named properties with constant values.
`prop_names` and `prop_values` are paired lists (same length).

**`prop_names = 'k'`** — declares a material property named `k`.

The name `k` is a string that identifies this property in the material system. Any
kernel or object that requests a property named `k` from the material system will get
the value declared here.

**`prop_values = '1.0'`** — the constant value of `k`.

With `k=1.0`, the heat equation becomes `dT/dt = div(grad T) + Q`, with unit thermal
conductivity. In a real problem you might set `k = 50.0` (steel, in SI units) or
`k = 0.6` (water).

**Why is this better than a hardcoded coefficient?**

Consider what happens when you want to model a plate made of two different materials
(Case 06 does exactly this). You can replace `GenericConstantMaterial` with a material
object that returns different values of `k` depending on which region of the mesh each
quadrature point is in — without changing the `MatDiffusion` kernel at all.

Or, for temperature-dependent conductivity: replace with a material that computes
`k = k0 * (1 + alpha * T)` — still without touching the kernel.

**`rho` and `cp` are absent — why?**

The input file comments explain: `rho=cp=1`. With unit values, `rho * cp * dT/dt`
simplifies to `dT/dt`. MOOSE's `TimeDerivative` kernel already implements `dT/dt`
(unit coefficient). To have a non-unit `rho*cp`, you would use `ADMatTimeDerivative`
or a similar kernel that reads a material property for the coefficient.

---

### Block: `[Postprocessors]`

```
[Postprocessors]
  [avg_temperature]
    type     = ElementAverageValue
    variable = T
  []

  [max_temperature]
    type          = ElementExtremeValue
    variable      = T
    value_type    = max
  []
[]
```

Postprocessors are scalar quantities computed from the solution at the end of each
timestep. They compress the full spatial field into a single number, making it easy
to track how the solution evolves over time.

**`[avg_temperature]` — ElementAverageValue**

Computes the spatial average of `T` over the entire domain:

```
avg_T = (1 / |Omega|) * integral( T  dV )
```

where `|Omega|` is the domain area (1.0 for the unit square).

This integral is evaluated using Gauss quadrature over all 400 elements. The result
is a single number printed to the console at each timestep and written to the CSV file.

At t=0: `avg_T = 0.0` (plate starts cold).
At t → infinity: `avg_T → 0.0127` approximately (the spatial average of the steady-state
solution for Q=1 on the unit square with zero BCs).

**`[max_temperature]` — ElementExtremeValue**

```
[max_temperature]
  type          = ElementExtremeValue
  variable      = T
  value_type    = max
[]
```

Finds the maximum value of `T` over all quadrature points in all elements. This tracks
the peak temperature in the plate (which always occurs near the center, due to symmetry).

`value_type = max` — alternatively `value_type = min` to track the minimum.

At t=0: `max_T = 0.0`.
At t → infinity: `max_T → 0.0737` approximately (the peak of the steady-state solution
at the domain center).

**What are postprocessors used for in practice?**

- Monitoring convergence to steady state (plot avg_T vs time; when it flattens, you
  have reached steady state)
- Safety criteria (is the peak temperature below a material limit?)
- Integral quantities (total energy in the domain)
- Automated verification (write expected values to CSV, compare with known results)

---

### Block: `[Executioner]`

```
[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type = ConstantDT
    dt   = 0.01
  []

  start_time = 0.0
  end_time   = 0.5

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
```

This is the most complex `[Executioner]` block so far. Several new parameters appear.

**`type = Transient`**

Switches from a one-shot steady solve to a time-marching loop. MOOSE repeatedly:
1. Advances time by `dt`
2. Assembles the residual and Jacobian for the new time level
3. Solves the nonlinear system `F(T^{n+1}) = 0`
4. Computes postprocessors and writes output
5. Repeats until `end_time` is reached

**`solve_type = 'PJFNK'`**

Same solver as before. For a transient, each time step involves one Newton solve.
Because the heat equation with a constant source is linear in T (the PDE contains no
`T^2` or `T * dT/dx` terms), Newton still converges in one iteration per time step.

**`[TimeStepper]` sub-block**

```
[TimeStepper]
  type = ConstantDT
  dt   = 0.01
[]
```

The `TimeStepper` controls how large each time step is.

- `type = ConstantDT` — fixed time step size. Simple and predictable.
- `dt = 0.01` — each time step advances t by 0.01 seconds.

With `start_time = 0.0` and `end_time = 0.5`, there will be exactly:

```
(0.5 - 0.0) / 0.01 = 50 time steps
```

**Why 0.01?** The time step should be small enough to resolve the transient dynamics.
The characteristic time scale `tau ≈ 0.1` seconds means a step of 0.01 (1/10th of tau)
gives good temporal resolution. Larger steps would miss details of the early transient.
Smaller steps would be accurate but slower.

For this problem (backward Euler, linear PDE), any `dt > 0` gives a stable solution
(backward Euler is unconditionally stable). The accuracy is O(dt) — halving `dt` halves
the time-integration error.

Alternative `TimeStepper` types include `IterationAdaptiveDT` (adjusts `dt` based on
Newton iteration count, introduced in Case 11) and `FunctionDT` (prescribed time
step schedule).

**`start_time = 0.0`**

The simulation begins at t=0. All initial conditions are applied at this time.
MOOSE uses `start_time` when evaluating time-dependent functions: a function specified
as `expression = 't'` would return 0.0 at the first step.

**`end_time = 0.5`**

The simulation stops after reaching t=0.5 seconds. The solution is output at every
time step (every 0.01 seconds), so the CSV file will have 51 rows (t=0 through t=0.5).

**`nl_rel_tol = 1e-8`**

The relative Newton convergence tolerance. The solver declares convergence when:

```
|R^{k+1}| / |R^0| < nl_rel_tol = 1e-8
```

where `R^k` is the residual at Newton iteration k and `R^0` is the initial residual.
This means the residual must be reduced by a factor of 10^8 from its starting value.

**`nl_abs_tol = 1e-10`**

The absolute Newton convergence tolerance. The solver also declares convergence when:

```
|R^{k+1}| < nl_abs_tol = 1e-10
```

regardless of the initial residual. This prevents the solver from stalling on very
small problems where the relative tolerance is trivially met immediately.

MOOSE declares convergence when **either** tolerance is satisfied. Setting both is
robust: `nl_rel_tol` handles typical cases, `nl_abs_tol` handles near-zero residuals.

Why are these tighter than the defaults? For a transient problem with 50 timesteps,
small errors accumulate. Tighter tolerances at each step reduce temporal drift. For
a linear problem like this, tighter tolerances cost nothing extra (Newton still
converges in 1 iteration to machine precision regardless).

---

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

Both output types are enabled.

**`exodus = true`**

For a transient problem, the Exodus file stores the solution field at **every time step**.
The file `case03_heat_transient_out.e` contains 51 time "frames" (t=0 through t=0.5
in 0.01 increments). In ParaView, you can animate through these frames to watch the
temperature field evolve.

**`csv = true`**

The CSV file `case03_heat_transient_out.csv` records postprocessor values at every
time step. With two postprocessors (`avg_temperature` and `max_temperature`), the CSV
has three columns and 51 rows:

```
time,avg_temperature,max_temperature
0,0,0
0.01,<value>,<value>
0.02,<value>,<value>
...
0.5,<value>,<value>
```

This time-history CSV is the primary artifact for analyzing transient behavior: plot
`avg_temperature` vs `time` to see the rise toward steady state, confirm that the
asymptotic value matches the theoretical prediction.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case03
../../test/moose_test-opt -i case03_heat_transient.i
```

### Step-by-Step Internally

**1. Startup**

MOOSE parses the input, creates all objects, and generates the mesh (400 elements,
441 nodes). The DofMap assigns 441 DOFs. Variable T is initialized to 0 everywhere.

**2. Initial Output**

At t=0, before any stepping, MOOSE computes postprocessors (both are 0.0) and writes
the initial state to the Exodus file (frame 0) and CSV (row 1: `0,0,0`).

**3. Time Step Loop (50 iterations)**

For each time step from 1 to 50:

  a. **Advance time**: t_new = t_old + dt = t_old + 0.01

  b. **Assemble residual**: Loop over 400 elements, calling three kernels at each:
     - `TimeDerivative`: computes `(T_new - T_old) / dt` at each quadrature point,
       integrates against test functions
     - `MatDiffusion`: computes `k * grad(T_new) . grad(phi_i)`, reads `k=1` from
       the material system
     - `BodyForce`: computes `1.0 * phi_i` (constant source)

  c. **Apply BCs**: Replace DOF equations for the 80 boundary nodes with T=0

  d. **Solve**: GMRES + BoomerAMG solves the linear system. For this linear problem,
     one Newton step drives the residual to machine precision.

  e. **Update solution**: T^{n+1} becomes the new T^n for the next step

  f. **Compute postprocessors**: Evaluate `avg_temperature` and `max_temperature`
     over the updated solution field

  g. **Write output**: Append a new time frame to the Exodus file, write a new row
     to the CSV file

**4. Termination**

After step 50 (t=0.5), the time loop ends. MOOSE reports timing statistics and exits.

### Typical Console Output

For a transient solve, MOOSE prints a summary at each time step:

```
Time Step 1, time = 0.01
 0 Nonlinear |R| = 2.500000e-01
      0 Linear |R| = 2.500000e-01
      1 Linear |R| = 3.721869e-17
 1 Nonlinear |R| = 2.775558e-17

Time Step 2, time = 0.02
 0 Nonlinear |R| = 2.467354e-01
      ...

...

Time Step 50, time = 0.5
 0 Nonlinear |R| = 4.312743e-05
      ...
 1 Nonlinear |R| = 1.234567e-17

Postprocessor Values:
  avg_temperature = 1.2682e-02
  max_temperature = 7.3707e-02
```

What the output tells you:

- `Time Step N, time = T` — the current time step number and time
- The initial nonlinear residual (before Newton) decreases over time: the solution
  changes less between steps as it approaches steady state, so the initial residual
  (which measures how far the old solution is from satisfying the new-time equations)
  gets smaller
- `1 Nonlinear |R| = ...` — after one Newton step, convergence to machine precision
- Postprocessor values at each step are printed (if `print_postprocessors` is true,
  which is the default)

---

## Output Files

### `case03_heat_transient_out.e`

The Exodus file stores **all 51 time frames** of the temperature field. File size will
be larger than the steady-state cases (roughly 50x larger for the same mesh).

**Visualizing in ParaView:**

1. File > Open, select `case03_heat_transient_out.e`, click OK, click Apply
2. The solution at t=0 (all zeros) appears — the plate is uniformly black/blue
3. Use the VCR controls (play/pause/step) in the toolbar to animate through time steps
4. Watch the temperature rise from 0, forming a dome shape centered at (0.5, 0.5)
5. By t=0.3-0.4, the shape barely changes between steps — you have reached near-steady-state

**Creating an animation:**

In ParaView:
1. Set the view angle and color map as desired
2. File > Save Animation
3. Choose number of frames (50), output format (AVI or PNG sequence)
4. The exported animation shows the transient temperature rise

**What to look for in the field:**

- At early times (t = 0.01-0.05): temperature barely above zero, concentrated in the
  very center, small and round
- At intermediate times (t = 0.1-0.3): temperature dome growing, spreading toward walls
- At late times (t = 0.4-0.5): dome reaches nearly its final shape — a smooth bell
  peaked at the center, tending to zero at all walls

### `case03_heat_transient_out.csv`

The CSV file records scalar postprocessor history. It has this structure:

```
time,avg_temperature,max_temperature
0,0,0
0.01,<avg>,<max>
0.02,<avg>,<max>
...
0.5,<avg>,<max>
```

**Column meanings:**

| Column | Meaning | Expected range |
|--------|---------|----------------|
| `time` | Simulation time (seconds) | 0.0 to 0.5 |
| `avg_temperature` | Spatial average of T over the domain | 0 to ~0.0127 |
| `max_temperature` | Maximum T over all quadrature points | 0 to ~0.0737 |

**Plotting the CSV:**

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('case03_heat_transient_out.csv')
plt.figure(figsize=(10, 5))
plt.plot(df['time'], df['avg_temperature'], label='Average T')
plt.plot(df['time'], df['max_temperature'], label='Max T')
plt.xlabel('Time (s)')
plt.ylabel('Temperature')
plt.title('Case 03: Transient Heat Equation')
plt.legend()
plt.grid(True)
plt.show()
```

The plot shows two curves that rise quickly at first, then flatten as they approach
their steady-state limits.

---

## Interpreting the Results

### What the Solution Looks Like

**Spatial field at t=0.5 (near steady state):**
- A smooth 2D bell curve centered at (0.5, 0.5)
- Maximum temperature approximately 0.0737 at the center
- Perfect 4-fold symmetry (the problem is symmetric about x=0.5 and y=0.5)
- Decreasing monotonically toward zero at all four walls
- Contour lines are roughly elliptical (actually they are rounded squares that become
  more rectangular far from the center)

**Time history of avg_temperature:**
- t=0: starts at exactly 0
- t=0.01-0.05: rapid initial rise (the plate is cold and heater is running full blast)
- t=0.1-0.3: rise slows as the plate warms and heat loss through walls increases
- t=0.4-0.5: very slow rise, almost flat — approaching the steady state of ~0.0127

**Time history of max_temperature:**
- Same qualitative behavior, but higher values
- Asymptotic value approximately 0.0737 at t → infinity

### Verifying Correctness

**Check 1: Steady-state limit**

At t=0.5, the simulation should be near (but not exactly at) steady state. The
average temperature should be close to the theoretical value:

```
T_avg_ss ≈ Q * A / (8 * pi^2 * k) ≈ 1.0 / (8 * pi^2) ≈ 0.01267
```

If your CSV shows `avg_temperature` at t=0.5 close to 0.0127, the simulation is correct.

**Check 2: Energy balance**

At every timestep, compute:
- Heat in: `Q * domain_area * dt = 1 * 1 * 0.01 = 0.01` (joules per step)
- Heat stored: `rho * cp * domain_area * d(avg_T)/dt * dt = 1 * delta_avg_T`
- Heat out: heat_in - heat_stored (via the walls)

At early times, nearly all heat is stored. At late times (near steady state), nearly
all heat goes out the walls. This balance is enforced by the PDE and should be
approximately satisfied by the numerical solution.

**Check 3: Symmetry**

The solution should be symmetric about x=0.5 and y=0.5. Check in ParaView:
- `u(0.3, 0.5) ≈ u(0.7, 0.5)` (left-right symmetry)
- `u(0.5, 0.3) ≈ u(0.5, 0.7)` (top-bottom symmetry)
- `u(0.3, 0.3) ≈ u(0.7, 0.3) ≈ u(0.3, 0.7) ≈ u(0.7, 0.7)` (fourfold symmetry)

### Physical Insight

Several physical principles are visible in this simulation:

1. **Thermal inertia**: the plate takes time to heat up because `rho*cp > 0`. An infinite
   thermal mass (`rho*cp → infinity`) would never heat up; zero thermal mass
   (`rho*cp → 0`) would reach steady state instantaneously.

2. **Diffusion smooths**: at t=0 the source is uniform, but the temperature distribution
   is not uniform — it is a smooth dome. Diffusion spreads heat from the center (where
   it accumulates) toward the cold walls. This smoothing is the hallmark of diffusion
   operators.

3. **Boundary effects**: the temperature near the walls (x≈0, x≈1, y≈0, y≈1) rises
   more slowly than at the center, even though the source `Q=1` is uniform. The cold
   walls drain heat from their neighborhood, keeping those regions cooler.

4. **Asymptotic approach**: the temperature never actually reaches steady state in finite
   time. It approaches exponentially: `|T(t) - T_ss| ~ exp(-t/tau)`. After `5*tau`,
   the error is `exp(-5) ≈ 0.7%` of the initial residual. This is why the simulation
   at t=0.5 (= 5*tau) is very close to but not exactly at steady state.

---

## Key Concepts Learned

- **Transient executioner**: `type = Transient` adds a time-marching outer loop;
  `start_time`, `end_time`, and `dt` control the time axis
- **TimeDerivative kernel**: implements `integral(phi_i * dT/dt dV)` in the weak form;
  works with backward Euler time integration (unconditionally stable, first-order accurate)
- **MatDiffusion kernel**: reads thermal conductivity `k` from the material system;
  more general than `Diffusion` because `k` can be spatially varying or nonlinear
- **Materials block**: `GenericConstantMaterial` declares named scalar material
  properties with constant values; kernels look up properties by name at quadrature points
- **BodyForce kernel**: adds a volumetric source term to the right-hand side; `value = 1.0`
  sets a constant uniform source
- **Postprocessors**: scalar summaries computed from the field solution at each timestep;
  appear in console output and CSV files; `ElementAverageValue` computes spatial averages,
  `ElementExtremeValue` tracks maxima or minima
- **CSV output**: time-history of scalar postprocessor values, one row per timestep;
  ideal for plotting with Python/matplotlib or any spreadsheet tool
- **ConstantDT TimeStepper**: fixed time step size; alternative to adaptive stepping
- **Newton convergence tolerances**: `nl_rel_tol` (relative) and `nl_abs_tol` (absolute)
  control when the Newton solver declares convergence at each timestep
- **Approach to steady state**: transient solutions approach steady state exponentially;
  the time constant is determined by the diffusivity and domain size
- **Exodus time frames**: for transient problems, the Exodus file stores all timesteps
  as frames that can be animated in ParaView

---

## Experiments to Try

### Experiment 1: Longer Simulation to Reach Steady State

Change `end_time` from 0.5 to 2.0 (and possibly increase `dt` to 0.05 to keep the
run fast):

```
[TimeStepper]
  type = ConstantDT
  dt   = 0.05
[]
start_time = 0.0
end_time   = 2.0
```

Run with 40 timesteps. The postprocessor values at t=2.0 should be very close to the
theoretical steady-state values (avg ≈ 0.01267, max ≈ 0.0737). Plot the CSV to see the
full S-curve approach to steady state.

Expected outcome: `avg_temperature` at t=2.0 matches the theoretical value to within
a fraction of a percent.

### Experiment 2: Change the Source Strength

Multiply or divide the source term:

```
# Stronger source:
[heat_source]
  type  = BodyForce
  variable = T
  value = 10.0
[]

# Weaker source:
[heat_source]
  type  = BodyForce
  variable = T
  value = 0.1
[]
```

Since the heat equation is linear, all temperatures scale linearly with Q. With Q=10,
the steady-state average temperature should be 10x larger (≈ 0.127). With Q=0.1,
it should be 10x smaller (≈ 0.00127). The time to reach steady state (determined by
`tau = rho*cp*L^2/k`) is independent of Q.

### Experiment 3: Change the Thermal Conductivity

Increase or decrease `k` in the `[Materials]` block:

```
# Higher conductivity (better conductor):
prop_values = '10.0'

# Lower conductivity (poorer conductor):
prop_values = '0.1'
```

Higher `k` means:
- The steady-state temperature is **lower** (heat escapes more easily)
- The time constant is **shorter** (approaches steady state faster)

Lower `k` means:
- The steady-state temperature is **higher** (heat stays trapped longer)
- The time constant is **longer** (approaches steady state more slowly)

The steady-state average temperature scales as `1/k` (confirmed by the theory).

### Experiment 4: Vary the Time Step Size

Try different fixed timestep sizes:

```
dt = 0.1    # 5 time steps to reach t=0.5 (very coarse in time)
dt = 0.001  # 500 time steps to reach t=0.5 (very fine in time)
```

Compare the `avg_temperature` at t=0.5 for all three choices (0.1, 0.01, 0.001).
Because backward Euler is first-order in time, halving `dt` should approximately halve
the time-integration error (the difference between the computed `avg_T` and the true
value). This is a temporal convergence study.

Note that the coarse `dt=0.1` gives only 5 CSV data points — the time history is
poorly resolved. The fine `dt=0.001` gives 500 points — very smooth but slower to run.
This trade-off is why adaptive time stepping (Case 11) is useful.

### Experiment 5: Non-Uniform Source

Replace the constant source with a spatially varying one using a `BodyForce` with a
function:

```
[Functions]
  [source_fn]
    type       = ParsedFunction
    expression = 'sin(pi*x) * sin(pi*y)'
  []
[]

[Kernels]
  ...
  [heat_source]
    type     = BodyForce
    variable = T
    function = source_fn    # use 'function' instead of 'value'
  []
[]
```

The source is now concentrated in the center of the plate (the sine function peaks at
(0.5, 0.5) and is zero on the walls). The steady-state temperature distribution will
be even more concentrated in the center than with the uniform source. Compare the
steady-state `max_temperature` values: they should be higher for the concentrated
source (same total energy, but deposited more centrally).
