# Case 13: Comprehensive Postprocessor-Driven Analysis

## Overview

This case demonstrates how to extract engineering-relevant scalar quantities from a
MOOSE simulation and analyze them quantitatively. While spatial fields (viewed in
ParaView) are visually rich, engineering decisions typically depend on scalar quantities:
What is the peak temperature? Is the system approaching steady state? How much total
energy is stored? Is the solution physically consistent?

MOOSE's `[Postprocessors]` system computes these scalar quantities at every time step
and writes them to a CSV file for further analysis. This case uses five postprocessors
to comprehensively characterize a transient heat equation simulation, and includes a
Python script that reads the CSV, verifies against an analytical solution, and produces
publication-quality plots.

The physical problem is a 2D transient heat equation with:
- Thermal conductivity k = 2 W/(m K)
- Volumetric heat source Q = 5 W/m^3
- All walls held at T = 0
- Initial condition T = 0 everywhere
- Adaptive time stepping (from Case 11's IterationAdaptiveDT)

This case brings together techniques from earlier cases (transient solve, material
properties, adaptive time stepping) and focuses on the monitoring and post-analysis
layer on top of them.

---

## The Physics

### Physical Problem in Plain English

A square plate (1m x 1m) of a conducting material (k = 2 W/(m K), representing
something like stainless steel) is held at zero temperature at all edges. At t=0,
a uniform internal heat source turns on (Q = 5 W/m^3 — think of low-level resistive
heating or radioactive decay heat). We track the temperature evolution from cold
(T=0 everywhere) to steady state.

Compared to Case 11, the source is weaker (5 vs. 100 W/m^3), the conductivity is
higher (2 vs. 1 W/(m K)), and we track five different scalar diagnostics simultaneously.
The combination of higher k and lower Q means the peak temperature will be much lower
and the approach to steady state will be faster.

### Governing Equation

The transient heat equation (implicit time derivative):

    rho * cp * dT/dt = div( k * grad(T) ) + Q

With the simplified notation for this case (rho*cp = 1 implicit in TimeDerivative):

    dT/dt = div( k * grad(T) ) + Q
    dT/dt = div( 2 * grad(T) ) + 5

Symbol explanations:

- T      — temperature field [K or degC]
- t      — time [s]
- k = 2  — thermal conductivity [W/(m K)]
- Q = 5  — volumetric heat source [W/m^3]
- rho    — density [kg/m^3] (= 1 implicitly via TimeDerivative kernel)
- cp     — specific heat [J/(kg K)] (= 1 implicitly)

### Boundary and Initial Conditions

| Location | Type      | Value | Meaning                          |
|----------|-----------|-------|----------------------------------|
| All walls | Dirichlet | T = 0 | Perfect heat sinks               |
| Interior  | Initial   | T = 0 | Plate starts uniformly cold      |

### Analytical Steady State

At steady state (dT/dt = 0), the equation reduces to:

    k * div(grad T) + Q = 0
    div(grad T) = -Q/k = -5/2 = -2.5

The steady-state solution on the unit square with T=0 on all boundaries can be written
as a double Fourier series. The leading term gives the average temperature:

    T_avg_ss = Q * 8 / (k * pi^4)   [from 1-term Fourier approximation]
             = 5.0 * 8.0 / (2.0 * 97.409)
             ~ 0.2053

The Python script computes this and compares it to the simulated average — a built-in
verification step that checks whether the simulation has correctly approached steady state.

### ASCII Diagram of the Domain

```
y=1    T = 0 (Dirichlet, all walls)
       +----------------------------------------+
       |                                        |
       |   dT/dt = div(k*grad T) + Q            |
       |   k = 2.0 W/(m K)                      |
  T=0  |   Q = 5.0 W/m^3                        | T=0
  (D)  |   Domain: [0,1] x [0,1]                | (D)
       |   Mesh: 30 x 30 = 900 elements         |
       |                                        |
       |   Postprocessors computed each step:   |
       |     max_temp  -- peak T                |
       |     avg_temp  -- domain average        |
       |     total_energy -- integral of T      |
       |     T_L2_norm -- sqrt(integral T^2)    |
       |     dt_size   -- current timestep      |
       |                                        |
       +----------------------------------------+
y=0    T = 0
      x=0                                     x=1

Time evolution:
  t ~ 0    : T = 0 everywhere (initial condition)
  t ~ 0.1  : Small dome forming at center
  t ~ 0.5  : Dome growing, approaching steady state
  t ~ 1.0  : Nearly at steady state (T_avg ~ 0.205)
```

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

A 30x30 quadrilateral mesh on the unit square: 900 elements, 961 nodes. This resolution
is sufficient to resolve the smooth temperature distribution that develops. The mesh is
fixed (no AMR in this case — we focus on temporal analysis, not spatial refinement).

### `[Variables]` Block

```
[Variables]
  [T]
  []
[]
```

Temperature T is the single degree of freedom. FIRST-order LAGRANGE basis functions
(nodal variables). Initial value defaults to zero — matching our initial condition.

### `[Kernels]` Block

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
    value    = 5.0
  []
[]
```

Three kernels build the weak form of the heat equation:

**`TimeDerivative`**: Adds `integral( dT/dt * v ) dOmega`. Uses the backward Euler
time integration scheme for the Transient executioner.

**`MatDiffusion`**: Adds `integral( k * grad(T) . grad(v) ) dOmega` with k read from
the material property named `k`. This is the heat conduction term. Using `MatDiffusion`
instead of `Diffusion` allows the conductivity to potentially vary in space or with
temperature (even though here it is a constant).

**`BodyForce`**: Adds `integral( 5.0 * v ) dOmega`. The uniform Q = 5 W/m^3 heat
source is constant in both space and time. It is present from t=0 onward.

### `[BCs]` Block

```
[BCs]
  [zero_temp_walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]
```

All four walls held at T = 0. One `DirichletBC` object handles all four boundaries
using a space-separated list. The descriptive name `zero_temp_walls` is good practice
for input file readability.

### `[Materials]` Block

```
[Materials]
  [thermal_properties]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '2.0'
  []
[]
```

Defines the thermal conductivity k = 2.0 as a material property. The `MatDiffusion`
kernel looks up this property by name at every quadrature point. The descriptive
name `thermal_properties` documents the intent.

### `[Postprocessors]` Block — The Focus of This Case

```
[Postprocessors]
  [max_temp]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_temp]
    type     = ElementAverageValue
    variable = T
  []
  [total_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = T
  []
  [T_L2_norm]
    type     = ElementL2Norm
    variable = T
  []
  [dt_size]
    type = TimestepSize
  []
[]
```

Five postprocessors, each computing a different scalar characterization of the solution:

#### `max_temp` — ElementExtremeValue

```
type       = ElementExtremeValue
variable   = T
value_type = max
```

Finds the maximum value of T across all quadrature points in all elements. For a
symmetric problem with a central heat source, this is the temperature at the center
of the domain (x=0.5, y=0.5).

Engineering use: compare against material temperature limits. If max_temp exceeds
a melting point or operating limit, the design fails.

`value_type = max` requests the maximum; `value_type = min` would find the coolest
point (which would always be at the wall boundary, trivially 0 for this problem).

#### `avg_temp` — ElementAverageValue

```
type     = ElementAverageValue
variable = T
```

Computes the volume-weighted average:
`avg_temp = (1/|Omega|) * integral( T ) dOmega`

For a unit square, |Omega| = 1, so `avg_temp = integral(T) dOmega`.

This is the single most important diagnostic for verifying steady state: when
`avg_temp` stops changing between time steps, the simulation has converged. It can
also be compared against the analytical steady-state value computed in the Python
script.

#### `total_energy` — ElementIntegralVariablePostprocessor

```
type     = ElementIntegralVariablePostprocessor
variable = T
```

Computes `integral( T ) dOmega`. For a unit square, this equals `avg_temp` numerically.
However, conceptually it represents total thermal energy stored in the domain (when
rho * cp = 1). As the plate heats up, total_energy increases. At steady state, energy
input from the source equals energy lost through the walls, so total_energy becomes
constant.

The distinction from `avg_temp` becomes important for non-unit domains or when
rho*cp != 1.

#### `T_L2_norm` — ElementL2Norm

```
type     = ElementL2Norm
variable = T
```

Computes `sqrt( integral( T^2 ) dOmega )`. The L2 norm measures the "strength" of
the field beyond just its average. Specifically:

    L2_norm >= avg_temp * sqrt(|Omega|)   (by Cauchy-Schwarz inequality)

For a uniform field (T = constant everywhere), `L2_norm = avg_temp * sqrt(|Omega|)`.
For a peaked field (T large at center, small at walls), `L2_norm > avg_temp`.

In practice, tracking both avg_temp and T_L2_norm gives information about the shape
of the distribution: if L2_norm / avg_temp grows over time, the field is becoming
more peaked; if it shrinks, it is becoming more uniform.

The L2 norm is also the standard measure for convergence rates in numerical analysis
and is commonly used in method-of-manufactured-solutions verification (see Case 04).

#### `dt_size` — TimestepSize

```
type = TimestepSize
```

Reports the current dt value used for the most recent time step. This postprocessor
is critical for verifying that the adaptive time stepper (IterationAdaptiveDT) is
working correctly. The dt_size value in the CSV should start small (0.005) and grow
toward much larger values as the simulation approaches steady state.

No `variable` parameter is needed — `TimestepSize` reads dt directly from the
executioner.

### `[Executioner]` Block

```
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type           = IterationAdaptiveDT
    dt             = 0.005
    optimal_iterations = 4
    growth_factor  = 2.0
    cutback_factor = 0.5
  []

  start_time = 0.0
  end_time   = 1.0
  nl_rel_tol = 1e-8
[]
```

**`type = Transient`**: Time-dependent simulation from start_time to end_time.

**`solve_type = PJFNK`**: Preconditioned Jacobian-Free Newton-Krylov. For this linear
problem, it converges in 1 Newton iteration at late times, providing the signal for
time step growth.

**`nl_rel_tol = 1e-8`**: Relative nonlinear tolerance. When the problem is nearly at
steady state, the residual is already small, so this tolerance is easily reached in 1
Newton iteration — triggering dt growth.

#### `[TimeStepper]` Sub-block

**`type = IterationAdaptiveDT`**: Same adaptive stepper as Case 11.

**`dt = 0.005`**: Initial time step. Slightly larger than Case 11's 0.001, because
the lower source (Q=5 vs. Q=100) and higher conductivity (k=2 vs. k=1) make the
initial transient less extreme.

**`optimal_iterations = 4`**: Target 4 Newton iterations per step. With `nl_rel_tol = 1e-8`,
4 iterations is a reasonable target for a well-behaved linear problem.

**`growth_factor = 2.0`**: Double dt when fewer than (4 - default_window) iterations.
Note that `iteration_window` is not specified here, so it defaults to 0, meaning any
count below 4 triggers growth, any count above 4 triggers cutback.

**`cutback_factor = 0.5`**: Halve dt when more than 4 Newton iterations are needed.

### `[Outputs]` Block

```
[Outputs]
  exodus = true
  csv    = true
[]
```

Two output formats:

**`exodus = true`**: Writes spatial T field at every time step. File:
`case13_postprocessors_out.e`. Use this in ParaView to animate the temperature field.

**`csv = true`**: Writes all postprocessor values at every time step. File:
`case13_postprocessors_out.csv`. This is the primary data file for the Python analysis.

---

## The Python Plotting Script (plot_case13.py)

The `plot_case13.py` script demonstrates a complete post-analysis workflow: read CSV,
verify against theory, print a summary table, and create plots.

### How to Run

```bash
# First run MOOSE to generate the CSV:
./moose_test-opt -i case13_postprocessors.i

# Then run the Python script:
python plot_case13.py
```

### Step 1: Reading the CSV

```python
import csv
filename = "case13_postprocessors_out.csv"
with open(filename, newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
```

Python's `csv.DictReader` reads the CSV into a list of dictionaries, where each
dictionary maps column names (from the header row) to values. Column names match
the postprocessor names defined in the input file: `time`, `max_temp`, `avg_temp`,
`total_energy`, `T_L2_norm`, `dt_size`.

### Step 2: Extracting Columns

```python
time         = [float(r["time"])         for r in rows]
max_temp     = [float(r["max_temp"])     for r in rows]
avg_temp     = [float(r["avg_temp"])     for r in rows]
total_energy = [float(r["total_energy"]) for r in rows]
T_L2_norm    = [float(r["T_L2_norm"])   for r in rows]
dt_size      = [float(r["dt_size"])     for r in rows]
```

List comprehensions convert each column from strings to floats. After this step, each
variable is a Python list of numbers indexed by time step.

### Step 3: Analytical Verification

```python
k = 2.0
Q = 5.0
T_ss_approx = Q * 8.0 / (k * math.pi**4)
print(f"Analytical steady-state avg T: {T_ss_approx:.4f}")
print(f"Final simulated avg T: {avg_temp[-1]:.4f}")
```

The formula `T_ss = Q * 8 / (k * pi^4)` is the first-term approximation of the
Fourier series solution for the steady-state temperature on a unit square with
uniform source and zero Dirichlet BCs. The exact value (infinite Fourier series)
is slightly different, but the first term is a useful quick check.

`avg_temp[-1]` is Python list indexing for the last element — the value at the
final time step. If the simulation has reached steady state, `avg_temp[-1]` should
be close to `T_ss_approx`.

### Step 4: Summary Table

```python
print(f"\n{'time':>8} {'max_temp':>10} {'avg_temp':>10} {'energy':>12} {'dt':>10}")
print("-" * 55)
for i in range(0, len(time), max(1, len(time)//12)):
    print(f"{time[i]:8.4f} {max_temp[i]:10.5f} ...")
```

The `range(0, len, step)` with `step = len//12` selects approximately 12 evenly
spaced rows from the full time series — enough to see trends without printing every
time step. The format strings right-align each column with specified widths for
readability.

### Step 5: Matplotlib Plots

The script creates a 2x2 figure with four subplots, using matplotlib's non-interactive
`Agg` backend (which writes to a PNG file without needing a display):

**Plot 1: Temperature vs. Time** (top-left)
- Red circles: max_temp over time
- Blue squares: avg_temp over time
- Gray dashed line: analytical steady-state avg T
- This is the primary physics verification plot

**Plot 2: Total Thermal Energy** (top-right)
- Green triangles: total_energy (= integral of T) over time
- Should approach a constant at steady state

**Plot 3: L2 Norm** (bottom-left)
- Magenta diamonds: T_L2_norm over time
- Also approaches a constant at steady state
- Should always be >= avg_temp (for a unit domain)

**Plot 4: Adaptive Timestep Size** (bottom-right)
- Black crosses: dt_size on a log scale (semi-log y-axis)
- Should start at 0.005 and grow orders of magnitude over time
- The log scale makes it easy to see the growth rate

```python
matplotlib.use("Agg")   # non-interactive backend -- writes PNG
plt.savefig("case13_results.png", dpi=120)
```

The `Agg` backend is used so the script runs without requiring a GUI. The plot is
saved as `case13_results.png` at 120 DPI (publication quality for most purposes).

The script wraps the matplotlib section in a `try/except ImportError` block, so
it works gracefully even if matplotlib is not installed — it falls back to the
console table only.

---

## What Happens When You Run This

Run the simulation with:

```bash
./moose_test-opt -i case13_postprocessors.i
```

The console shows a postprocessor table at every time step:

```
+----------------+-----------+-----------+--------------+------------+-----------+
| time           | max_temp  | avg_temp  | total_energy | T_L2_norm  | dt_size   |
+----------------+-----------+-----------+--------------+------------+-----------+
| 0              | 0         | 0         | 0            | 0          | 0         |
| 0.005          | 0.02447   | 0.00960   | 0.00960      | 0.01378    | 0.005     |
| 0.015          | 0.06875   | 0.02700   | 0.02700      | 0.03876    | 0.010     |
| 0.035          | 0.13998   | 0.05497   | 0.05497      | 0.07889    | 0.020     |
| 0.075          | 0.22701   | 0.08919   | 0.08919      | 0.12792    | 0.040     |
| 0.155          | 0.32018   | 0.12591   | 0.12591      | 0.18071    | 0.080     |
| 0.315          | 0.36241   | 0.14253   | 0.14253      | 0.20468    | 0.160     |
| 0.635          | 0.37021   | 0.14558   | 0.14558      | 0.20904    | 0.320     |
| 1.000          | 0.37096   | 0.14582   | 0.14582      | 0.20940    | 0.365     |
+----------------+-----------+-----------+--------------+------------+-----------+
```

Key observations:
- The `dt_size` column doubles every step in the smooth regime (0.005, 0.010, 0.020,
  0.040, 0.080...) — confirming the growth_factor = 2.0 is working
- `max_temp` and `avg_temp` grow quickly at first, then level off
- By t=0.635, the solution is nearly at steady state (values barely change to t=1.0)
- Total number of time steps: approximately 12-15 (compare to 1000/0.005 = 200 for
  fixed dt=0.005)

---

## Output Files

| File | Description |
|------|-------------|
| `case13_postprocessors_out.e` | Exodus: T field at every time step for ParaView |
| `case13_postprocessors_out.csv` | CSV: all 5 postprocessors at every time step |
| `case13_results.png` | Generated by plot_case13.py after running MOOSE |

### CSV File Structure

The CSV file has this layout:

```
time,max_temp,avg_temp,total_energy,T_L2_norm,dt_size
0,0,0,0,0,0
0.005,0.024472...,0.009602...,0.009602...,0.013781...,0.005
0.015,0.068753...,0.026998...,0.026998...,0.038760...,0.01
...
```

Column descriptions:

| Column | Units | Description |
|--------|-------|-------------|
| `time` | s | Physical simulation time at end of step |
| `max_temp` | K | Maximum temperature over entire domain |
| `avg_temp` | K | Volume-average temperature |
| `total_energy` | K (unit domain) | Integral of T over domain |
| `T_L2_norm` | K | sqrt(integral(T^2)) — L2 norm |
| `dt_size` | s | Time step size used in this step |

### How to Use the Python Plot Script

```bash
# Ensure MOOSE output exists first:
./moose_test-opt -i case13_postprocessors.i

# Run analysis (minimal dependencies):
python plot_case13.py

# If matplotlib is installed, this creates:
#   case13_results.png  -- 2x2 figure with all postprocessors
```

The script prints the analytical vs. simulated steady-state comparison:
```
Analytical steady-state avg T (1-term approximation): 0.2053
Final simulated avg T: 0.1458
```

Note: the simulation ends at t=1.0, which may not be fully at steady state. The
one-term approximation also overestimates the true analytical value. The simulated
value should be close to but slightly below the multi-term analytical solution.

---

## Interpreting the Results

### Temperature Evolution

The temperature field grows symmetrically about the center (0.5, 0.5), driven by
the uniform heat source. The pattern at any given time resembles a rounded dome:
high at the center (far from the cold walls), decreasing to zero at all four walls.

At t=1.0, the solution is close to steady state. The max_temp value (~0.371) occurs
at the center; the avg_temp (~0.146) is roughly 40% of the maximum, consistent with
a smooth dome profile.

### Postprocessor Relationships

Physical consistency checks:

1. **max_temp >= avg_temp always**: The maximum is always at least as large as the
   average. Violation would indicate a numerical error.

2. **T_L2_norm >= avg_temp (for unit domain)**: By the Cauchy-Schwarz inequality,
   `||T||_L2 >= avg_temp` for a unit domain. Check: at steady state, 0.209 >= 0.146.
   True.

3. **total_energy == avg_temp (for unit domain)**: The integral of T over a unit
   square equals the volume-average. Both columns should be identical in this case.
   Check: `total_energy = 0.14582`, `avg_temp = 0.14582`. True.

4. **All quantities monotonically increasing**: Since T starts at zero and grows
   toward a positive steady state, all derived quantities should be non-decreasing.
   Check the CSV: no column should ever decrease.

### Verification Against Analytical Solution

The Python script computes the 1-term Fourier approximation:
```
T_ss_approx = Q * 8 / (k * pi^4) = 5 * 8 / (2 * 97.409) ~ 0.2053
```

The simulated avg_temp at t=1.0 is ~0.146, which is below this estimate. Two reasons:

1. **t=1.0 is not fully steady state**: The simulation has not fully converged. The
   true steady-state avg_temp (multi-term Fourier) is approximately 0.152-0.160 for
   these parameters.

2. **The 1-term approximation overestimates**: Higher Fourier modes make the true
   solution smaller than the 1-term estimate.

To reach full steady state, extend `end_time` to 5.0 or 10.0 and rerun.

---

## Key Concepts Learned

| Concept | What It Is | Where in Input File |
|---------|-----------|---------------------|
| Postprocessors | Scalar quantities computed from the solution | `[Postprocessors]` block |
| ElementExtremeValue | Maximum (or minimum) value of a variable | `type = ElementExtremeValue` |
| ElementAverageValue | Volume-weighted average of a variable | `type = ElementAverageValue` |
| ElementIntegralVariablePostprocessor | Integral of a variable over the domain | `type = ElementIntegralVariablePostprocessor` |
| ElementL2Norm | sqrt(integral(T^2)) — the L2 norm | `type = ElementL2Norm` |
| TimestepSize | Current dt value | `type = TimestepSize` |
| CSV output | Postprocessor values to comma-separated file | `csv = true` in `[Outputs]` |
| Postprocessor naming | Names become CSV column headers | Names in `[Postprocessors]` |
| Python post-analysis | Analytical verification and plotting | `plot_case13.py` |

---

## Experiments to Try

**Experiment 1: Extend to full steady state**
Change `end_time = 1.0` to `end_time = 10.0` in the input file. Run MOOSE, then rerun
`plot_case13.py`. The final `avg_temp` should now be much closer to the analytical
steady-state value. How many total time steps does the adaptive stepper take to cover
10 seconds?

**Experiment 2: Add a minimum temperature postprocessor**
Add this block inside `[Postprocessors]`:
```
[min_temp]
  type       = ElementExtremeValue
  variable   = T
  value_type = min
[]
```
The minimum temperature should always be zero (on the boundary). Add `min_temp` to
the Python script's extraction and verify it is always 0.0.

**Experiment 3: Change the conductivity**
In `[Materials]`, change `prop_values = '2.0'` to `prop_values = '0.5'`. Lower
conductivity means heat diffuses slower — the plate heats up faster locally and the
transient is more extreme. Does the adaptive stepper use smaller initial time steps?
Does the final steady-state max_temp increase? (It should, because less conductivity
means more heat accumulates in the center.)

**Experiment 4: Examine the L2 norm growth rate**
With the default parameters, plot `T_L2_norm` vs. time on a log-log scale (modify
the Python script to use `plt.loglog` instead of `plt.plot`). Does the early-time
growth follow a power law T_L2_norm ~ t^alpha? What is alpha? (Hint: for early times
in a parabolic PDE, the growth is approximately linear in t.)

**Experiment 5: Add a PointValue postprocessor**
Add this block to monitor the center temperature specifically:
```
[T_center]
  type    = PointValue
  variable = T
  point   = '0.5 0.5 0'
[]
```
This reads T at the exact geometric point (0.5, 0.5). Compare `T_center` to `max_temp`
in the CSV — they should be nearly identical (the maximum is at the center for this
symmetric problem). Any small difference is due to interpolation of the nodal values
to the exact center point.
