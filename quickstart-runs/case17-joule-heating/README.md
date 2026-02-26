# Case 17: Joule Heating — Electric Current Generates Heat

## Overview

This case introduces **electro-thermal coupling**: a flowing electric current
dissipates energy as heat, raising the temperature of the conductor. The two
phenomena are linked in one direction — the electric field determines the heat
source, but the temperature does not feed back into the electrical conductivity
(one-way coupling). This is the simplest physically meaningful multi-physics
problem involving electricity and heat.

Two new MOOSE features appear here:

1. **Reusing a kernel for a different physics**: `ADHeatConduction` solves
   both the heat equation and the Laplace equation for voltage by pointing it
   at different material properties.
2. **ElectromagneticHeatingMaterial + ADJouleHeatingSource**: A material object
   computes the volumetric power density Q = sigma * |grad(V)|^2, which a
   kernel then injects as a body force into the temperature equation.

---

## The Physics

### Ohm's Law and the Laplace Equation

When a voltage difference is applied across a conductor, free electrons drift in
response to the electric field E = -grad(V). The current density J (amperes per
square metre) follows Ohm's law:

```
J = sigma * E = -sigma * grad(V)
```

where sigma is the electrical conductivity [S/m]. Conservation of charge in a
steady conductor (no accumulation) requires that the divergence of J is zero:

```
div(J) = 0   =>   -div(sigma * grad(V)) = 0
```

With uniform sigma this is simply the Laplace equation for V:

```
Laplacian(V) = 0
```

The solution on a 2D rectangle with V=10 on the left and V=0 on the right is a
linear ramp: V(x) = 10 * (1 - x/L), independent of y. The electric field
E = -dV/dx = 10/L points uniformly from left to right, and the current density
J = sigma * 10/L is uniform throughout the domain.

### Joule Dissipation

Every element of the conductor resists the current flow, converting electrical
energy into heat. The volumetric power density (Joule heating) is:

```
Q = J . E = sigma * |grad(V)|^2     [W/m^3]
```

In the uniform-field case considered here, Q is constant throughout the
domain:

```
Q = sigma * (V_left - V_right)^2 / L^2
  = 1e6 * (10)^2 / (2)^2
  = 25,000,000 W/m^3  (25 MW/m^3)
```

This large power density is intentional — it makes the temperature rise fast
enough to observe in a 5-second transient.

### Heat Equation with Joule Source

The temperature T satisfies the transient heat equation with Q as a body force:

```
rho * cp * dT/dt = div(k * grad(T)) + Q
```

Symbol definitions:

- rho = 8000 kg/m^3 — density (typical steel)
- cp  = 500 J/(kg K) — specific heat
- k   = 50 W/(m K)  — thermal conductivity
- Q   = sigma * |grad(V)|^2 — Joule heat source

Boundary conditions:

- T = 300 K on left and right boundaries: the electrodes act as perfect heat
  sinks, clamping the temperature at the electrode faces.
- Zero flux (insulated) on top and bottom: no heat leaves through the sides.

Initial condition: T = 300 K everywhere (uniform, equal to the electrode
temperature).

### Expected Analytical Behaviour

Because Q is uniform and the boundary conditions for T are symmetric (T=300 K
at both x=0 and x=L), the temperature profile at any time t is symmetric
about the midplane x=L/2. The steady-state temperature profile is parabolic:

```
T_ss(x) = 300 + Q * x * (L - x) / (2 * k)
```

At the midplane x = L/2 = 1 m this gives:

```
T_ss(L/2) = 300 + Q * L^2 / (8 * k)
           = 300 + 25e6 * 4 / (8 * 50)
           = 300 + 25,000 K  (very large!)
```

In practice the simulation is stopped at t = 5 s, long before steady state is
reached. During transient heating the peak temperature at the midplane grows
roughly as:

```
T_peak(t) ~ 300 + Q * t / (rho * cp)
           = 300 + 25e6 * t / (8000 * 500)
           = 300 + 6.25 * t   [K]
```

At t = 5 s: T_peak ≈ 300 + 31.25 = 331 K. The actual MOOSE result will be
slightly lower because heat conduction carries some energy to the cooled
electrodes, reducing the peak.

### ASCII Domain Diagram

```
y=1  T: insulated (zero flux)      V: insulated (zero flux)
     +------------------------------------------+
     |                                          |
     |   -div(sigma*grad(V)) = 0               |
     |   rho*cp*dT/dt = div(k*grad(T)) + Q    |
     |                                          |
V=10 |   sigma=1e6 S/m, k=50 W/mK             | V=0
T=300|   rho=8000 kg/m3, cp=500 J/kgK         | T=300
(Dir)|                                          | (Dir)
     |   Q = sigma*|grad(V)|^2 = 25 MW/m^3    |
     |   uniform throughout domain              |
     |                                          |
     +------------------------------------------+
y=0  T: insulated (zero flux)      V: insulated (zero flux)
     x=0                                     x=2

Mesh: 40 x 20 uniform quadrilateral elements
      800 elements, 861 nodes
      2 x 861 = 1722 total DOFs (V and T)
```

---

## Input File Walkthrough

The input file is `case17_joule_heating.i`.

### Block: `[Mesh]`

```
[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    nx    = 40
    ny    = 20
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1
  []
[]
```

A 2:1 aspect ratio rectangle (2 m x 1 m). With 40 x 20 elements the element
size is 0.05 m in both directions. The coarser y resolution (20 vs 40) is
acceptable because the steady-state solution for this case varies only in x.

### Block: `[Variables]`

```
[Variables]
  [V]
    initial_condition = 0.0
  []
  [T]
    initial_condition = 300.0
  []
[]
```

Two first-order Lagrange unknowns. V is given a zero initial condition (the
Laplace equation will immediately place the correct linear profile at the first
Newton iteration). T starts at 300 K to match the electrode boundary conditions.

### Block: `[Kernels]`

Four kernels implement the two PDEs.

**For V (electric potential):**

```
[V_diff]
  type                 = ADHeatConduction
  variable             = V
  thermal_conductivity = electrical_conductivity
[]
```

`ADHeatConduction` implements -div(k * grad(u)), reading the conductivity from
a material property. By setting `thermal_conductivity = electrical_conductivity`
we redirect it to read the electrical conductivity instead. This reuse of a heat
kernel for an electrostatic equation works because both PDEs have the same
mathematical form (a Laplacian with a material-property coefficient).

There is no time derivative for V because the Laplace equation is steady — the
electric field responds instantaneously to changes in boundary conditions.

**For T (temperature):**

```
[T_time]
  type     = ADHeatConductionTimeDerivative
  variable = T
[]
[T_diff]
  type     = ADHeatConduction
  variable = T
[]
[T_joule]
  type         = ADJouleHeatingSource
  variable     = T
  heating_term = electric_field_heating
[]
```

- `ADHeatConductionTimeDerivative` contributes rho * cp * dT/dt. It reads
  the material properties named `density` and `specific_heat` automatically.
- `ADHeatConduction` contributes -div(k * grad(T)). It reads `thermal_conductivity`.
- `ADJouleHeatingSource` adds the Joule body force Q. It reads the material
  property `electric_field_heating` computed by `ElectromagneticHeatingMaterial`.

### Block: `[BCs]`

```
[V_left]  type = ADDirichletBC  variable = V  boundary = left   value = 10.0
[V_right] type = ADDirichletBC  variable = V  boundary = right  value = 0.0
[T_left]  type = ADDirichletBC  variable = T  boundary = left   value = 300.0
[T_right] type = ADDirichletBC  variable = T  boundary = right  value = 300.0
```

`ADDirichletBC` is the AD version of `DirichletBC`. Using AD-consistent objects
throughout avoids mixed AD/non-AD assembly issues. Top and bottom boundaries
have no entries — MOOSE applies the natural (zero-flux Neumann) condition by
default.

### Block: `[Materials]`

```
[thermal]
  type        = ADGenericConstantMaterial
  prop_names  = 'thermal_conductivity specific_heat density'
  prop_values = '50.0              500.0        8000.0'
[]
[electrical]
  type        = ADGenericConstantMaterial
  prop_names  = 'electrical_conductivity'
  prop_values = '1e6'
[]
[joule_material]
  type                       = ElectromagneticHeatingMaterial
  electric_field             = V
  electric_field_heating_name = electric_field_heating
  electrical_conductivity    = electrical_conductivity
  formulation                = time
  solver                     = electrostatic
[]
```

`ElectromagneticHeatingMaterial` is the central coupling object. It reads the
gradient of V at every quadrature point, computes sigma * |grad(V)|^2, and
stores the result in the material property `electric_field_heating`. The
`ADJouleHeatingSource` kernel then reads this property and adds it to the T
residual. The parameters `formulation = time` and `solver = electrostatic`
select the correct physics variant (steady-state electrostatics, time domain).

### Block: `[Preconditioning]`

```
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]
```

`SMP` (Single Matrix Preconditioner) with `full = true` assembles the complete
Jacobian including all off-diagonal V-T coupling blocks. This is important
because the Joule source depends on grad(V), so the T residual has a derivative
with respect to V that must be captured for fast Newton convergence.

### Block: `[Executioner]`

```
[Executioner]
  type       = Transient
  solve_type = NEWTON
  dt         = 0.25
  end_time   = 5.0
[]
```

`NEWTON` (vs. `PJFNK`) uses the exact assembled Jacobian rather than a
finite-difference approximation. With AD objects throughout, the Jacobian is
computed exactly and cheaply via dual-number arithmetic. Fixed time stepping
with dt = 0.25 s gives 20 steps to reach t = 5 s.

---

## Running the Simulation

```bash
cd quickstart-runs/case17-joule-heating

# Run with the combined module application (includes heat_transfer module)
/path/to/combined-opt -i case17_joule_heating.i
```

The run produces:

- `case17_joule_heating_out.e` — Exodus file with V and T fields at all 20
  timesteps plus the initial condition frame
- `case17_joule_heating_out.csv` — CSV file with max_T, avg_T, and max_V vs.
  time (21 rows including t=0)

Expected runtime: a few seconds on a single core.

---

## Expected Results

### Electric Potential (V field)

V is linear in x and constant in y: V(x) = 10 * (1 - x/2). It reaches steady
state immediately (no time dependence). The maximum is 10 V on the left face
and the minimum is 0 V on the right face. The postprocessor `max_V` should
read approximately 10.0 throughout the run.

### Temperature (T field)

At t=0, T=300 K everywhere. Joule heating adds ~25 MW/m^3 uniformly. Heat
diffuses toward the cooled electrodes at x=0 and x=2. The temperature profile
becomes parabolic, symmetric about x=1.

At t=5 s the postprocessors should show approximately:

| Postprocessor | Approximate value |
|---------------|-------------------|
| max_T         | ~330 K            |
| avg_T         | ~320 K            |
| max_V         | ~10.0 V           |

`max_T` grows roughly linearly in time (6-7 K/s) because the conduction path
to the electrodes is long (1 m) and thermal diffusivity alpha = k/(rho*cp) =
50/(8000*500) = 1.25e-5 m^2/s is moderate — the conductor heats up much
faster than it can cool through the electrodes.

The `avg_T` grows slightly more slowly than `max_T` because the electrode
boundary conditions hold the temperature at 300 K near x=0 and x=2, pulling
the spatial average below the peak.

### Visualization in ParaView

1. Open `case17_joule_heating_out.e` in ParaView and click Apply.
2. Select the `V` field and step through time — the V distribution should look
   identical at every frame (steady Laplace solution).
3. Select the `T` field and animate. Watch the temperature build from a flat
   300 K surface to a parabolic arch centred on x=1.
4. Use Filters > Plot Over Line (from (0,0.5,0) to (2,0.5,0)) to extract the
   midline T profile. At late time it should be a symmetric parabola with the
   peak at x=1 and T=300 K at both ends.

---

## Key Takeaways

- **One-way coupling via a material object**: `ElectromagneticHeatingMaterial`
  bridges the electric and thermal physics. It reads the V solution and writes
  a heat source property that the temperature equation consumes. This material-
  mediated coupling is a common MOOSE pattern.

- **Kernel reuse across physics**: `ADHeatConduction` implements any scalar
  Laplacian-type equation. Setting `thermal_conductivity = electrical_conductivity`
  turns it into the charge-conservation kernel for V. The same mathematical
  kernel serves two different physics.

- **AD objects for exact Jacobians**: Using `AD` prefixed objects (ADDirichletBC,
  ADHeatConduction, ADJouleHeatingSource, ADGenericConstantMaterial) lets MOOSE
  assemble the exact Jacobian via dual-number automatic differentiation. No
  hand-coded Jacobian entries are needed for either the diagonal or off-diagonal
  blocks.

- **SMP full preconditioning**: With two coupled variables the off-diagonal
  Jacobian block (dR_T/dV) is non-trivial because the Joule source depends on
  grad(V). Setting `full = true` in SMP ensures this cross-variable derivative
  is included in the preconditioner, enabling fast Newton convergence.

- **Uniform Joule heating in a rectangular conductor**: When sigma is constant
  and the electric potential is purely 1D (linear in x), the heat source Q is
  spatially uniform. The resulting T profile is parabolic — the classic result
  for a heated rod with cooled ends.

---

## Experiments to Try

### Experiment 1: Observe the Parabolic Temperature Profile

Run the simulation and use ParaView's Plot Over Line filter along y=0.5
(the domain midline). At t=5 s the T curve should be a symmetric parabola.
Fit the data to T = 300 + A*x*(2-x). Compare the coefficient A with the
analytical prediction A = Q/(2*k) = 25e6/(2*50) = 250,000. At t=5 s the
transient profile won't reach this steady-state value, but the parabolic
shape should be visible.

### Experiment 2: Increase the Voltage

Change `V_left` from 10 to 20 V. The electric field doubles, so the Joule
source quadruples (Q scales as V^2). The temperature rise rate should be four
times faster. Verify by comparing max_T at t=5 s to the original run.

Expected outcome: max_T grows from ~330 K to ~420 K at t=5 s.

### Experiment 3: Add Temperature-Dependent Conductivity

Replace the constant electrical conductivity with one that decreases with
temperature (metals become less conductive when hot). Remove the constant
`[electrical]` material block and instead use `ADElectricalConductivity` from
the heat_transfer module (which models copper). This turns the problem into
two-way coupling: hotter material conducts less, reducing the Joule source.
This negative feedback slows the temperature rise and the system reaches a
true coupled steady state.

### Experiment 4: Insulate One Electrode

Remove the `T_right` boundary condition (comment out or delete the block). The
right electrode no longer acts as a heat sink. All the Joule heat must exit
through the left electrode. At steady state the temperature profile becomes
linear (no cooling from the right). Observe how max_T rises much higher than
in the original two-sided cooling case.

### Experiment 5: Add a Constriction with Two Subdomains

Use `SubdomainBoundingBoxGenerator` to create a narrow region (e.g., a thin
strip in the middle of the domain at lower electrical conductivity). Current
crowding in the constriction produces locally higher Joule heating. This is the
effect that causes fuses to blow and circuit board traces to overheat at narrow
points.
