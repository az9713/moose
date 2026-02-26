# Case 21: Thermo-Mechanical Bimetallic Strip

## Overview

A bimetallic strip is one of the most elegant demonstrations of differential thermal
expansion: bond two metals together, heat them, and the strip bends because each metal
tries to expand by a different amount but cannot separate. The result is a curved shape
whose curvature is proportional to the temperature change and the mismatch in thermal
expansion coefficients.

This case models that phenomenon in MOOSE using:

- **Two mesh subdomains** — steel (bottom, block 0) and aluminum (top, block 1)
- **Block-restricted materials** — each block gets its own elasticity tensor and thermal
  eigenstrain, with physically correct properties for each metal
- **Prescribed temperature field** — temperature is an `AuxVariable` set to 500 K
  everywhere; no heat equation is solved because the heating is uniform
- **SolidMechanics QuasiStatic action** — creates displacement variables, residuals,
  and stress output automatically via a high-level `[Physics]` block

Strip geometry: 10 m long, 1 m tall, with the interface at y = 0.5 m.
The left end is pinned (zero displacement). The right end is free.

| Property              | Steel (block 0) | Aluminum (block 1) |
|-----------------------|-----------------|--------------------|
| Young's modulus       | 200 GPa         | 70 GPa             |
| Poisson's ratio       | 0.3             | 0.33               |
| Thermal expansion     | 12e-6 /K        | 23e-6 /K           |
| Stress-free temp      | 300 K           | 300 K              |
| Final temperature     | 500 K           | 500 K              |
| Free thermal strain   | 0.24%           | 0.46%              |

---

## The Physics

### Differential Thermal Expansion

When a material is heated by deltaT, it wants to expand in every direction by a strain
equal to alpha * deltaT, where alpha is the coefficient of thermal expansion (CTE).
This is the **free thermal strain** — what the material would do if it were unconstrained.

For the two metals heated by deltaT = 200 K:

```
Steel:    epsilon_free = 12e-6 * 200 = 0.0024  (0.24%)
Aluminum: epsilon_free = 23e-6 * 200 = 0.0046  (0.46%)
```

Aluminum wants to expand roughly twice as much as steel in the longitudinal direction.
But the two layers are bonded at the interface: they must deform together. Neither can
achieve its free expansion. Instead, a compromise deformation emerges — the strip bends.

### Why the Strip Bends (and in Which Direction)

Consider the bonded interface. Aluminum (top) wants to be longer than steel (bottom).
Since they cannot separate, the interface is forced to a common length somewhere between
the two free expansions. The aluminum layer ends up shorter than it wants to be
(compressive stress along x) and the steel layer ends up longer than it wants to be
(tensile stress along x). This self-equilibrating stress state — compression above the
interface, tension below — produces a **bending moment** that curves the strip.

The direction: aluminum is on top and is compressed, steel is on bottom and is in
tension. A beam in bending curves toward the compressed side. So the strip bends
**downward** (toward the steel, away from the aluminum). The free tip at x=10 m moves
in the negative y direction.

```
Before heating (straight strip):
  y=1  +--aluminum (top)----------------------------------+
  y=0.5| - - - - - - - - interface - - - - - - - - - - - |
  y=0  +--steel (bottom)----------------------------------+
       x=0                                              x=10

After heating (bends toward steel side):
  y=1  +-aluminum-_
  y=0.5|interface   \__
  y=0  +-steel---------\___
       x=0              x=10  <-- tip deflects downward (neg y)
```

### Timoshenko's Bimetallic Strip Formula

For a two-layer bimetallic strip under uniform heating, Timoshenko (1925) derived the
exact curvature under the assumption of Euler-Bernoulli beam theory:

```
kappa = 6 * (alpha_2 - alpha_1) * deltaT * (1 + m)^2
        ---------------------------------------------------
        h * [ 3*(1+m)^2 + (1 + m*n)*(m^2 + 1/(m*n)) ]

where:
  alpha_1 = CTE of layer 1 (steel, bottom)
  alpha_2 = CTE of layer 2 (aluminum, top)
  deltaT  = temperature change (200 K)
  h       = total thickness (1.0 m)
  m       = h_2/h_1 = thickness ratio = 0.5/0.5 = 1
  n       = E_2/E_1 = modulus ratio = 70e9/200e9 = 0.35
```

For equal thickness (m=1), the formula simplifies to:

```
kappa = 6 * (alpha_2 - alpha_1) * deltaT
        ------------------------------------------
        h * [ 3 + 2*(1+n) + 1/n ]   ... (m=1 case)
```

Substituting values:
- alpha_2 - alpha_1 = (23 - 12)e-6 = 11e-6 /K
- deltaT = 200 K
- h = 1.0 m
- n = 70/200 = 0.35

```
Numerator:   6 * 11e-6 * 200 = 0.01320
Denominator: 1.0 * [3 + 2*(1+0.35) + 1/0.35]
           = 1.0 * [3 + 2.70 + 2.857]
           = 1.0 * 8.557
kappa      = 0.01320 / 8.557 = 0.001543 m^{-1}
```

For a cantilever of length L=10 m, the tip deflection under uniform curvature is:

```
delta_tip = kappa * L^2 / 2 = 0.001543 * 100 / 2 = 0.0772 m
```

The MOOSE result should be close to this value (it will differ somewhat because
Timoshenko's formula assumes pure bending with no end effects, and the 2D plane-stress
solution includes Poisson effects).

### The Eigenstrain Approach

MOOSE implements thermal expansion via the **eigenstrain** mechanism. An eigenstrain is
a strain that is imposed on the material (rather than derived from displacement
gradients). The thermal eigenstrain is:

```
epsilon_thermal = alpha * (T - T_ref) * I
```

where I is the identity tensor (isotropic expansion in all directions). This eigenstrain
is subtracted from the total strain to get the **elastic strain**, which drives stress:

```
epsilon_elastic = epsilon_total - epsilon_thermal
sigma = C : epsilon_elastic
```

When a material heats uniformly and is unconstrained, epsilon_total = epsilon_thermal,
so epsilon_elastic = 0 and sigma = 0 (no stress in a freely expanding body). In the
bimetallic strip, the constraint from the bonded interface forces epsilon_total to
deviate from epsilon_thermal, producing nonzero elastic strain and stress.

---

## Input File Walkthrough

### `[Mesh]`

The mesh pipeline builds a 40x8 strip (10 m by 1 m) and then relabels the top four rows
of elements as block 1 (aluminum). The interface at y=0.5 falls exactly on element
boundaries because ny=8 gives row edges at y = 0, 0.125, 0.25, 0.375, 0.5, ...

```
[gen]          -- 40x8 uniform mesh, all elements in block 0 (steel)
[top_block]    -- SubdomainBoundingBoxGenerator: elements with centroid
                  y > 0.5 become block 1 (aluminum)
```

### `[GlobalParams]`

```
displacements = 'disp_x disp_y'
```

This global declaration tells the SolidMechanics system the names of the displacement
variables. Every object in the system (kernels, materials, postprocessors) reads this
rather than having the names repeated in each block.

### `[Physics/SolidMechanics/QuasiStatic]`

The `QuasiStatic` action is the high-level entry point for static or quasi-static
solid mechanics in MOOSE. It automatically:

- Creates the `disp_x` and `disp_y` nodal variables (`add_variables = true`)
- Adds the momentum balance kernels (divergence of stress)
- Adds a `ComputeSmallStrain` material that computes total strain from displacements
- Requests output of derived quantities (`vonmises_stress`, `stress_xx`, `stress_yy`)

The `eigenstrain_names = 'thermal_eigenstrain'` parameter tells the strain calculator
to subtract the named eigenstrain(s) when computing elastic strain.

`strain = SMALL` selects the small-strain (linear) kinematic assumption, appropriate
here because the tip deflection (~0.08 m) is much smaller than the strip length (10 m).

### `[AuxVariables]` — Prescribed Temperature

```
[T]
  initial_condition = 500
[]
```

Temperature is not solved; it is prescribed. `AuxVariable` creates a field variable
that is set directly (not through a residual equation). Setting `initial_condition = 500`
fills every node with T = 500 K before the solve begins. Because the executioner is
`Steady` and no `AuxKernels` update T, this value is fixed throughout the solve.

The thermal expansion materials reference this `T` variable. With T_ref = 300 K,
every element sees deltaT = 200 K.

### `[Materials]`

Six material objects are defined, three per block:

1. `ComputeIsotropicElasticityTensor` — fills the 4th-order elasticity tensor C from
   Young's modulus and Poisson's ratio using the standard isotropic formula.

2. `ComputeThermalExpansionEigenstrain` — computes the thermal eigenstrain tensor
   at each quadrature point as `alpha * (T - T_ref) * I`. The `eigenstrain_name`
   must match what was declared in the `[Physics]` block.

3. `ComputeLinearElasticStress` — computes `sigma = C : epsilon_elastic`. This is
   the Hookean stress calculation; it reads the elasticity tensor and the elastic strain
   (total strain minus eigenstrain) and returns the stress tensor.

All six objects carry `block = 0` or `block = 1` to prevent them from being applied
to the wrong region.

### `[BCs]`

Both displacement components are set to zero on the left boundary, pinning the entire
left edge. All other boundaries are traction-free (natural BCs, zero normal stress),
which is the default when no BC is specified.

### `[Postprocessors]`

- `max_disp_y` / `min_disp_y` — extreme vertical displacements over all nodes.
  `min_disp_y` captures the downward tip deflection.
- `tip_disp_y` — vertical displacement at the geometric midpoint of the right edge
  (x=10, y=0.5). This is the most direct comparison with Timoshenko's formula.
- `max_vonmises` — peak von Mises stress, expected near the bonded interface.

---

## Running the Simulation

This case requires the `combined` application (or any MOOSE app built with the
`SOLID_MECHANICS` module enabled).

```bash
cd quickstart-runs/case21-bimetallic-strip

# With the combined app (adjust path as needed):
combined-opt -i case21_bimetallic_strip.i

# With mpirun for parallel:
mpirun -n 4 combined-opt -i case21_bimetallic_strip.i
```

The solve is a single linear Newton step (the problem is linear: small strain,
linear elasticity, prescribed temperature). Convergence should be immediate.

Output files produced:
- `case21_bimetallic_strip_out.e` — Exodus mesh + displacement + stress fields
- `case21_bimetallic_strip_out.csv` — postprocessor values (tip deflection, max stress)

---

## Expected Results

### Displacement Field

The strip bends downward (negative y direction). The deformed shape resembles an arc,
with maximum displacement at the free tip (x=10 m) and zero displacement along the
pinned left edge.

Approximate values for deltaT = 200 K:

| Quantity          | Timoshenko estimate | Expected MOOSE range |
|-------------------|---------------------|----------------------|
| Tip deflection    | -0.077 m            | -0.06 to -0.09 m     |
| Max disp_x        | ~0.024 m (axial)    | 0.02 to 0.03 m       |
| Peak von Mises    | ~500 MPa (steel)    | 200 to 600 MPa       |

The tip deflection will be in the negative y direction (`min_disp_y` will be the
relevant postprocessor). The Timoshenko estimate of -0.077 m is for pure bending
with no Poisson effects; the 2D result will be slightly different.

### Stress Distribution

The stress pattern is the distinguishing signature of the bimetallic effect:

- **Along x (axial stress)**: tension in the steel layer, compression in the aluminum
  layer. The magnitude increases away from the neutral axis. Peak stress at the top
  surface of the aluminum and the bottom surface of the steel.
- **At the bonded interface (y=0.5)**: the axial stress jumps discontinuously because
  the elastic moduli differ (E_steel = 200 GPa vs E_aluminum = 70 GPa). The strain is
  continuous across the interface (bonded), but stress = E * strain differs.
- **Von Mises stress**: concentrates at the interface and at the corners of the fixed
  left boundary (stress concentration from the pin constraint).

### Physical Interpretation

The bimetallic strip converts a temperature change into a mechanical displacement with
no moving parts. This principle underlies:

- **Thermostats** — the strip opens or closes an electrical contact when temperature
  crosses a threshold
- **Thermal circuit breakers** — bimetallic element trips the breaker under overcurrent
  heating
- **Temperature-compensating mechanisms** — the strip cancels thermal expansion in
  precision instruments

In all cases, the useful output is the deflection at the free end, which is linearly
proportional to temperature change (for small deflections). The MOOSE simulation
captures this linear relationship exactly.

---

## Key Takeaways

- **Eigenstrain** is MOOSE's mechanism for prescribed strains (thermal, swelling,
  plasticity); it enters the constitutive update as a subtraction from total strain.
- **Block-restricted materials** allow each mesh region to carry independent material
  properties while sharing the same set of governing equations (kernels).
- **AuxVariable for prescribed fields** avoids solving unnecessary PDEs — when
  temperature is uniform and known, declare it as an `AuxVariable` with an
  `initial_condition` rather than adding a heat equation.
- **The QuasiStatic action** is the recommended entry point for solid mechanics in
  MOOSE; it wires together elasticity, strain, stress, and output automatically.
- **Differential thermal expansion** produces self-equilibrating stress states
  (no external loads), but the stresses can be very large — hundreds of MPa for
  engineering metals over moderate temperature changes.
- **Timoshenko's formula** provides a clean analytical check for the tip deflection,
  confirming that the MOOSE result is physically reasonable.
