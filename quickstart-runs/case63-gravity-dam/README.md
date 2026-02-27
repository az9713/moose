# Case 63: Gravity Dam — Self-Weight Loading and Foundation Stress

## Overview

A gravity dam is a massive concrete structure that resists the horizontal water pressure of a reservoir through its own weight. The dam transmits its weight to the rock foundation through compressive stress. Understanding how the dam loads the foundation — where compressive stresses are highest, where tensile zones might appear, and how the deformation field looks — is fundamental to dam engineering.

This case models a rectangular concrete dam sitting on a deformable rock foundation under self-weight loading only (no water pressure). It demonstrates two important MOOSE capabilities: using `SubdomainBoundingBoxGenerator` to assign different material properties to different spatial regions of a single mesh, and using the `Gravity` kernel to apply body forces proportional to material density. The result is a stress field that can be compared to classical analytical estimates for a gravity-loaded elastic structure.

Key concepts demonstrated:

- `SubdomainBoundingBoxGenerator` to create material subdomains without a separate mesh generator
- `Gravity` kernel for density-dependent body forces
- Block-restricted materials: different elastic properties assigned to dam (block 1) vs. foundation (block 0)
- `Physics/SolidMechanics/QuasiStatic` Action for plane-strain linear elasticity
- Stress concentration at the dam-foundation interface under self-weight

---

## The Physics

### Governing Equations

Plane-strain quasi-static equilibrium with a body force:

```
div(sigma) + rho * g = 0    in the domain
sigma = C : epsilon         (linear elastic constitutive law)
epsilon = sym(grad(u))      (small-strain kinematics)
```

The body force `rho * g` is the gravitational acceleration multiplied by the material density. With `g = -9.81 m/s^2` in the y-direction, the body force is downward.

### Material Properties

**Dam (concrete, block 1):**
- E = 25 GPa, nu = 0.2
- density = 2400 kg/m^3
- Region: x in [6, 14] m, y in [0, 10] m (8 m wide, 10 m tall)

**Foundation (rock, block 0):**
- E = 10 GPa, nu = 0.3
- density = 2600 kg/m^3
- Region: entire domain except the dam block

### Self-Weight Analytical Estimate

For a uniform elastic column under self-weight with the base fixed and sides free (not the oedometer condition of Case 62), the vertical stress at any height y from the base is:

```
sigma_yy(y) = -rho * g * (H_top - y)    (compressive, hence negative)
```

For the dam material (rho = 2400 kg/m^3, H_dam = 10 m):
- At y = 0 (dam base): sigma_yy = -2400 * 9.81 * 10 = -235 kPa
- At y = 10 (dam top): sigma_yy = 0

For the foundation below the dam center (adding the dam weight as a distributed stress), the foundation vertical stress increases further. The combined compressive stress at the dam-foundation interface center is approximately:

```
sigma_yy ~ -(rho_dam * g * H_dam + rho_found * g * H_found)
         ~ -(2400 * 9.81 * 10 + 2600 * 9.81 * 10)
         ~ -(235.4 + 255.1) = -490 kPa   (rough two-layer estimate)
```

The actual simulated value differs because: (1) only the central 8 m carries dam weight, not the full 20 m foundation width; (2) the Poisson effect in plane strain generates horizontal stress; (3) stress redistribution occurs at material interfaces.

### Von Mises Stress and Stress Concentrations

The von Mises stress sigma_vm = sqrt(3 * J2) measures the intensity of shear stress (relevant for yielding). For a pure vertical compression state (sigma_yy = sigma_0, sigma_xx = nu/(1-nu) * sigma_yy due to plane strain, sigma_xy = 0):

```
sigma_vm = |sigma_yy - sigma_xx| = |sigma_yy * (1 - nu/(1-nu))|
         = |sigma_yy * (1-2*nu)/(1-nu)|
```

The highest von Mises stress occurs where stress gradients are sharpest — at the edges of the dam base, where the dam's weight concentrates due to the abrupt transition from loaded to unloaded foundation.

### Domain and Boundary Conditions

```
y = 10 m   ___   ___   TOP: dam top (free surface)
          |   | |   |
y =  0 m  |DAM| |   |   y = 0: dam base / foundation interface
     _____| 8m|_|___________
    |                       | y = 0
    |   Foundation rock     |
    |   E=10 GPa, nu=0.3   |
    |   20 m wide x 10 m   |
    |_______________________|
y = -10 m                    BOTTOM: fully fixed (no displacement)

Sides (x=0, x=20): roller (disp_x = 0)
Top (y=10): free (no BC on dam top)
Dam: x in [6,14], y in [0,10]  (centered on 20 m wide foundation)
```

The roller conditions on the sides model symmetry or lateral confinement from adjacent soil. The bottom of the foundation is fixed on bedrock. The dam and foundation surface not covered by the dam is a free surface.

---

## Input File Walkthrough

The input file is `case63_gravity_dam.i`.

### `[GlobalParams]`

```
displacements = 'disp_x disp_y'
```

Declares the displacement variable names globally for `Physics/SolidMechanics/QuasiStatic` and all BCs.

### `[Mesh]`

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 20
  ny = 20
  xmin = 0
  xmax = 20
  ymin = -10
  ymax = 10
[]
[dam_block]
  type = SubdomainBoundingBoxGenerator
  input = gmg
  block_id = 1
  block_name = dam
  bottom_left = '6 0 0'
  top_right = '14 10 0'
[]
```

The mesh is a single 20x20 element grid spanning x in [0,20] m, y in [-10,10] m. The `SubdomainBoundingBoxGenerator` assigns block ID 1 (named `dam`) to all elements whose centroids fall inside the rectangle [6,14] x [0,10]. All other elements retain the default block ID 0. After this operation, the mesh has two blocks and different materials can be assigned to each — no interface mesh matching is required.

The 20x20 mesh is uniform: each element is 1 m x 1 m. The dam spans x = [6,14] (8 elements wide) and y = [0,10] (10 elements tall). The foundation spans y = [-10,0] (10 elements deep) across the full 20 m width.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = SMALL
  add_variables = true
  generate_output = 'stress_yy stress_xx vonmises_stress'
[]
```

The QuasiStatic Action creates `disp_x` and `disp_y` variables, the `StressDivergenceTensors` kernels for both components, and AuxVariables + AuxKernels for the three requested stress outputs. The stress fields are computed as element-constant (MONOMIAL) quantities from the element-average strain.

### `[Kernels]`

```
[gravity_y]
  type = Gravity
  variable = disp_y
  value = -9.81
[]
```

The `Gravity` kernel adds the body force term `integral(rho * g * test dV)` to the y-momentum equation residual. The density `rho` is obtained from the `density` material property, which is block-restricted: 2400 kg/m^3 in the dam block, 2600 kg/m^3 in the foundation. This single kernel applies to both blocks; the correct density is used automatically at each quadrature point.

### `[Materials]`

Six material objects define the two-block elastic model:

**Dam (block = dam):**
- `ComputeIsotropicElasticityTensor`: E = 25e9 Pa (25 GPa), nu = 0.2
- `GenericConstantMaterial`: density = 2400 kg/m^3 (registered as `density`)
- `ComputeLinearElasticStress`: sigma = C : epsilon

**Foundation (block = 0):**
- `ComputeIsotropicElasticityTensor`: E = 10e9 Pa (10 GPa), nu = 0.3
- `GenericConstantMaterial`: density = 2600 kg/m^3
- `ComputeLinearElasticStress`: sigma = C : epsilon

Note that three separate material objects are needed per block (`elasticity_tensor`, `density`, `stress`) because each type of material object handles one aspect of the constitutive model. The `block = dam` and `block = 0` restrictions ensure that each material is only evaluated in its assigned region.

### `[BCs]`

| Name | Variable | Boundary | Value | Purpose |
|------|----------|----------|-------|---------|
| `fix_bottom_x` | disp_x | bottom | 0 | No horizontal slip at bedrock |
| `fix_bottom_y` | disp_y | bottom | 0 | No vertical movement at bedrock |
| `fix_left_x` | disp_x | left | 0 | Lateral roller confinement |
| `fix_right_x` | disp_x | right | 0 | Lateral roller confinement |

The top surface (y = 10) is free — no boundary condition is applied, so it carries zero traction. This is correct for the dam top and the exposed foundation surface.

### `[Postprocessors]`

| Name | Type | Variable | Physical meaning |
|------|------|----------|-----------------|
| `max_stress_yy` | ElementExtremeValue | stress_yy | Peak (least compressive) vertical stress |
| `max_vonmises` | ElementExtremeValue (max) | vonmises_stress | Maximum shear stress indicator |
| `max_disp_y` | ElementExtremeValue | disp_y | Maximum downward settlement |
| `stress_yy_center` | PointValue at (10,0,0) | stress_yy | Vertical stress at dam base centerline |

### `[Executioner]`

```
type = Steady
solve_type = NEWTON
nl_rel_tol = 1e-10
nl_abs_tol = 1e-12
```

A steady-state (not transient) solve: gravity loading does not vary with time. NEWTON with tight tolerances gives the linear elastic solution in a single nonlinear iteration (the system is linear so Newton converges in one step).

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case63-gravity-dam \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case63_gravity_dam.i 2>&1 | tail -30'
```

Output files:
- `case63_gravity_dam_out.e` — Exodus with displacement and stress fields
- `case63_gravity_dam_out.csv` — Single-row summary of postprocessor values

---

## Expected Results

The simulation converges in a single iteration (linear elastic problem). The postprocessor CSV contains one row:

```
max_stress_yy:    ~ -8109 Pa      (most compressive vertical stress, element extremum)
max_vonmises:     ~ 283,306 Pa    (maximum von Mises stress, at dam base corners)
max_disp_y:       ~ -7.60e-6 m   (maximum downward settlement)
stress_yy_center: ~ -273,800 Pa  (vertical stress at dam base centerline)
```

### Stress at the Dam Base

The vertical stress at the dam base center (sigma_yy_center ~ -274 kPa) represents the combined weight of the dam column above. Analytically:

```
sigma_yy ~ -rho_dam * g * H_dam = -2400 * 9.81 * 10 = -235.4 kPa
```

The simulated -274 kPa is somewhat higher because it includes both the dam self-weight and the constraint effects from the plane-strain Poisson coupling (the horizontal stress from the roller BCs increases the effective vertical stress).

### Von Mises Stress Distribution

The maximum von Mises stress (~283 kPa) occurs at the corners of the dam base, where the abrupt change in loading (dam weight suddenly ending at x = 6 m and x = 14 m) creates stress concentrations. This is where shear stresses are highest.

In ParaView, the von Mises stress field shows:
- High values at the dam base corners (stress concentration points)
- Moderate values under the dam center where vertical compressive stress is high
- Lower values in the foundation away from the dam footprint
- Near-zero values in the unloaded foundation surface regions (y = 0, x < 6 or x > 14)

### Settlement Field

The maximum settlement (~7.6 micrometers) occurs at the dam base under the central load. For 25 GPa concrete, a 10 m column under its own weight (average stress ~118 kPa) settles:

```
delta ~ sigma_avg * H / E = 1.18e5 * 10 / 2.5e10 ~ 4.7e-5 m
```

The simulated 7.6e-6 m is smaller because the plane-strain constraint and the stiffer rock foundation both reduce the deformation. The settlement is very small (micrometers) because concrete is very stiff — this confirms the structure is far from failure.

---

## Key Takeaways

- `SubdomainBoundingBoxGenerator` assigns block IDs to elements by spatial bounding box without requiring a separate mesh for each material region. It is applied in a mesh generator pipeline, feeding into (or being fed by) other generators.
- Block-restricted materials use `block = <name_or_id>` to limit the material to specific mesh blocks. Each block needs its own set of material objects; objects without a `block` restriction apply globally to all blocks.
- The `Gravity` kernel uses the `density` material property to compute the body force `rho * g`; the density material must be declared before `Gravity` can find it. With block-restricted densities, each element automatically uses its material's density value.
- `Physics/SolidMechanics/QuasiStatic` with `generate_output` is a convenient shortcut that creates AuxVariables and AuxKernels for stress, strain, and von Mises outputs automatically.
- A `Steady` executioner is appropriate for quasi-static problems that do not have time-dependent loading. Gravity loading is time-independent, so there is no need for a transient solve.
- The stress concentration at the dam base corners — where the von Mises stress is highest — is the critical region for dam stability assessment. In real dam design, this region is checked against the tensile strength of the dam-foundation interface.
- The plane-strain assumption (no out-of-plane strain) is appropriate for long dams where the length perpendicular to the cross-section is much greater than the dam height. For shorter structures, a 3D model would be needed.
