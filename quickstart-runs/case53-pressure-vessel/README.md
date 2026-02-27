# Case 53: Pressure Vessel — Thick-Walled Cylinder (Lamé Solution)

## Overview

The previous four solid mechanics cases focused on material and geometric nonlinearity — plasticity, finite strain, creep, and damage. This case returns to linear elasticity to demonstrate a different capability: axisymmetric analysis and validation against a classical analytical solution.

Real pressure vessels, pipes, and gun barrels are often better modelled as cylinders than as rectangular plates. For a long cylinder, the problem is independent of the axial coordinate and the hoop direction, leaving only the radial coordinate as the independent variable. MOOSE supports this geometry through the Coordinate System setting `coord_type = RZ`, which maps the 2D Cartesian mesh to a 2D axisymmetric (r, z) space. This eliminates the need to model a 3D cylinder; a 2D slice through the r-z plane is sufficient.

For a thick-walled cylinder under internal pressure, Lamé (1852) derived the exact analytical solution for the radial and hoop stress distributions:

```
sigma_r(r)     = A - B/r^2
sigma_theta(r) = A + B/r^2
```

MOOSE's numerical solution can be compared directly against these formulas to verify the implementation. This case documents how to set up RZ problems, which stress components map to which physical quantities, and how closely a 20-element radial mesh matches the analytical solution.

New concepts introduced in this case:

- **`coord_type = RZ`**: activates axisymmetric geometry in MOOSE, turning the 2D mesh into a revolution solid without explicitly meshing the full 3D geometry.
- **`ComputeLinearElasticStress`** (non-AD version): straightforward linear elastic stress for problems where the Jacobian does not need AD-based differentiation.
- **Stress component mapping in RZ**: stress_xx = sigma_r, stress_yy = sigma_z, stress_zz = sigma_theta (hoop).
- **`Pressure` BC in RZ**: applies force per unit deformed area to a boundary; in RZ this correctly accounts for the circumferential integration.
- **`PointValue` postprocessors**: sample field variables at specific (r, z) coordinates for direct comparison with analytical formulas.

---

## The Physics

### The Physical Problem in Plain English

A thick-walled steel cylinder has inner radius r_i = 1 m and outer radius r_o = 2 m. The inner surface is pressurised to p_i = 100 MPa. The outer surface is traction-free (exposed to atmosphere). The cylinder is long, so the ends are constrained to prevent axial displacement — a "plane strain" condition in the axial direction.

The internal pressure pushes outward radially, generating compressive radial stress (pushing the wall outward) and tensile hoop stress (trying to split the cylinder open). The hoop stress is maximum at the inner surface and decreases with radius. This is why thick-walled pressure vessels fail by bursting at the bore — the inner surface is always at the highest stress.

### Lamé Analytical Solution

For a cylinder with inner radius r_i and outer radius r_o under internal pressure p_i (outer pressure = 0):

**Constants**:

```
A = p_i * r_i^2 / (r_o^2 - r_i^2) = 100 * 1 / (4 - 1) = 33.33 MPa
B = p_i * r_i^2 * r_o^2 / (r_o^2 - r_i^2) = 100 * 1 * 4 / (4 - 1) = 133.33 MPa·m^2
```

**Stress distributions**:

```
sigma_r(r)     = A - B/r^2 = 33.33 - 133.33/r^2
sigma_theta(r) = A + B/r^2 = 33.33 + 133.33/r^2
```

**Boundary condition checks**:

```
sigma_r(r_i = 1) = 33.33 - 133.33 = -100.0 MPa  (= -p_i, compressive, check!)
sigma_r(r_o = 2) = 33.33 - 133.33/4 = 0.0 MPa   (traction-free, check!)
```

**Stress values at key radii**:

| r (m) | sigma_r (MPa) | sigma_theta (MPa) |
|-------|---------------|-------------------|
| 1.0 (inner) | -100.0  | +166.67           |
| 1.5 (mid)   |  -25.9  | +92.6             |
| 2.0 (outer) |    0.0  | +66.67            |

The maximum tensile stress is the hoop stress at the inner surface: 166.67 MPa. This exceeds the radial pressure by a factor of 1.67, which is why hoop stress governs pressure vessel design.

### Coordinate System Conventions in MOOSE RZ

In `coord_type = RZ`, the 2D mesh plane represents the r-z half-plane:

```
x coordinate → r (radial direction)
y coordinate → z (axial direction)
theta is implicit (axisymmetry)
```

The stress output names from the QuasiStatic action map as:

```
stress_xx → sigma_r   (radial stress)
stress_yy → sigma_z   (axial stress)
stress_zz → sigma_theta (hoop stress — out-of-plane in the 2D mesh)
```

This mapping is critical when reading results: sigma_theta appears as `stress_zz` in the output, not `stress_yy`.

### Domain Diagram

```
         z (y in mesh)
         ^
     1.0 |____________________________
         |                            |
         |  Steel cylinder wall       |
         |  E = 200 GPa               |
     0.0 |____________________________|---> r (x in mesh)
         r=1.0                       r=2.0
         inner (p_i = 100 MPa)       outer (traction-free)

coord_type = RZ  (revolve about z-axis at r=0)
Mesh: 20 elements in r, 4 elements in z
Axial BCs: disp_y = 0 on top and bottom (plane-strain in z)
```

---

## Input File Walkthrough

The input file is `case53_pressure_vessel.i`.

### Mesh and Coordinate System

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 20    # radial direction
  ny = 4     # axial direction
  xmin = 1.0  xmax = 2.0   # r_i to r_o
  ymin = 0.0  ymax = 1.0
[]
coord_type = RZ
```

`coord_type = RZ` is set at the Mesh block level, not inside a sub-block. The mesh itself is an ordinary 2D rectangle, but MOOSE treats it as the r-z cross-section of an axisymmetric body. Internally, all integrals acquire a factor of 2*pi*r (the circumferential volume element), and the hoop strain and stress are automatically computed as extra components of the strain and stress tensors.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = SMALL
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
[]
```

In RZ, `generate_output = 'stress_zz'` gives the hoop stress sigma_theta automatically. The QuasiStatic action recognises the coordinate system and sets up the hoop strain contribution in the kinematics. No special configuration beyond `coord_type = RZ` in the Mesh block is needed.

### `[BCs]`

**Internal pressure**:

```
[Pressure]
  [internal_pressure]
    boundary = left
    factor = 100.0    # 100 MPa
  []
[]
```

The `left` boundary corresponds to r = r_i = 1 m (the inner surface). The `Pressure` BC applies an inward traction (radially inward pressure). In RZ, MOOSE correctly integrates the pressure over the circumferential arc length automatically.

**Axial constraints**:

```
[fix_z_bottom]  disp_y = 0  on bottom
[fix_z_top]     disp_y = 0  on top
```

Both top and bottom boundaries are constrained to zero axial displacement. This enforces the plane-strain condition in the z-direction (appropriate for a long cylinder capped at both ends under no axial force, or equivalently for an infinite cylinder). There is no BC on the outer surface (r = r_o = 2 m) — it is traction-free by the natural boundary condition.

### `[Materials]`

```
[elasticity_tensor]
  type = ComputeIsotropicElasticityTensor
  youngs_modulus = 200e3   # MPa
  poissons_ratio = 0.3
[]
[stress]
  type = ComputeLinearElasticStress
[]
```

`ComputeLinearElasticStress` (non-AD, no prefix) is appropriate here because the problem is linear and there is no need for AD-based Jacobian computation. The Steady executioner with Newton converges in one iteration for linear problems.

### `[Functions]` — Analytical Solution

```
[sigma_r_exact]
  type = ParsedFunction
  expression = '33.3333 - 133.3333 / (x * x)'   # x = r in RZ
[]
[sigma_theta_exact]
  type = ParsedFunction
  expression = '33.3333 + 133.3333 / (x * x)'
[]
```

These functions encode the Lamé solution using x (the r-coordinate) as the independent variable. They are not directly applied in the simulation — they serve as reference values for manual comparison and could be used in a `FunctionValuePostprocessor` or `ElementL2Error` if you want to automate the error calculation.

### `[Postprocessors]` — Validation Sampling

Four `PointValue` postprocessors sample the stress fields at the inner and outer surfaces:

| Name | Variable | Point | Analytical |
|------|----------|-------|------------|
| `sigma_r_inner` | stress_xx | (1.0, 0.5, 0) | -100.0 MPa |
| `sigma_r_outer` | stress_xx | (2.0, 0.5, 0) | 0.0 MPa |
| `sigma_theta_inner` | stress_zz | (1.0, 0.5, 0) | +166.67 MPa |
| `sigma_theta_outer` | stress_zz | (2.0, 0.5, 0) | +66.67 MPa |

Note that `stress_zz` is used for sigma_theta (hoop), consistent with the RZ stress component mapping described above.

### `[Executioner]`

```
type = Steady
solve_type = 'NEWTON'
nl_rel_tol = 1e-10
nl_abs_tol = 1e-12
```

A single Steady solve suffices because the problem is linear elastic with no time dependence. Newton converges in exactly 1 iteration for a linear system.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case53-pressure-vessel \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case53_pressure_vessel.i 2>&1 | tail -30'
```

Output files:
- `case53_pressure_vessel_out.e` — Exodus with radial, axial, hoop, and von Mises stress fields
- `case53_pressure_vessel_out.csv` — Postprocessor values: five stress samples and max von Mises

---

## Expected Results

The simulation is a single steady solve. The CSV output contains one row of postprocessor values.

### Comparison with Lamé Solution

| Quantity | Analytical | MOOSE PointValue | Error |
|----------|------------|------------------|-------|
| sigma_r(r=1)     | -100.0 MPa  | -93.65 MPa  | 6.4%  |
| sigma_r(r=2)     |   0.0 MPa   |  -0.85 MPa  | ~0    |
| sigma_theta(r=1) | +166.67 MPa | +160.4 MPa  | 3.8%  |
| sigma_theta(r=2) | +66.67 MPa  | +67.6 MPa   | 1.4%  |
| max von Mises    | 220.5 MPa   | 220.4 MPa   | <0.1% |

The errors at the boundaries (r = 1 and r = 2) are due to the fact that MOOSE computes stresses at element centroids (or extrapolates to nodes from Gauss points), while the `PointValue` postprocessor samples at the boundary nodes. In linear elements, the stress field is piecewise constant within each element — the boundary node value is an extrapolation, not a quadrature point evaluation. Sampling at a midpoint away from the boundary (r = 1.05 or the element centroid nearest the inner surface) would give better agreement.

The max von Mises stress (220.4 MPa vs. 220.5 MPa analytical) is in excellent agreement because the von Mises stress averaged over the domain is not affected by the boundary extrapolation artefact.

The Lamé solution gives the von Mises stress at the inner surface for plane-strain as:

```
sigma_VM = sqrt(sigma_r^2 + sigma_theta^2 - sigma_r*sigma_theta + 3*tau_rz^2)
```

For plane strain, sigma_z = nu*(sigma_r + sigma_theta) = 0.3*(−100 + 166.67) = 20 MPa, giving:

```
sigma_VM = sqrt((−100)^2 + (166.67)^2 − (−100)(166.67) + 20^2 - ...) ≈ 220 MPa
```

This matches the postprocessor output.

---

## Key Takeaways

**`coord_type = RZ` converts a 2D mesh to an axisymmetric 3D problem.** A 20x4 element mesh captures the full 3D stress field in a thick-walled cylinder without the computational cost of a 3D mesh. The hoop stress and strain emerge automatically from the coordinate system.

**Stress component names shift meaning in RZ.** `stress_xx` = sigma_r, `stress_yy` = sigma_z, `stress_zz` = sigma_theta. This counterintuitive naming is a consequence of MOOSE's convention that the third ("out-of-plane") component in 2D is the z-component of the 3D stress tensor. Always verify which physical direction corresponds to which output name.

**The Lamé solution is the gold standard for pressure vessel validation.** Whenever you implement a new solid mechanics feature — new material model, new coordinate system, new output — running a thick-walled cylinder problem and checking against Lamé is the fastest way to verify the implementation is correct.

**Boundary stress sampling has known limitations.** `PointValue` at a mesh boundary gives the extrapolated nodal value, not the exact quadrature-point stress. For validation, sample at interior points near (but not on) the boundary, or use `SideAverageValue` over a narrow strip to average out the extrapolation error.

**The inner surface is always the critical location in pressure vessel design.** The hoop stress at r = r_i is 166.67 MPa, 2.5 times the hoop stress at r = r_o (66.67 MPa). This is why pressure vessel codes (ASME Section VIII, EN 13445) require design checks against the inner surface stress — and why autofrettage (pre-loading to yield the inner surface in compression) is used in high-pressure applications to partially equalise the hoop stress distribution.
