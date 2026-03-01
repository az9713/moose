# Case 103: Kelvin Polarization Force — Dielectric Fluid Rise in a Capacitor

## Overview

When a parallel-plate capacitor is partially immersed in a dielectric liquid, the Kelvin polarization body force pulls the liquid up into the gap between the plates. The liquid rises until the electrostatic force is balanced by gravity, reaching an equilibrium height:

```
h = (eps - eps0) * V² / (2 * rho * g * a²)
```

where ε is the liquid permittivity, V is the applied voltage, ρ is the liquid density, g is gravitational acceleration, and a is the plate spacing. This phenomenon is governed by MIT 6.641 Lectures 12 and 15 on the Kelvin force density and Maxwell stress tensor.

This case focuses on the electrostatic part of the problem: computing the electric potential distribution in a capacitor cross-section with a horizontal dielectric/air interface at y = 1 m. The domain [0,1] × [0,2] represents a unit-width capacitor with the left plate at Phi = 1 V and the right plate grounded. The lower half (y < 1) contains a dielectric liquid with ε_r = 3; the upper half (y > 1) is air with ε_r = 1.

The key physics is that the normal component of D = εE must be continuous across the dielectric-air interface, while the tangential component of E is continuous. For this geometry (plates parallel to y, interface perpendicular to y), the electric field is horizontal (E_x = -dPhi/dx) and the interface condition D_y^(dielectric) = D_y^(air) is satisfied trivially because there is no y-variation of the potential in the bulk of either region.

MOOSE solves nabla·(ε nabla Phi) = 0 using a two-block mesh with different permittivity values in each block, connected by a `StitchMeshGenerator`. This is the standard pattern for multi-material electrostatics in MOOSE.

---

## The Physics

**Governing equation** (Gauss's law, no free charge):

```
nabla . (eps * nabla Phi) = 0
```

In weak form: integral( eps * nabla Phi . nabla psi ) dV = 0

**Boundary conditions:**

| Boundary | Name | Condition | Value |
|----------|------|-----------|-------|
| Left plate (x = 0) | plate_high | Dirichlet | Phi = 1 V |
| Right plate (x = 1) | plate_ground | Dirichlet | Phi = 0 V |
| Top (y = 2) | top_bc | Natural Neumann | d Phi/dn = 0 (fringe) |
| Bottom (y = 0) | bottom_bc | Natural Neumann | d Phi/dn = 0 (fringe) |

The natural Neumann condition at top and bottom represents the electric field being parallel to those boundaries (no fringe field exits through the top or bottom of the modelled region).

**Material properties:**

| Region | y range | Block | Permittivity eps_r |
|--------|---------|-------|-------------------|
| Dielectric liquid | 0 to 1 | 0 | 3.0 |
| Air | 1 to 2 | 1 | 1.0 |

**Domain geometry:**

- Full domain: [0, 1] × [0, 2], two stacked blocks joined at y = 1
- Dielectric block: 30 × 20 QUAD4 elements (y = 0 to 1)
- Air block: 30 × 20 QUAD4 elements (y = 1 to 2)
- Total: 1200 elements across two subdomains

**Kelvin force density (for reference):**

The polarization force density in the dielectric is:

```
f = (eps - eps0) / 2 * nabla(E . E) = (eps - eps0) / 2 * nabla(|E|²)
```

At the sharp dielectric/air interface this becomes a surface force per unit area:

```
T_y = (eps_liquid - eps_air) / 2 * E_x²
```

which is positive (upward) when eps_liquid > eps_air. This surface traction, integrated over the interface, must balance the weight of the raised fluid column.

---

## Input File Walkthrough

### Mesh Construction

```
[Mesh]
  [dielectric]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 20
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
  [dielectric_id]
    type = SubdomainIDGenerator
    input = dielectric
    subdomain_id = 0
  []
  [air]
    type = GeneratedMeshGenerator
    ...
    ymin = 1
    ymax = 2
    subdomain_id = (before rename)
  []
  [air_id]
    type = SubdomainIDGenerator
    input = air
    subdomain_id = 1
  []
  [combine]
    type = StitchMeshGenerator
    inputs = 'dielectric_id air_id'
    stitch_boundaries_pairs = 'top bottom'
  []
  [rename]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'left right top bottom'
    new_boundary = 'plate_high plate_ground top_bc bottom_bc'
  []
[]
```

Two `GeneratedMeshGenerator` blocks create the dielectric (block 0) and air (block 1) regions separately. `SubdomainIDGenerator` assigns subdomain IDs before stitching. `StitchMeshGenerator` connects the top boundary of the dielectric mesh to the bottom boundary of the air mesh, creating a conforming interface at y = 1. Finally `RenameBoundaryGenerator` replaces the auto-generated boundary names with descriptive ones for the BCs.

### Variables and Kernels

```
[Variables]
  [phi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [gauss_law]
    type        = ADMatDiffusion
    variable    = phi
    diffusivity = permittivity
  []
[]
```

`ADMatDiffusion` implements integral(kappa * nabla phi . nabla psi) dV where kappa is the material property `permittivity`. When `permittivity` takes different values on block 0 (3.0) and block 1 (1.0), MOOSE automatically evaluates the correct permittivity at each quadrature point using the block-wise material assignment. The permittivity discontinuity at y = 1 is handled naturally — no special interface kernel is needed.

### Materials

```
[Materials]
  [dielectric_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '${eps_r}'
    block       = 0
  []
  [air_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '1.0'
    block       = 1
  []
[]
```

Block-restricted materials are the standard MOOSE mechanism for multi-material problems. The `block` parameter restricts the material to subdomain 0 (dielectric) and subdomain 1 (air) respectively. Both use the AD variant to match `ADMatDiffusion`.

### Auxiliary Variable

```
[AuxVariables]
  [E_squared]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [E_sq_aux]
    type              = ParsedAux
    variable          = E_squared
    expression        = 'phi * phi'
    coupled_variables = 'phi'
  []
[]
```

`E_squared` stores phi² as a proxy for visualising the potential field. A more accurate gradient-based |E|² would require `VariableGradientComponent` auxiliary kernels (one for each component) combined in a separate ParsedAux. The current proxy is sufficient to visualise the potential distribution; the true electric field magnitude can be computed in ParaView using the Gradient filter applied to the phi field.

### Postprocessors

| Postprocessor | Location | Expected value |
|--------------|----------|----------------|
| `avg_phi_dielectric` | block 0 | ~0.5 V (linear potential average) |
| `avg_phi_air` | block 1 | ~0.5 V (linear potential average) |
| `phi_interface` | (0.5, 1.0) | 0.5 V (midpoint between plates) |
| `phi_dielectric_mid` | (0.5, 0.5) | 0.5 V (linear in x) |
| `phi_air_mid` | (0.5, 1.5) | 0.5 V (linear in x) |
| `max_phi` | domain max | 1.0 V (at plate_high) |

### Executioner

```
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]
```

A single steady-state Newton solve. The problem is linear (nabla·(ε nabla Phi) = 0 with constant ε per block), so Newton converges in one iteration. BoomerAMG is efficient for the SPD Laplacian system.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case103-dielectric-fluid-rise \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case103_dielectric_fluid_rise.i 2>&1 | tail -30'
```

---

## Expected Results

**Potential distribution:**

For infinite parallel plates with no y-variation, the potential would be Phi = V(1 - x/a) = 1 - x, linear in x, independent of y and of ε_r. This is because nabla·(ε nabla Phi) = ε nabla²Phi = ε d²Phi/dx² = 0 is satisfied by any linear function in x, regardless of ε. The permittivity cancels and does not affect Phi.

This means all the postprocessors return Phi ≈ 0.5 V at x = 0.5, consistent with a linear profile:

```
Phi(x, y) = 1 - x    everywhere in domain
E_x = -dPhi/dx = 1 V/m    (uniform horizontal field)
E_y = 0                    (no vertical field in bulk)
```

**Interface condition:**

At the dielectric/air interface (y = 1), continuity of the tangential E (horizontal) and normal D (vertical) components is automatically satisfied:
- E_x is continuous: both sides have E_x = 1 V/m
- D_y is continuous: D_y = 0 on both sides (no vertical field)

**Kelvin force:**

The horizontal electric field E_x = 1 V/m is uniform across the interface. The surface Kelvin force per unit area is:

```
T_y = (eps_liquid - eps_air) / 2 * E_x²
    = (3 - 1) * eps0 / 2 * 1²
    = eps0  [N/m²]
```

For V = 1 V and a = 1 m (normalised), the equilibrium rise height is:

```
h = (eps_r - 1) * eps0 * V² / (2 * rho * g * a²) = 2 * eps0 / (2 * rho * g)
```

With physical values (eps₀ = 8.85×10⁻¹² F/m, ρ = 1000 kg/m³, g = 9.81 m/s²) this is a very small height — real demonstrations use high voltages (kV range) to produce measurable rise.

**Visualisation:**

Opening `case103_dielectric_fluid_rise_out.e` in ParaView shows a uniform horizontal potential gradient. The phi contours are vertical lines at x = 0, 0.25, 0.5, 0.75, 1.0. No visual discontinuity appears at the y = 1 interface because the potential and electric field are continuous there for this geometry.

---

## Key Takeaways

- **StitchMeshGenerator for multi-material domains**: Two independently generated meshes can be joined at a shared boundary, creating a conforming multi-block mesh without manual meshing tools. This is the standard MOOSE workflow for layered or composite material problems.
- **Block-restricted materials**: Assigning different `ADGenericConstantMaterial` objects to different subdomains is the correct way to implement spatially varying material properties in MOOSE electrostatics problems.
- **ADMatDiffusion for nabla·(eps nabla Phi)**: A single kernel with a block-wise permittivity property handles the full multi-material Gauss's law; no special interface treatment is needed at the dielectric/air boundary.
- **Potential independence from permittivity**: For infinite parallel plates, the potential distribution Phi = 1 - x is independent of ε_r because the linear function satisfies nabla²Phi = 0 identically. The permittivity affects the electric energy density and surface forces but not the field distribution itself in this geometry.
- **Kelvin force mechanism**: The polarization force density f = (ε - ε₀)/2 * nabla(|E|²) acts at the dielectric/air interface as a surface traction, pulling the higher-ε liquid into the stronger-field region — a macroscopic consequence of dipole alignment in the dielectric.
- **Connection to MIT 6.641 Lectures 12 and 15**: This is Zahn's parallel-plate capacitor fluid-rise problem, a classic demonstration of electromagnetic body forces on dielectric materials. The equilibrium height formula provides a practical measurement technique for liquid permittivity.
