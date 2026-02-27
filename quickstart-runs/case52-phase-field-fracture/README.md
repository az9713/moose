# Case 52: Phase-Field Fracture — Notched Specimen Under Tension

## Overview

Classical fracture mechanics uses discrete crack representations: a sharp crack tip with a stress singularity, governed by stress intensity factors (K_I, K_II, K_III) or the energy release rate G. These approaches work well for simple geometries with single, pre-defined crack paths but struggle with crack branching, merging, and nucleation in complex domains — the crack geometry must be tracked explicitly, which becomes intractable in three dimensions.

Phase-field fracture avoids discrete crack tracking entirely. The sharp crack is replaced by a smooth damage field c that varies continuously from 0 (intact material) to 1 (fully broken material) over a length scale l. The damage field is governed by a variational principle: the fracture evolves to minimize the total energy, which is the sum of elastic strain energy and fracture surface energy. Cracks nucleate, propagate, branch, and merge automatically wherever energetically favorable.

This case simulates a 2D rectangular domain with a pre-existing edge notch under vertical tension. The notch creates a stress concentration that drives crack propagation across the specimen. The load-displacement curve shows a linear elastic rise, a peak at crack onset, and then softening as the crack propagates and the specimen loses load-carrying capacity.

New concepts introduced in this case:

- **`[Modules/PhaseField/Nonconserved]` action**: sets up the Allen-Cahn equation for the damage variable c, including TimeDerivative, ACBulk (bulk driving force), and ACInterface (gradient penalty) kernels.
- **`ComputeLinearElasticPFFractureStress`**: computes the degraded elastic stress `sigma = g(c) * C : epsilon` and provides the elastic energy density split needed for the damage driving force.
- **Spectral decomposition**: splits the elastic energy into tensile and compressive parts so that damage is driven only by tensile strain — compressive regions do not crack.
- **`DerivativeSumMaterial`**: assembles the total free energy F from separate elastic and fracture energy contributions, providing the derivatives needed by the Allen-Cahn kernels.
- **`PhaseFieldFractureMechanicsOffDiag`**: adds the off-diagonal Jacobian coupling between displacement and damage, enabling the fully coupled Newton solve.

---

## The Physics

### The Physical Problem in Plain English

A 1 mm wide by 0.5 mm tall specimen has an edge notch that runs from the left edge (x = 0) to the centre (x = 0.5) along the bottom. The bottom right half (the "noncrack" boundary from x = 0.5 to x = 1) is fixed in y. The top is pulled upward at a rate of 1e-3 mm per time unit. The damage field c starts at zero everywhere and evolves as the elastic energy builds up.

The notch creates a stress concentration at its tip (x = 0.5, y = 0). As the displacement increases, the tensile stress at the notch tip eventually provides enough driving force to overcome the fracture resistance gc. The damage variable grows from 0 toward 1, and the material stiffness in the damaged region degrades as g(c) = (1-c)^2. The crack propagates rightward and upward across the specimen, and the reaction force at the top boundary peaks and then falls as the effective cross-section is reduced by the damage band.

### Governing Equations

**Degraded elastic energy**:

```
psi_e(epsilon, c) = g(c) * psi_e^+(epsilon) + psi_e^-(epsilon)
g(c) = (1-c)^2 + eta    (quadratic degradation, eta = 1e-6 for numerical stability)
```

The spectral decomposition splits the elastic energy into positive-eigenvalue (tensile) part psi_e^+ and negative-eigenvalue (compressive) part psi_e^-. Only psi_e^+ is degraded, so compression does not drive damage.

**Allen-Cahn damage evolution**:

```
(1/M) * dc/dt = -dF/dc + div(kappa * grad(c))
```

where the total free energy density is:

```
F = psi_e + psi_fracture
psi_fracture = gc / (2*l) * c^2    (local fracture energy)
```

The material parameter combinations:

```
Mobility:  M = L = 1 / (gc * visco)      with visco = 1e-4
Gradient:  kappa = gc * l
```

This specific form of M and kappa ensures that the stationary damage profile across an isolated crack is an exponential of width l, and the total energy release per unit crack area is exactly gc.

**Fracture parameters**:

```
gc = 2.7e-3 kN/mm   (critical energy release rate)
l  = 0.04 mm         (regularization length; crack width ~ 4*l = 0.16 mm)
E  = 210 GPa (approximately, from lambda=80, mu=120 in symmetric_isotropic fill)
nu = 0.3
```

### Domain Diagram

```
         top: disp_y = 1e-3 * t  (slow tension, disp_x = 0)
         ____________________________
        |                            |
        |  2D specimen               |
        |  nx = 40, ny = 20          |
 notch  |  E ~ 210 GPa, nu = 0.3    |
=======+|  gc = 2.7e-3, l = 0.04    |
  x=0   x=0.5                      x=1
                noncrack boundary fixed in y
                (x from 0.5 to 1.0)

Notch: c = 0 imposed on left half of bottom (no BC; crack nucleates from tip)
Fixed: disp_y = 0 on right half of bottom (noncrack nodeSet)
```

---

## Input File Walkthrough

The input file is `case52_phase_field_fracture.i`.

### Mesh Generation

```
[gen] -- GeneratedMeshGenerator: 40x20 elements, [0,1] x [0,0.5]
[noncrack] -- BoundingBoxNodeSetGenerator
  new_boundary = noncrack
  bottom_left = '0.5 0 0'
  top_right   = '1 0 0'
```

`BoundingBoxNodeSetGenerator` creates a named node set for the right half of the bottom boundary. This is used to apply the `disp_y = 0` BC only to the supported region, leaving the notch region (left half of bottom) free — simulating a notch that is open and unloaded.

### `[Modules/PhaseField/Nonconserved]`

```
[c]
  free_energy = F
  kappa       = kappa_op
  mobility    = L
[]
```

This action automatically creates three kernels for the damage variable c:
- `TimeDerivative`: the (1/M) * dc/dt term
- `ACBulk`: the driving force term -dF/dc (derivative of the local energy with respect to c)
- `ACInterface`: the gradient penalty term div(kappa * grad(c)) (regularization)

The `free_energy = F` parameter points to the `DerivativeSumMaterial` named `F`, which provides the first and second derivatives of the total energy with respect to c for the ACBulk kernel.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[mech]
  add_variables = true
  strain = SMALL
  additional_generate_output = 'stress_yy vonmises_stress'
  save_in = 'resid_x resid_y'
[]
```

`save_in = 'resid_x resid_y'` saves the mechanical residual force vector to AuxVariables. The reaction force postprocessor `NodalSum` on `resid_y` at the top boundary gives the total upward force on the loading surface — the y-axis of the load-displacement curve.

### `[Kernels]` — Off-Diagonal Coupling

```
[solid_x]
  type = PhaseFieldFractureMechanicsOffDiag
  variable = disp_x
  component = 0
  c = c
[]
[solid_y]
  type = PhaseFieldFractureMechanicsOffDiag
  variable = disp_y
  component = 1
  c = c
[]
```

These kernels add the off-diagonal Jacobian terms d(R_mech)/dc. Without them, Newton would use a block-diagonal preconditioner that ignores the mechanical-damage coupling, which would converge poorly or not at all near the peak load where the coupling is strongest.

### `[Materials]`

The material system has six objects working together:

**`GenericConstantMaterial`** (`pf_constants`): declares `gc_prop = 2.7e-3`, `l = 0.04`, and `visco = 1e-4` as material properties available to other materials by name.

**`ParsedMaterial`** (`define_mobility`): computes `L = 1 / (gc_prop * visco)`.

**`ParsedMaterial`** (`define_kappa`): computes `kappa_op = gc_prop * l`.

**`ComputeElasticityTensor`** (`elasticity_tensor`): builds the stiffness tensor from Lame parameters lambda = 80 GPa and mu = 120 GPa using `fill_method = symmetric_isotropic`. This gives approximately E = 210 GPa and nu = 0.3.

**`ComputeLinearElasticPFFractureStress`** (`damage_stress`): the core fracture-mechanics stress object. It:
- Computes the spectral decomposition of the strain tensor
- Computes psi_e^+ (tensile) and psi_e^- (compressive) energy parts
- Multiplies psi_e^+ by g(c) (degradation) to get the damaged stress
- Deposits `elastic_energy`, `degradation`, and `local_fracture_energy` as material properties for use by the free energy assembly

**`DerivativeParsedMaterial`** (`degradation`): implements g(c) = (1-c)^2*(1-eta) + eta. The `derivative_order = 2` instructs the object to compute and store dc and d2c derivatives automatically, which the ACBulk kernel needs for the Newton Jacobian. `disable_fpoptimizer = true` and `enable_jit = false` are workarounds for known issues with JIT compilation in certain Docker environments.

**`DerivativeParsedMaterial`** (`local_fracture_energy`): implements psi_fracture = gc * c^2 / (2*l). Again provides c and cc derivatives automatically.

**`DerivativeSumMaterial`** (`fracture_driving_energy`): sums `elastic_energy` and `local_fracture_energy` into a single material `F` with consistent derivatives, as required by the `[PhaseField/Nonconserved]` action.

### `[Preconditioning]`

```
[smp]
  type = SMP
  full = true
[]
```

Single Matrix Preconditioner with `full = true` includes all off-diagonal blocks in the preconditioner. This is required for the fully coupled mechanics-damage system to converge with PJFNK — without it, the Newton iterations diverge near the crack tip where the off-diagonal coupling is strong.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case52-phase-field-fracture \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case52_phase_field_fracture.i 2>&1 | tail -30'
```

Output files:
- `case52_phase_field_fracture_out.e` — Exodus with damage field c, stress, and displacement
- `case52_phase_field_fracture_out.csv` — Time series of reaction force, max damage, and top displacement

---

## Expected Results

The simulation runs 10 time steps from t = 0 to t = 10, applying a total displacement of 0.01 mm to the top boundary.

### Load-Displacement Curve

The reaction force (reaction_force_y in the CSV) first rises linearly, peaks, and then softens as the crack propagates:

| Time | top_disp_y (mm) | max_damage | reaction_force_y (kN/mm) |
|------|-----------------|------------|---------------------------|
| 1    | 0.001           | 0.000      | 0.284                     |
| 3    | 0.003           | 0.143      | 0.769                     |
| 5    | 0.005           | 0.544      | 0.942                     |
| 6    | 0.006           | 0.799      | 0.876                     |
| 8    | 0.008           | ~1.0       | 0.559                     |
| 10   | 0.010           | ~1.0       | 0.316                     |

The reaction force peaks at approximately 0.94 kN/mm around t = 5, when max_damage ≈ 0.54. After that, the crack band is fully formed and the load-carrying capacity decreases as the damaged region widens. By t = 8, max_damage exceeds 1 (a known quirk of the Allen-Cahn formulation without a constraint — the damage slightly overshoots 1 near the crack tip before settling).

### Damage Field Topology

In the Exodus file, the damage field c shows:
- A narrow band of c → 1 propagating from the notch tip (x = 0.5, y = 0) upward and rightward
- The undamaged region (c ≈ 0) on the left side of the crack and in the upper right of the specimen
- A smooth transition zone of width ~4l = 0.16 mm between intact and damaged material

The crack does not propagate into the compressed region below the notch because spectral decomposition prevents compressive energy from driving damage.

---

## Key Takeaways

**Phase-field fracture avoids explicit crack tracking.** The damage variable c is just another nodal DOF on the existing mesh. Crack initiation, propagation, branching, and arrest are all automatic consequences of the energy minimisation — no remeshing, no level-set functions, no crack-tip enrichment.

**The regularization length l is a physical parameter, not just a numerical one.** The correct l for a material is related to the material's intrinsic length scale: l ≈ (E * Gc) / sigma_c^2 where sigma_c is the tensile strength. Choosing l too large gives artificially diffuse cracks; choosing l too small requires mesh refinement (the mesh must resolve l, typically requiring 4-5 elements across the damage band width of ~4l).

**Spectral decomposition is essential for realistic crack patterns.** Without it, compressive regions also accumulate damage, producing unphysical behaviour (compression-driven fracture). The `decomposition_type = strain_spectral` option ensures that only the positive (tensile) strain energy drives the damage evolution.

**`DerivativeParsedMaterial` with `derivative_order = 2` is the standard way to define free energies in MOOSE.** The ACBulk kernel needs first and second derivatives of F with respect to c; the `DerivativeParsedMaterial` computes these symbolically from the expression string. This avoids hand-coding derivative materials.

**The SMP full preconditioner is mandatory for coupled mechanics-damage problems.** The off-diagonal coupling terms between displacement and damage become dominant near the crack tip. Without `full = true`, the Krylov solver sees an incomplete Jacobian and stagnates.

**Phase-field fracture scales naturally to 3D.** The same input structure — replace 2D mesh with 3D, replace 2D BCs with 3D BCs — produces full 3D crack patterns without algorithmic changes. This is the primary advantage over LEFM-based approaches for complex geometries.
