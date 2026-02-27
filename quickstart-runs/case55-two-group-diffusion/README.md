# Case 55: 2-Group Neutron Diffusion — Fast/Thermal with Fission

## Overview

This case extends the one-group criticality problem of Case 54 to a two-energy-group model, which is the minimum necessary to capture the physics of a thermal (water-moderated) reactor. In a thermal reactor, fast neutrons born from fission must slow down (moderate) to thermal energies before they are likely to cause further fissions. The two-group model tracks fast neutrons (group 1) and thermal neutrons (group 2) separately, with a downscattering transfer term coupling them.

The result is a coupled 2x2 eigenvalue system: both flux variables phi1 (fast) and phi2 (thermal) are solved simultaneously, and the k_eff eigenvalue reflects the combined fission economy of the two groups. For the 40 cm slab with these cross sections, the system is supercritical (k_eff > 1), meaning a critical reactor would need to be smaller or would need control rods to reduce reactivity.

This case introduces `CoupledForce` as the kernel for inter-variable coupling. The downscatter source in the thermal equation (Sigma_s12 * phi1) and the fission source from thermal fissions into the fast group (nu_Sigma_f2 * phi2 -> phi1) are both handled as `CoupledForce` kernels operating across the two flux variables. Understanding the sign conventions for `CoupledForce` versus `Reaction` is critical to setting up multi-group diffusion correctly.

## The Physics

### Governing Equations

The two-group diffusion equations are:

**Group 1 (fast):**
```
-D1 * laplacian(phi1) + (Sigma_a1 + Sigma_s12) * phi1 = (1/k) * [nu_Sf1 * phi1 + nu_Sf2 * phi2]
```

**Group 2 (thermal):**
```
-D2 * laplacian(phi2) + Sigma_a2 * phi2 = Sigma_s12 * phi1
```

In words:
- Fast neutrons (phi1) are produced by fission from both groups and lost by absorption and downscattering to thermal energies.
- Thermal neutrons (phi2) are produced by downscattering from the fast group and lost only by absorption (thermal fission is included in the absorption cross section as a partial contribution).
- The removal cross section for the fast group is Sigma_r1 = Sigma_a1 + Sigma_s12 = 0.01 + 0.02 = 0.03 /cm.
- All fission neutrons are born fast (chi1 = 1, chi2 = 0), so the fission source appears only in the group 1 equation.

### Energy Group Structure

| Parameter | Value | Units | Physical meaning |
|-----------|-------|-------|-----------------|
| D1 | 1.5 | cm | Fast group diffusion coefficient |
| D2 | 0.4 | cm | Thermal group diffusion coefficient |
| Sigma_a1 | 0.01 | /cm | Fast absorption cross section |
| Sigma_a2 | 0.08 | /cm | Thermal absorption cross section |
| Sigma_s12 | 0.02 | /cm | Fast-to-thermal downscatter cross section |
| Sigma_r1 | 0.03 | /cm | Fast removal (= Sigma_a1 + Sigma_s12) |
| nu_Sigma_f1 | 0.005 | /cm | Fast fission production (small) |
| nu_Sigma_f2 | 0.10 | /cm | Thermal fission production (dominant) |
| chi1 | 1.0 | — | All fission neutrons born fast |
| chi2 | 0.0 | — | No fission neutrons born thermal |

The dominance of thermal fissions (nu_Sigma_f2 >> nu_Sigma_f1) reflects the physics of a thermal reactor: the fast fission factor is a small correction, while most energy comes from thermal fissions.

### Boundary Conditions

Vacuum boundary conditions are applied to both flux variables at x = 0 and x = L:

```
phi1(0) = phi1(L) = 0    (fast flux, vacuum)
phi2(0) = phi2(L) = 0    (thermal flux, vacuum)
```

Only phi1 requires `EigenDirichletBC` (the fission source only drives phi1), but it is good practice to apply it consistently. The top and bottom surfaces (y-direction) use natural boundary conditions.

### Domain and Mesh

- Geometry: 1D slab, x in [0, 40] cm, modeled as 2D quasi-1D (80 x 2 elements)
- Element size: 0.5 cm in x, fine enough for the fast group (migration length M1 ~ sqrt(D1/Sigma_r1) ~ 7 cm)
- Expected k_eff ~ 1.34 (supercritical by a significant margin)

## Input File Walkthrough

### `[Variables]`

```
[phi1]   # fast group flux
[]
[phi2]   # thermal group flux
[]
```

Two coupled variables are declared. MOOSE assembles a single coupled nonlinear system with both variables simultaneously. The Jacobian includes cross-variable entries from the `CoupledForce` kernels.

### `[Kernels]`

The kernels for each group are:

**Group 1 kernels:**

`diff1` (`MatDiffusion`, variable = phi1): diffusion term -D1 * laplacian(phi1)

`removal1` (`Reaction`, rate = 0.03, variable = phi1): removal term +(Sigma_a1 + Sigma_s12) * phi1. Positive rate means this is a loss.

`fission11` (`Reaction`, rate = -0.005, variable = phi1, `extra_vector_tags = 'eigen'`): fast-fission source -(1/k) * nu_Sf1 * phi1. Negative rate because fission is a source (gain). Tagged eigen to go into the B-matrix.

`fission21` (`CoupledForce`, v = phi2, coef = 0.10, variable = phi1, `extra_vector_tags = 'eigen'`): thermal-fission source -(1/k) * nu_Sf2 * phi2 appearing in group 1 equation. `CoupledForce` residual contribution is -coef * v * test, so coef = +nu_Sf2 = 0.10 gives -nu_Sf2 * phi2 * test (a source, i.e., negative residual). Tagged eigen.

**Group 2 kernels:**

`diff2` (`MatDiffusion`, variable = phi2): diffusion term -D2 * laplacian(phi2)

`absorption2` (`Reaction`, rate = 0.08, variable = phi2): absorption term +Sigma_a2 * phi2 (loss)

`scatter12` (`CoupledForce`, v = phi1, coef = 0.02, variable = phi2): downscatter source -Sigma_s12 * phi1 appearing in group 2 equation. `CoupledForce` residual = -coef * v * test, so coef = +0.02 gives -0.02 * phi1 * test (a source). This is NOT tagged eigen — it is a fixed physical coupling, not scaled by 1/k.

### `[Materials]`

```
[group_diffusion]
  type = GenericConstantMaterial
  prop_names  = 'D1   D2'
  prop_values = '1.5  0.4'
[]
```

Both diffusion coefficients are constant throughout the homogeneous slab. For a heterogeneous reactor (fuel pins + moderator), these would be spatially varying or block-specific.

### `[BCs]`

Vacuum conditions for both groups. A subtlety: only phi1 has `EigenDirichletBC` applied in this input (since the fission source only directly drives phi1), but phi2 still has ordinary `DirichletBC` to enforce zero flux at the boundaries.

### `[Executioner]`

```
type = Eigenvalue
solve_type = PJFNK
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre    boomeramg'
```

Identical to Case 54. SLEPc finds the dominant eigenvalue of the coupled 2x2 system.

### `[VectorPostprocessors]`

```
[flux_profiles]
  type = LineValueSampler
  variable = 'phi1 phi2'
  start_point = '0 0.25 0'
  end_point = '${L} 0.25 0'
  num_points = 81
  sort_by = x
[]
```

Both flux variables are sampled simultaneously along the centerline. The output CSV will contain columns for x, phi1, and phi2, allowing comparison of fast and thermal flux shapes.

### `[Postprocessors]`

`PointValue` postprocessors sample both fluxes at the slab center (x = 20 cm). The ratio phi2_center / phi1_center is the thermal-to-fast flux ratio, which characterizes how well-moderated the reactor is.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case55-two-group-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case55_two_group_diffusion.i 2>&1 | tail -30'
```

## Expected Results

The solver converges to k_eff ~ 1.342. This means the 40 cm slab is about 34% supercritical — it would need to be made smaller (or loaded with control rods) to achieve criticality.

The flux profiles will show:
- **phi1 (fast flux)**: broad cosine-like distribution peaking at center, relatively flat due to the large fast diffusion length (D1 = 1.5 cm)
- **phi2 (thermal flux)**: also a cosine shape but proportionally larger than phi1 at every point, since most neutrons thermalize before being absorbed

At the center (x = 20 cm), expect approximately:
```
phi1_center ~ 1.0   (normalized by SLEPc)
phi2_center ~ 3-5   (thermal flux substantially larger than fast)
```

The thermal-to-fast flux ratio (phi2/phi1) reflects the dominance of thermal reactions in a water-moderated reactor. In a fast reactor, this ratio would be much less than 1.

Physical interpretation of k_eff = 1.342: the infinite medium multiplication is k_inf = nu_Sigma_f2 / Sigma_a2 ~ 0.10 / 0.08 ~ 1.25 (thermal group alone). The 40 cm slab loses relatively few neutrons to leakage because the slab is much larger than the migration length, hence the moderate supercriticality.

## Key Takeaways

- Two-group diffusion introduces inter-variable coupling: `CoupledForce` handles source terms from one group appearing in another group's equation.
- The sign convention for `CoupledForce` is critical: the residual contribution is -coef * v * test. To add a source (negative residual contribution), use a positive coef. This is counterintuitive but consistent with MOOSE's residual formulation.
- Fission sources are tagged `extra_vector_tags = 'eigen'`; downscatter sources are not — they are fixed physical couplings that appear in the A-matrix.
- The thermal fission cross section (nu_Sigma_f2 = 0.10) dominates the fast fission cross section (nu_Sigma_f1 = 0.005) by a factor of 20, correctly capturing that thermal reactors depend on slow neutrons for most of their energy.
- k_eff > 1 for a 40 cm slab means the reactor would need to be subcritical-sized or controlled; k_eff < 1 would mean no self-sustaining reaction is possible at any size.
- The two-group model is the simplest practical model for thermal reactors. Real design codes use 2-69 energy groups depending on accuracy requirements.
