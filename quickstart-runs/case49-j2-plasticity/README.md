# Case 49: J2 Plasticity — Uniaxial Tension with Isotropic Hardening

## Overview

Cases 14 through 21 used linear elasticity: stress is always proportional to strain, and the material returns to its original shape when the load is removed. Real engineering materials do not behave this way beyond a certain stress level. When the applied stress exceeds the yield stress, metals undergo permanent plastic deformation — the material deforms irreversibly, and the stress-strain relationship becomes nonlinear.

This case implements J2 (von Mises) plasticity with isotropic hardening on a 2D plane-strain bar pulled in tension. The yield criterion is based on the second invariant of the deviatoric stress tensor (J2), which for a uniaxial load reduces to the familiar condition that yielding occurs when the axial stress reaches sigma_y. After yielding, the flow stress increases linearly with accumulated plastic strain through an isotropic hardening modulus H.

The result is the classical bilinear stress-strain curve: a steep elastic slope (E = 200 GPa) followed by a shallower hardening slope (H = 1000 MPa) after yield. MOOSE's radial return algorithm enforces this consistently at every integration point on every time step.

New concepts introduced in this case:

- **`IsotropicPlasticityStressUpdate`**: implements J2 plasticity via the radial return mapping algorithm, the standard method for enforcing plastic consistency at the quadrature-point level.
- **`ComputeMultipleInelasticStress`**: the host stress object that orchestrates one or more inelastic models and assembles the final stress from elastic and inelastic contributions.
- **Incremental strain**: plasticity requires the strain to be computed incrementally (step-by-step) rather than from total displacement, because plastic strain history must be tracked.
- **`PiecewiseLinear` hardening function**: the yield stress as a function of accumulated plastic strain, defining the bilinear response analytically.

---

## The Physics

### The Physical Problem in Plain English

A 1 mm wide by 5 mm tall steel bar (modelled as a 2D plane-strain slice) is gripped at the bottom and pulled upward at the top. The top displacement increases linearly with time: at t = 1, the top has moved 0.05 mm — a nominal axial strain of 0.05/5 = 1%.

Steel yields at 250 MPa. Before yield, stress increases steeply with strain at the elastic modulus E = 200 GPa. Once the stress reaches 250 MPa, further deformation generates plastic strain with a much gentler slope governed by the hardening modulus H = 1000 MPa. The total strain splits as:

```
epsilon_total = epsilon_elastic + epsilon_plastic
sigma = E * epsilon_elastic         (elastic law)
sigma = sigma_y + H * epsilon_plastic  (hardening law after yield)
```

### Governing Equations

**Equilibrium** (quasi-static, no body forces):

```
div(sigma) = 0    in [0,1] x [0,5]
```

**Yield surface** (J2 / von Mises):

```
f(sigma) = sqrt(3/2 * s:s) - (sigma_y + H * epsilon_p) = 0
```

where s is the deviatoric stress tensor and epsilon_p is the accumulated equivalent plastic strain.

**Bilinear analytical solution** for uniaxial tension:

```
Elastic regime  (epsilon < epsilon_y):
  sigma = E * epsilon
  epsilon_y = sigma_y / E = 250 / 200000 = 0.00125

Plastic regime  (epsilon > epsilon_y):
  sigma = sigma_y + H * (epsilon - epsilon_y)
        = 250 + 1000 * (epsilon - 0.00125)    (MPa)
```

**Hardening function** (PiecewiseLinear):

```
epsilon_p:    0.0    0.10
sigma_y:    250.0   350.0   MPa
=> H = (350 - 250) / 0.10 = 1000 MPa
```

### Domain Diagram

```
         top: disp_y = 0.05*t  (tension, t goes 0→1)
         ______________________________
        |                              |
        |  Steel bar                   |  disp_x = 0 (left, symmetry)
        |  E  = 200 GPa                |
        |  nu = 0.3                    |
        |  sigma_y = 250 MPa           |
        |  H = 1000 MPa                |
        |                              |
        |______________________________|
         bottom: disp_x = disp_y = 0  (fixed)

Mesh: 5 x 25 elements, [0,1] x [0,5] mm, plane strain
```

---

## Input File Walkthrough

The input file is `case49_j2_plasticity.i`.

### `[GlobalParams]`

```
displacements = 'disp_x disp_y'
```

Declares the displacement variable names once so that the QuasiStatic action and all BCs pick them up automatically.

### `[Functions]`

Two functions drive the simulation:

```
[top_pull]
  type = ParsedFunction
  expression = '0.05 * t'   # 0.05 mm max displacement at t=1
[]
[hardening_func]
  type = PiecewiseLinear
  x = '0.0    0.10'    # plastic strain
  y = '250.0  350.0'   # yield stress (MPa)
[]
```

`top_pull` provides the displacement-controlled loading: the top edge moves at a constant rate of 0.05 mm per unit time. `hardening_func` provides the yield stress as a function of accumulated plastic strain to `IsotropicPlasticityStressUpdate`. The slope between the two points is H = (350 - 250) / 0.10 = 1000 MPa.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = SMALL
  incremental = true
  add_variables = true
  generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_yy elastic_strain_yy'
[]
```

`incremental = true` is the critical setting for plasticity. The strain is computed from the difference between the current and previous displacement gradient, not from the displacement gradient alone. This is required for the return-mapping algorithm to correctly track accumulated plastic strain history.

`generate_output` creates AuxVariables and AuxKernels for five field quantities: axial stress, lateral stress, von Mises stress, plastic strain component, and elastic strain component. These are reported by postprocessors and written to the Exodus file.

### `[BCs]`

| Name | Variable | Boundary | Value | Purpose |
|------|----------|----------|-------|---------|
| `fix_bottom_y` | disp_y | bottom | 0 | Prevents rigid-body translation |
| `fix_bottom_x` | disp_x | bottom | 0 | Prevents rigid-body rotation |
| `pull_top_y` | disp_y | top | `top_pull` function | Displacement-controlled loading |
| `fix_left_x` | disp_x | left | 0 | Symmetry (plane-strain bar) |

The left edge symmetry condition (`disp_x = 0`) models the interior of a wider bar. The top edge is displacement-controlled rather than load-controlled because load control becomes ill-conditioned near yield (the tangent stiffness drops sharply), while displacement control remains well-posed throughout.

### `[Materials]`

Three material objects form the plasticity model:

**`ComputeIsotropicElasticityTensor`**: builds the 4th-order stiffness tensor from E = 200e3 MPa and nu = 0.3. Note units are MPa throughout this input file (consistent with mm and kN).

**`IsotropicPlasticityStressUpdate`**: the J2 return-mapping algorithm. At each quadrature point and time step, it:
1. Computes the trial stress assuming the entire strain increment is elastic.
2. Checks whether the trial stress violates the yield surface.
3. If yielded, solves a 1D nonlinear equation (the consistency condition) to find the plastic strain increment that returns the stress to the yield surface.
4. Reports the final stress, elastic strain, and plastic strain.

The `hardening_function` parameter accepts the PiecewiseLinear function object, which provides sigma_y(epsilon_p) during the consistency iteration.

**`ComputeMultipleInelasticStress`**: coordinates the inelastic models (here just `plasticity`) with the elasticity tensor. The `tangent_operator = elastic` option uses the elastic (not algorithmic) tangent, which is simpler and sufficient for convergence with PJFNK.

### `[Postprocessors]`

| Name | Type | Variable | Meaning |
|------|------|----------|---------|
| `avg_stress_yy` | ElementAverageValue | stress_yy | Average axial stress (tracks bilinear curve) |
| `max_vonmises` | ElementExtremeValue (max) | vonmises_stress | Peak von Mises stress |
| `avg_plastic_strain_yy` | ElementAverageValue | plastic_strain_yy | Average plastic strain (zero until yield) |
| `avg_elastic_strain_yy` | ElementAverageValue | elastic_strain_yy | Average elastic strain (locks near sigma_y/E after yield) |
| `top_disp_y` | SideAverageValue | disp_y | Top face displacement (loading parameter) |

### `[Executioner]`

```
type = Transient
solve_type = 'PJFNK'
dt = 0.02
end_time = 1.0
```

Fifty time steps ramp the displacement from 0 to 0.05 mm. PJFNK (Preconditioned Jacobian-Free Newton-Krylov) is used because `tangent_operator = elastic` provides an approximate rather than exact Jacobian. LU factorization with MUMPS provides the preconditioner.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case49-j2-plasticity \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case49_j2_plasticity.i 2>&1 | tail -30'
```

Output files:
- `case49_j2_plasticity_out.e` — Exodus with stress, strain, and displacement fields
- `case49_j2_plasticity_out.csv` — Time series of postprocessor values

---

## Expected Results

The simulation runs 50 time steps. The key output is the relationship between `avg_stress_yy` and total strain (`top_disp_y / 5`), which should trace the bilinear elastoplastic curve.

### Elastic Phase (t = 0 to ~0.125)

Before yield, avg_stress_yy increases linearly with displacement:

```
sigma = E * (disp_y / L) = 200000 MPa * (disp_y / 5)
```

At t = 0.1 (disp_y = 0.005 mm, epsilon = 0.001), avg_stress_yy = 222 MPa — below yield.
At t = 0.12 (disp_y = 0.006 mm, epsilon = 0.0012), avg_stress_yy = 266 MPa — yield onset.

The plastic strain is identically zero throughout this phase.

### Plastic Phase (t > ~0.125)

After yield, avg_plastic_strain_yy grows while avg_elastic_strain_yy stabilises near sigma_y / E = 250 / 200000 = 0.00125:

| Time | avg_stress_yy (MPa) | avg_plastic_strain_yy | avg_elastic_strain_yy |
|------|---------------------|-----------------------|-----------------------|
| 0.14 | 283.0               | 1.34e-4               | 1.27e-3               |
| 0.50 | 295.1               | 3.75e-3               | 1.25e-3               |
| 1.00 | 303.7               | 8.72e-3               | 1.28e-3               |

The elastic strain component stays approximately constant at 0.00125 (= sigma_y / E) while the plastic strain carries the additional deformation. The total strain at t = 1 is 0.05 / 5 = 0.01, so epsilon_p ≈ 0.01 - 0.00125 = 0.00875, matching the CSV.

The hardening slope in the plastic region is:

```
d(sigma) / d(epsilon_total) = H_effective
```

For uniaxial plane-strain tension, the effective modulus differs slightly from H alone because of the Poisson constraint, but the dominant behaviour follows the analytical H = 1000 MPa.

The max_vonmises postprocessor tracks slightly above avg_stress_yy because von Mises stress at element corners differs from the element average under the plane-strain constraint.

---

## Key Takeaways

**J2 plasticity is the standard metal yielding model.** The von Mises yield surface in stress space is a cylinder aligned with the hydrostatic axis. Loading inside the cylinder is elastic; reaching the surface triggers plastic flow in the direction normal to the surface (associative flow rule).

**The radial return algorithm is unconditionally stable.** Unlike explicit plasticity schemes, the radial return map solves an implicit nonlinear equation at each quadrature point to enforce the consistency condition exactly. This allows large time steps without numerical instability.

**`incremental = true` is mandatory for plasticity.** Without incremental strain, MOOSE computes strain from the total displacement gradient. Plastic models require the strain increment at each step to perform the return mapping correctly.

**`ComputeMultipleInelasticStress` is the entry point for inelastic models.** Whether you have plasticity, creep, damage, or combinations, this object orchestrates the inelastic updates. Adding a second model (e.g., creep) simply means adding another entry to `inelastic_models`.

**Displacement control is preferred over load control for softening problems.** Near yield, the tangent stiffness drops sharply. Load-controlled problems become singular at the limit load; displacement-controlled problems remain well-posed and can trace post-peak behaviour.
