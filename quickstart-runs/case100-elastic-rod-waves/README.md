# Case 100: Elastic Wave on Thin Rod — Longitudinal Resonance with Sinusoidal Driving

## Overview

This case simulates longitudinal elastic wave propagation in a thin aluminium rod driven at near-resonance frequency. A sinusoidal displacement of 0.1 mm amplitude at approximately 250 Hz is applied to the left end; the right end is free (zero-stress Neumann condition). The rod dimensions are 10 m × 0.5 m, giving a quasi-1D geometry suitable for demonstrating standing wave formation and resonance buildup in the solid mechanics module.

The governing equation for axial displacement ξ is the classical wave equation:

```
d²xi/dt² = (E/rho) * d²xi/dx²
```

For aluminium (E = 70 GPa, ρ = 2700 kg/m³) the phase velocity is v_p = sqrt(E/ρ) ≈ 5092 m/s. The fundamental resonant frequency of a rod with one end driven and one end free is:

```
f1 = vp / (4L) = 5092 / 40 ≈ 127 Hz    (quarter-wave resonance, fixed-free)
f1 = vp / (2L) = 5092 / 20 ≈ 254 Hz    (half-wave resonance, driven-free)
```

The driving frequency of 250 Hz is close to the half-wave resonance at 254 Hz, causing the displacement amplitude at the free end to build up over successive reflections. This case is drawn from MIT 6.641 Lecture 16 (Elastodynamics and Wave Propagation) and employs the `Physics/SolidMechanics/Dynamic` action with Newmark-β time integration — the same framework as Case 20, now with a sinusoidal rather than impulsive source.

The simulation runs for 10 ms (approximately 2.5 wave transit times), long enough to observe the beginning of standing wave formation and the resonant amplitude amplification at the free end.

---

## The Physics

**Governing equation** (longitudinal wave on thin rod):

```
rho * d²xi/dt² = E * d²xi/dx²
```

In weak form with Newmark-β integration, the solid_mechanics Dynamic action assembles the full inertial + stiffness system automatically.

**Boundary conditions:**

| Boundary | Variable | Condition | Value |
|----------|----------|-----------|-------|
| Left (x = 0) | disp_x | Sinusoidal Dirichlet | 1e-4 sin(1570.8 t) m |
| Right (x = 10) | disp_x | Free (natural Neumann) | sigma_xx · n = 0 |
| All edges | disp_y | Dirichlet = 0 | Enforces 1D wave motion |

**Material properties:**

| Property | Symbol | Value | Units |
|----------|--------|-------|-------|
| Young's modulus | E | 70.0 × 10⁹ | Pa |
| Poisson's ratio | nu | 0 | — |
| Density | rho | 2700 | kg/m³ |
| Wave speed | v_p | sqrt(E/rho) ≈ 5092 | m/s |

Poisson's ratio is set to zero to enforce pure 1D longitudinal wave behaviour without transverse coupling.

**Domain geometry:**

- Rod: [0, 10] × [0, 0.5], 100 × 5 QUAD4 elements
- Element size: 0.1 m longitudinally = v_p × dt / 5 (well resolved)
- Time step: 2×10⁻⁵ s, giving ~50 steps per driving period (T = 4 ms at 250 Hz)

---

## Input File Walkthrough

### Mesh

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100
  ny   = 5
  xmin = 0.0
  xmax = 10.0
  ymin = 0.0
  ymax = 0.5
[]
```

One hundred elements along the rod axis provide 10 elements per wavelength at 250 Hz (wavelength = v_p / f ≈ 20 m; 10 m rod captures half a wavelength). Five elements in y are sufficient for a 2D representation of the essentially 1D problem.

### Physics Block (Dynamic Solid Mechanics)

```
[Physics]
  [SolidMechanics]
    [Dynamic]
      [all]
        add_variables  = true
        newmark_beta   = 0.25
        newmark_gamma  = 0.5
        strain         = SMALL
        density        = 2700
        generate_output = 'stress_xx vonmises_stress'
      []
    []
  []
[]
```

The `Dynamic` action automatically creates `disp_x` and `disp_y` variables, stress divergence kernels, and the Newmark-β inertial kernel. The parameters beta = 0.25 and gamma = 0.5 yield the unconditionally stable trapezoidal (Crank-Nicolson) time integrator, which introduces no numerical damping. `generate_output = 'stress_xx vonmises_stress'` adds auxiliary variables for the axial stress and von Mises stress, useful for visualising the stress wave pattern.

### Functions and Boundary Conditions

```
[Functions]
  [sinusoidal_drive]
    type       = ParsedFunction
    expression = '1.0e-4 * sin(1570.8 * t)'
  []
[]
```

The angular frequency 1570.8 rad/s corresponds to exactly 250 Hz. The amplitude 1×10⁻⁴ m (0.1 mm) is small relative to the rod length (1/100,000 of L), confirming the small-strain assumption.

```
[BCs]
  [drive_left]
    type     = FunctionDirichletBC
    variable = disp_x
    boundary = left
    function = sinusoidal_drive
  []
  [fix_y_bottom] ... [fix_y_right]
    type     = DirichletBC
    variable = disp_y
    value    = 0.0
  []
[]
```

The right end has no explicit BC on `disp_x`, so MOOSE applies the natural (zero-traction) Neumann condition — physically correct for a free end where the stress is zero.

### Materials

```
[Materials]
  [elasticity]
    type           = ComputeIsotropicElasticityTensor
    youngs_modulus  = 70.0e9
    poissons_ratio  = 0.0
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]
```

`ComputeIsotropicElasticityTensor` fills the 4th-order elasticity tensor C_ijkl from the two Lamé-equivalent parameters (E, ν). With ν = 0 the tensor reduces to the 1D elastic modulus E in the x-direction. `ComputeLinearElasticStress` uses the small-strain assumption (Cauchy stress = C : epsilon).

### Postprocessors and Outputs

Three `PointValue` postprocessors monitor axial displacement at x = 0 (driven end), x = 5 m (mid-rod), and x = 10 m (free end). These time-series values in the CSV output show the standing wave buildup and the phase relationship between the driving point and the free end. Exodus output is written every 5 time steps to keep file size manageable.

### Executioner

```
[Executioner]
  type     = Transient
  solve_type = PJFNK
  dt       = 2.0e-5
  end_time = 0.01
[]
```

PJFNK (Preconditioned Jacobian-Free Newton-Krylov) with LU factorisation. The time step 2×10⁻⁵ s gives a Courant number v_p × dt / dx = 5092 × 2e-5 / 0.1 ≈ 1.0, right at the advective Courant limit (although Newmark-β is unconditionally stable, keeping C ≤ 1 ensures good temporal accuracy).

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case100-elastic-rod-waves \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case100_elastic_rod_waves.i 2>&1 | tail -30'
```

---

## Expected Results

**Wave parameters:**

| Quantity | Formula | Value |
|----------|---------|-------|
| Phase velocity v_p | sqrt(E/rho) | 5092 m/s |
| Transit time (one-way) | L / v_p | 1.96 ms |
| Driving period T | 1 / f | 4.0 ms |
| Half-wave resonance | v_p / (2L) | 254.6 Hz |
| Detuning |f_drive - f_res| | 4.6 Hz |

**Time series at the free end (x = 10 m):**

The displacement at the free end starts at zero. After each transit time (~2 ms) a reflected wave arrives back at the driven end, then re-reflects and returns to the free end. Because 250 Hz is close to (but not exactly at) the 254.6 Hz resonance, the free-end amplitude grows but does not grow without bound over the 10 ms simulation window. By t = 10 ms approximately 5 transit times have elapsed, and the free-end amplitude should be roughly 2–4× the driven amplitude (amplitude ratio = 1 / |sin(pi L f / v_p)|, which is large near resonance).

**Standing wave pattern:**

In the Exodus output, plotting `disp_x` at any fixed time near the end of the simulation should show a spatial pattern close to sin(pi x / (2L)) (the mode 1 shape for a fixed-free resonator), with an antinode at the free end and a node at the driven end.

**Stress wave:**

`stress_xx` in the Exodus file shows the compressive and tensile wave fronts travelling along the rod. Peak stress amplitude is sigma_max = rho * v_p * omega * A ≈ 2700 × 5092 × 1571 × 1e-4 ≈ 2.2 GPa (transiently, before significant resonance buildup), well within the illustrative range even though it would exceed aluminium's actual yield stress for demonstration purposes.

---

## Key Takeaways

- **Physics/SolidMechanics/Dynamic action**: A single block replaces the manual assembly of displacement variables, stress-divergence kernels, and inertial kernels — this is the recommended approach for all dynamic solid mechanics problems in MOOSE.
- **Newmark-β (beta=0.25, gamma=0.5)**: The trapezoidal rule variant is unconditionally stable and second-order accurate in time with no numerical dissipation, ideal for undamped wave propagation.
- **Near-resonance driving**: Driving a rod at a frequency close to a natural frequency causes progressive amplitude growth over successive wave transits, demonstrating resonant energy storage in elastic structures.
- **Poisson's ratio = 0**: Setting ν = 0 decouples the longitudinal and transverse degrees of freedom, enforcing a pure 1D wave without lateral Poisson contraction effects — the standard thin-rod approximation.
- **Free-end Neumann condition**: No explicit BC on `disp_x` at the right boundary is the natural way to enforce a free end in MOOSE; the FEM residual automatically produces zero-traction at the boundary.
- **Connection to MIT 6.641 Lecture 16**: This is the canonical elastodynamics example from Zahn's course, illustrating how mechanical waves in solids obey the same mathematical structure as electromagnetic waves in transmission lines.
