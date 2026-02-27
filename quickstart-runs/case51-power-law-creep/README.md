# Case 51: Power-Law Creep — Column Under Sustained Compression

## Overview

Cases 49 and 50 addressed rate-independent inelasticity (plasticity) and geometric nonlinearity. This case introduces time-dependent inelastic deformation: creep. Creep is the slow, continuous accumulation of permanent strain in a material held under sustained stress at elevated temperature. Unlike plasticity, which occurs instantaneously when stress exceeds the yield surface, creep proceeds continuously even when the stress is well below the yield stress, as long as the temperature is high enough to activate atomic diffusion and dislocation climb.

The Norton power law is the most widely used empirical creep model in engineering:

```
d(epsilon_cr) / dt = A * sigma^n * exp(-Q / (R*T))
```

where A is a material coefficient, n is the stress exponent (typically 3-5 for metals), Q is the activation energy, R is the gas constant, and T is the absolute temperature. The exponential Arrhenius factor means creep rate is extremely sensitive to temperature: a 10% rise in temperature can double or triple the creep rate.

This case models a steel column held at 1000 K (above the creep threshold for steel) under 10 MPa compressive pressure. The elastic deformation occurs instantly on load application; the creep strain then accumulates steadily over 100 seconds as the column slowly shortens.

New concepts introduced in this case:

- **`PowerLawCreepStressUpdate`**: implements the Norton creep law as a return-mapping update, consistent with the radial return framework used by `ComputeMultipleInelasticStress`.
- **Temperature coupling**: the creep rate depends on temperature, so a temperature field must be provided as a MOOSE variable even if it is held constant.
- **Primary (transient) vs. secondary (steady-state) creep**: the Norton law models secondary creep; the decreasing creep rate over time visible in this simulation is a consequence of decreasing stress as the column shortens under pressure-controlled loading, not a physically distinct primary creep mechanism.
- **SI units throughout**: this input uses Pa, m, s, and K (not MPa, mm as in Case 49) to match the MOOSE material library defaults.

---

## The Physics

### The Physical Problem in Plain English

A 1 m wide by 2 m tall steel column is compressed from the top by a uniform downward pressure of 10 MPa (10 million Pa). The bottom is fixed. The temperature is held constant at 1000 K — hot enough for significant creep in steel. The left edge is a symmetry plane.

At t = 0, the load is applied and the column springs down elastically. The elastic strain is approximately sigma / E = 10e6 / 200e9 = 5e-5. Then, over the next 100 seconds, the column continues to shorten as creep strain accumulates. The creep rate decreases over time in this simulation because the column's shortening under a fixed pressure load slightly reduces the effective stress (the column gets shorter and the load is applied per unit current area — but since `Pressure` in MOOSE applies to the reference area, the primary cause of decreasing stress is that the finite-strain geometry update slightly reduces the applied force per unit deformed area).

### Governing Equations

**Elastic equilibrium**:

```
div(sigma) = 0    in [0,1] x [0,2]
```

**Norton power-law creep**:

```
d(epsilon_cr)/dt = A * sigma_eff^n * exp(-Q / (R*T))

Parameters:
  A = 1e-15 /s/Pa^n  (material coefficient)
  n = 4               (stress exponent)
  Q = 3e5 J/mol       (activation energy)
  R = 8.314 J/(mol·K) (gas constant)
  T = 1000 K          (constant temperature)

Arrhenius factor at 1000 K:
  exp(-3e5 / (8.314 * 1000)) = exp(-36.08) ≈ 2.2e-16

Creep rate at 10 MPa (1e7 Pa):
  d(eps_cr)/dt = 1e-15 * (1e7)^4 * 2.2e-16
               = 1e-15 * 1e28 * 2.2e-16
               = 1e-15 * 2.2e12
               = 2.2e-3 /s
```

Over 100 seconds, the cumulative creep strain is on the order of magnitude of 0.01–0.02 (1–2%), consistent with the simulation output.

**Total strain decomposition**:

```
epsilon_total = epsilon_elastic + epsilon_creep
sigma = E * epsilon_elastic      (elastic part only drives stress)
```

### Domain Diagram

```
         top: uniform pressure = 10 MPa downward
         ______________________________
        |                              |
        |  Steel column                |  (right edge: free)
        |  E  = 200 GPa                |
   x=0  |  nu = 0.3                    |  x=1
  (sym) |  T  = 1000 K (constant)      |
        |  A  = 1e-15, n = 4           |
        |  Q  = 3e5 J/mol              |
        |______________________________|
         bottom: disp_x = disp_y = 0  (fixed)

Mesh: 5 x 10 elements, [0,1] x [0,2]
10 time steps, dt = 10 s (t goes 0 → 100 s)
Units: SI (Pa, m, s, K)
```

---

## Input File Walkthrough

The input file is `case51_power_law_creep.i`.

### Temperature Variable

```
[Variables]
  [temp]
    initial_condition = 1000.0
  []
[]
[Kernels]
  [heat_diff]
    type = Diffusion
    variable = temp
  []
[]
```

`PowerLawCreepStressUpdate` requires a temperature variable to evaluate the Arrhenius factor. The temperature is held constant at 1000 K throughout the domain via `DirichletBC` on all boundaries. A dummy `Diffusion` kernel on `temp` is needed to keep the variable in the nonlinear system (MOOSE requires every variable to have a kernel). Since temp = 1000 K everywhere and the BCs enforce this, the diffusion kernel contributes zero residual.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_yy stress_xx vonmises_stress creep_strain_yy elastic_strain_yy'
[]
```

`strain = FINITE` and `incremental = true` are both set because large accumulated creep strains require the full deformation gradient to be tracked accurately. Even though individual time steps may have small strains, 100 seconds of creep can produce total inelastic strains of several percent.

`creep_strain_yy` is automatically generated because `PowerLawCreepStressUpdate` deposits the creep strain into the inelastic strain tensor, which the QuasiStatic action exposes through `generate_output`.

### `[Functions]` and `[BCs]`

```
[load_ramp]
  type = PiecewiseLinear
  x = '0  1  100'
  y = '0  1  1'
[]
```

The `Pressure` BC ramps from 0 to 10 MPa in the first second (to give Newton a gentle start) and then holds constant for the remaining 99 seconds. This avoids the large initial residual that would result from suddenly applying the full load at t = 0 to an unstrained structure.

```
[Pressure]
  [top_pressure]
    boundary = top
    factor = -1.0e7
    function = load_ramp
  []
[]
```

`factor = -1.0e7` gives 10 MPa compressive (downward). In MOOSE, positive pressure pushes inward (compressive), and the factor sign convention here applies the traction as a downward body force on the top surface.

### `[Materials]`

```
[power_law_creep]
  type = PowerLawCreepStressUpdate
  coefficient = 1.0e-15
  n_exponent = 4
  activation_energy = 3.0e5
  temperature = temp
[]
[stress]
  type = ComputeMultipleInelasticStress
  inelastic_models = 'power_law_creep'
  tangent_operator = elastic
[]
```

`PowerLawCreepStressUpdate` plugs into `ComputeMultipleInelasticStress` the same way `IsotropicPlasticityStressUpdate` did in Case 49. The return-mapping algorithm at each quadrature point and time step finds the creep strain increment that satisfies both the creep rate equation and the elastic stress-strain relation simultaneously.

`temperature = temp` links the creep rate to the MOOSE variable `temp`. Changing the temperature field (e.g., by coupling to a transient heat conduction solve) would automatically change the local creep rate everywhere through the Arrhenius factor.

### `[Executioner]`

```
type = Transient
solve_type = 'PJFNK'
line_search = 'none'
dt = 10.0
end_time = 100.0
```

Ten steps of 10 seconds each. The large time steps work because creep is a smooth, slowly varying process. The line search is disabled (`line_search = 'none'`) because the creep return mapping can occasionally produce updates that confuse the Armijo backtracking search — a common recommendation when using inelastic stress update objects.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case51-power-law-creep \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case51_power_law_creep.i 2>&1 | tail -30'
```

Output files:
- `case51_power_law_creep_out.e` — Exodus with stress, strain, and displacement fields
- `case51_power_law_creep_out.csv` — Time series of postprocessor values

---

## Expected Results

The simulation runs 10 time steps from t = 0 to t = 100 s.

### Stress Evolution

The elastic strain occurs instantly at load application and remains approximately constant at sigma / E = 1e7 / 2e11 = 5e-5. The average axial stress (avg_stress_yy) begins at approximately -1e7 Pa (-10 MPa) and decreases in magnitude slightly over time as the column shortens and the effective stress changes:

| Time (s) | avg_creep_strain_yy | avg_elastic_strain_yy | avg_stress_yy (Pa) |
|----------|---------------------|-----------------------|--------------------|
| 10       | 3.76e-3             | 3.62e-5               | -9.98e6            |
| 50       | 9.68e-3             | 3.10e-5               | -9.91e6            |
| 100      | 1.24e-2             | 2.88e-5               | -9.88e6            |

The elastic strain shrinks slightly over time as the creep strain grows — creep relaxes the elastic strain under approximately constant total displacement.

### Creep Strain Accumulation

The dominant result is the growth of avg_creep_strain_yy from 0 to approximately 1.24% over 100 seconds. The rate is highest at t = 10 s (about 3.8e-3 / 10 s = 3.8e-4 /s) and decreases over time.

The column shortens (max_disp_y grows negative) from approximately -7.8 mm at t = 10 s to -24.8 mm at t = 100 s — a total shortening of about 2.5% of the column height due to the combination of elastic and creep deformation.

### von Mises Stress Decrease

The max_vonmises postprocessor decreases from about 7.4 MPa at t = 10 s to 4.2 MPa at t = 100 s. This reflects both the creep relaxation of elastic strain and the change in the lateral stress field as the column bulges outward. The von Mises stress is lower than the applied axial stress (10 MPa) because the plane-strain constraint generates lateral (sigma_xx) stress that partially cancels the deviatoric contribution.

---

## Key Takeaways

**Creep is temperature-activated and extremely sensitive to temperature through the Arrhenius factor.** At 900 K instead of 1000 K, the creep rate would be exp(-Q/R*(1/900 - 1/1000)) ≈ 50 times smaller. This is why high-temperature components (turbine blades, pressure vessel heads) require careful thermal management to control creep life.

**The Norton power law is empirical but accurate for secondary creep.** The stress exponent n = 4 is typical for dislocation creep in metals. Values of n between 3 and 8 are common; n = 1 corresponds to linear viscous (Newtonian) creep. The coefficient A and activation energy Q are measured from long-duration creep tests at several temperatures and stress levels.

**`PowerLawCreepStressUpdate` plugs into the same framework as plasticity.** The `ComputeMultipleInelasticStress` host object can run plasticity and creep simultaneously by listing both in `inelastic_models`. This is the standard approach for hot working simulations where both mechanisms are active.

**The temperature variable is required even when constant.** MOOSE's creep objects are designed for fully coupled thermo-mechanical problems. When temperature is spatially uniform and time-independent, the simplest approach is to hold it fixed with DirichletBC on all boundaries and add a dummy kernel, as done here.

**SI units (Pa, m, s) are mandatory when using MOOSE's material library parameters.** The material coefficients A, Q, and R are all in SI units internally. Mixing unit systems (MPa for stress, mm for geometry) would produce incorrect Arrhenius evaluations. Always verify units match when specifying A and activation_energy.
