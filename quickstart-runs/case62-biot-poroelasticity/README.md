# Case 62: Biot Poroelasticity — Coupled Consolidation Column

## Overview

Terzaghi's consolidation theory (Case 59) treats the solid skeleton as rigid — only the fluid moves. Biot poroelasticity (1941) extends this to account for the deformation of the solid skeleton coupled to pore fluid flow. When a load is applied, the soil skeleton deforms and fluid flows; the two processes are coupled because skeleton compression drives fluid out, and fluid pressure provides effective stress support to the skeleton.

This case solves the classic 1D Biot consolidation problem with explicit coupling between displacements and pore pressure. Three primary variables are solved simultaneously: horizontal displacement `disp_x`, vertical displacement `disp_y`, and `porepressure`. Unlike Cases 59-61 which used the `PorousFlowBasicTHM` or `PorousFlowUnsaturated` Actions, this case writes the coupling kernels explicitly — giving direct visibility into how the mechanical and flow equations are coupled at the kernel level.

Key concepts demonstrated:

- Explicit Biot poroelastic coupling via `StressDivergenceTensors` + `PorousFlowEffectiveStressCoupling` + `PorousFlowMassVolumetricExpansion`
- `PorousFlowEffectiveStressCoupling` kernel: adds the `-alpha * grad(p)` term to the momentum equation
- `PorousFlowMassVolumetricExpansion` kernel: adds the `alpha * d(eps_v)/dt` volumetric strain coupling to the fluid mass balance
- Analytical Biot settlement verification: delta = sigma * H / M (oedometric modulus)
- `PorousFlowVolumetricStrain` material: tracks the volumetric strain tensor needed for coupling

---

## The Physics

### Biot Governing Equations

The coupled Biot system consists of two sets of equations:

**Mechanical equilibrium** (with effective stress):
```
div(sigma') - alpha * grad(p) = 0
```

where `sigma'` is the effective (drained) stress tensor and the `alpha * grad(p)` term is the pore pressure contribution to total stress. In weak form, `PorousFlowEffectiveStressCoupling` adds the pressure-gradient coupling to the displacement equations.

**Fluid mass balance** (with volumetric strain coupling):
```
S * d(p)/dt + alpha * d(eps_v)/dt + div(-k/mu * grad(p)) = 0
```

where `eps_v = trace(strain)` is the volumetric strain. The `alpha * d(eps_v)/dt` term means that when the skeleton compresses (`d(eps_v)/dt < 0`), fluid is squeezed out even if the pressure gradient is zero. This is handled by `PorousFlowMassVolumetricExpansion`.

### Material Parameters

| Parameter | Value | Units | Physical meaning |
|-----------|-------|-------|-----------------|
| E | 100 MPa | Pa | Young's modulus (soft soil) |
| nu | 0.25 | — | Poisson's ratio |
| alpha | 0.6 | — | Biot coefficient |
| k | 1e-14 m^2 | m^2 | Permeability (tight clay) |
| phi | 0.3 | — | Porosity |
| Applied load | 10 kPa | Pa | Surface compressive stress |

### Analytical Solutions

**Initial pore pressure** (undrained loading):
When the load is applied instantaneously, the undrained response gives an initial pore pressure:
```
p_initial = alpha * sigma_applied = 0.6 * 10000 = 6000 Pa
```
This matches the `initial_condition = 6000` in the input file.

**Final settlement** (fully drained):
The oedometric (1D constrained) modulus M is:
```
M = E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
  = 100e6 * 0.75 / (1.25 * 0.5)
  = 120 MPa
```

The total settlement at full consolidation:
```
delta = sigma * H / M = 10000 * 10 / 120e6 = 8.333e-4 m = 0.833 mm
```

This is the value that `disp_y_top` converges to at late time.

**Consolidation coefficient**:
```
c_v = k / (mu * (S + alpha^2 / M))
```

With S from the Biot modulus (~2e-10 /Pa) and alpha^2/M = 0.36/120e6 = 3e-9 /Pa, the fluid coupling dominates. c_v ~ k / (mu * alpha^2/M) = 1e-14 / (0.001 * 3e-9) ~ 3.3e-3 m^2/s.

Characteristic time: t_c = H^2 / c_v = 100 / 3.3e-3 ~ 30000 s. The simulation runs to 2.5e7 s (many time constants), so full consolidation is achieved.

### Domain and Boundary Conditions

```
y = 10 m  ----  TOP: p = 0 (drained) + NeumannBC disp_y = -10 kPa (load)
          |    |
          | Biot poroelastic column  |
          | E = 100 MPa, nu = 0.25  |
          | alpha = 0.6             |
          | k = 1e-14 m^2           |
          |    |
y = 0 m   ----  BOTTOM: disp_x = 0, disp_y = 0 (fully fixed)

Left/right: disp_x = 0 (roller — plane strain / oedometer condition)

Mesh: 3 x 30 elements, [0,1] x [0,10] m
```

The roller conditions on left and right sides enforce the 1D oedometer geometry: no lateral displacement is allowed, so the only deformation is vertical compression.

---

## Input File Walkthrough

The input file is `case62_biot_poroelasticity.i`.

### `[GlobalParams]`

```
PorousFlowDictator = dictator
displacements = 'disp_x disp_y'
biot_coefficient = 0.6
```

All three global parameters are declared here: the PorousFlow dictator name, the displacement variable names (used by solid mechanics objects), and the Biot coefficient (used by several PorousFlow materials and the coupling kernels). Declaring `biot_coefficient` globally means it does not need to be repeated in each material block.

### `[Variables]`

Three primary variables are solved simultaneously:
- `porepressure` — initial 6000 Pa (undrained response to 10 kPa load)
- `disp_x` — zero initial condition (no initial horizontal displacement)
- `disp_y` — zero initial condition

### `[Kernels]`

The kernel list explicitly shows the Biot coupling. Six kernels total:

**Mechanics kernels (two equations, one per component):**

`StressDivergenceTensors` (components 0 and 1): the standard solid mechanics weak form `integral(sigma : grad(test) dV)`. These handle the elastic response of the solid skeleton.

`PorousFlowEffectiveStressCoupling` (components 0 and 1): adds the `alpha * p * I` correction to make the total stress effective. In the weak form this contributes `-alpha * p * div(test)` to the momentum equation residual. This is the term that couples pressure to the mechanical response.

**Flow kernels (one equation):**

`PorousFlowMassTimeDerivative`: time derivative of fluid mass `d(phi * rho * S)/dt`. Standard storage term.

`PorousFlowAdvectiveFlux`: Darcy advective flux `div(-rho * k/mu * grad(p))`. Standard Darcy term.

`PorousFlowMassVolumetricExpansion`: the `alpha * d(eps_v)/dt` coupling term. This is what makes this a full Biot formulation rather than just pressurized Darcy flow — skeleton compression drives additional fluid out even without a pressure gradient.

### `[UserObjects]`

```
[dictator]
  type = PorousFlowDictator
  porous_flow_vars = 'porepressure'
  number_fluid_phases = 1
  number_fluid_components = 1
[]
```

When using explicit kernels rather than the `PorousFlowBasicTHM` Action, the `PorousFlowDictator` must be declared manually. It registers `porepressure` as the single PorousFlow variable. The displacement variables are not listed here — they are handled by the solid mechanics system, not PorousFlow.

### `[Materials]`

The materials include both solid mechanics objects and PorousFlow objects:

**Solid mechanics:**
- `ComputeIsotropicElasticityTensor`: 4th-order stiffness tensor from E = 100 MPa, nu = 0.25
- `ComputeSmallStrain`: small-strain kinematics
- `ComputeLinearElasticStress`: constitutive law sigma = C : epsilon

**PorousFlow:**
- `PorousFlowTemperature`: provides temperature to PorousFlow materials (uses default 293 K)
- `PorousFlow1PhaseFullySaturated`: sets S = 1 everywhere (fully saturated)
- `PorousFlowMassFraction`: tracks fluid component fractions (single component here)
- `PorousFlowSingleComponentFluid`: links the fluid properties object to the PorousFlow material hierarchy
- `PorousFlowPorosity` with `mechanical = true`, `fluid = true`: porosity changes with volumetric strain and pressure — essential for the Biot coupling
- `PorousFlowPermeabilityConst`: k = 1e-14 m^2 (tight clay)
- `PorousFlowRelativePermeabilityConst`: k_r = 1 (always fully saturated)
- `PorousFlowVolumetricStrain`: computes the volumetric strain tensor for `PorousFlowMassVolumetricExpansion`
- `PorousFlowEffectiveFluidPressure`: computes alpha * p for the effective stress correction
- `PorousFlowConstantBiotModulus`: storage from Biot coefficient, solid compliance, fluid bulk modulus

### `[BCs]`

| Name | Variable | Boundary | Value | Purpose |
|------|----------|----------|-------|---------|
| `fix_bottom_x` | disp_x | bottom | 0 | No horizontal movement at base |
| `fix_bottom_y` | disp_y | bottom | 0 | No vertical movement at base |
| `fix_left_x` | disp_x | left | 0 | Roller — lateral confinement |
| `fix_right_x` | disp_x | right | 0 | Roller — lateral confinement |
| `top_load` | disp_y | top | -1e4 (NeumannBC) | 10 kPa downward applied load |
| `top_drained` | porepressure | top | 0 | Drained drainage surface |

The `NeumannBC` on `disp_y` at the top applies a traction (force per area). The value -1e4 Pa is negative because the force points downward (compressive). This is how a distributed surface load is applied in MOOSE solid mechanics.

### `[Executioner]`

```
dt = 5e5         # 500000 s per step
end_time = 2.5e7 # ~290 days total
```

Fifty time steps of 5e5 s each. The very long time steps are chosen because the consolidation time scale (t_c ~ 30000 s) is much smaller than the time step — the system consolidates rapidly. After the first few steps, the pore pressure has fully dissipated and only the mechanical equilibrium remains. The simulation then coasts to 2.5e7 s to confirm the steady final settlement.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case62-biot-poroelasticity \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case62_biot_poroelasticity.i 2>&1 | tail -30'
```

Output files:
- `case62_biot_poroelasticity_out.e` — Exodus with displacement and pore pressure fields
- `case62_biot_poroelasticity_out.csv` — Time series of p_base, p_mid, disp_y_top, p_avg

---

## Expected Results

### Pore Pressure Dissipation

The pore pressure dissipates almost entirely within the first two time steps (by t = 1e6 s), reflecting the fast consolidation with c_v >> H^2 / end_time:

| Time (s) | p_base (Pa) | p_mid (Pa) | p_avg (Pa) |
|----------|------------|-----------|-----------|
| 0        | 6000.0     | 6000.0    | 6000.0    |
| 5e5      | 51.8       | 38.9      | 34.6      |
| 1e6      | 1.65       | 1.18      | 1.06      |
| 1.5e6+   | ~0.01      | ~0.008    | ~0.007    |

After t = 1.5e6 s, pore pressure is effectively zero — full consolidation has occurred.

### Settlement

The vertical displacement at the top converges to the analytically predicted value:

```
disp_y_top at t = 5e5 s:  -8.316e-4 m  (~99.8% of final settlement)
disp_y_top at t >= 1.5e6 s:  -8.3333e-4 m  (converged, matches analytical)
```

Analytical prediction: delta = sigma * H / M = 10000 * 10 / (120e6) = 8.333e-4 m = 0.833 mm. The simulation matches this to 6 significant figures, confirming the Biot formulation is correctly implemented.

### Physical Interpretation

The sequence of events:
1. t = 0: Load applied. Pore pressure instantly rises to 6 kPa (alpha * sigma). Displacement is zero (undrained skeleton cannot settle yet).
2. t = 0 to ~1e5 s: Fluid drains through top surface. Pore pressure drops. Skeleton begins to compress.
3. t > 1e6 s: Fully drained. Pore pressure ~ 0 everywhere. Full settlement achieved.

---

## Key Takeaways

- Biot poroelasticity couples solid deformation and fluid flow bidirectionally: loads change pore pressure (undrained response), and pore pressure changes drive settlement (drained response). Terzaghi's simpler theory gets the same settlement but misses the initial undrained response.
- `PorousFlowEffectiveStressCoupling` adds the `alpha * grad(p)` term to the momentum equation — this is the pressure-to-mechanics coupling direction.
- `PorousFlowMassVolumetricExpansion` adds the `alpha * d(eps_v)/dt` term to the fluid mass balance — this is the mechanics-to-pressure coupling direction. Without this term, skeleton compression would not drive fluid out.
- The initial pore pressure equals alpha * sigma_applied in the undrained limit; the Biot coefficient alpha (0 to 1) controls how much of the applied load is initially carried by the fluid. For alpha = 1 (very soft or saturated skeleton), the fluid carries the entire load initially.
- The analytical settlement delta = sigma * H / M_oed is the key verification target. The oedometric modulus M includes the Poisson constraint from the roller BCs; it is higher than the Young's modulus E for nu > 0.
- Using explicit kernels rather than an Action (`PorousFlowBasicTHM`) requires declaring the `PorousFlowDictator` manually and listing all coupling kernels individually. This is more verbose but provides complete transparency into the formulation.
- `PorousFlowVolumetricStrain` is a material object (not a kernel) that computes the trace of the strain tensor; it must be present whenever `PorousFlowMassVolumetricExpansion` is used.
