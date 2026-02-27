# Case 61: Unsaturated Flow — Wetting Front in Soil Column

## Overview

When rain falls on dry soil, water infiltrates downward through the unsaturated (vadose) zone — the region between the surface and the water table where pore space is only partially filled with water. This process is governed by Richards' equation, which extends Darcy's law to partially saturated conditions by introducing two nonlinear constitutive relationships: the saturation-pressure relationship (the water retention curve) and the relative permeability function.

This case models a 2 m vertical soil column that is initially dry (suction pressure -50 kPa) and then saturated at the surface (ponded water, p = 0). A sharp wetting front advances downward through the column over 1 hour as gravity and capillary suction draw water downward. The problem demonstrates how MOOSE's `PorousFlowUnsaturated` Action handles the highly nonlinear Richards equation with van Genuchten-Mualem constitutive relations.

Key concepts demonstrated:

- `PorousFlowUnsaturated` Action for variably saturated (Richards') flow
- Van Genuchten saturation-pressure model (`van_genuchten_alpha`, `van_genuchten_m`)
- Corey relative permeability model (`relative_permeability_type = COREY`)
- Gravity-driven infiltration with capillary suction driving the wetting front
- Nonlinear convergence challenges in unsaturated flow

---

## The Physics

### Richards' Equation

For single-phase unsaturated flow, Richards' equation is:

```
d(phi * rho * S) / dt + div(rho * k * k_r / mu * (grad(p) - rho * g)) = 0
```

where:
- `phi = 0.35` — porosity
- `rho = 1000 kg/m^3` — fluid density
- `S(p)` — saturation as a function of pressure (van Genuchten)
- `k = 1e-11 m^2` — saturated intrinsic permeability
- `k_r(S)` — relative permeability (Corey model)
- `mu = 0.001 Pa.s` — dynamic viscosity
- `g = (0, -9.81, 0)` — gravitational acceleration

The key difference from fully saturated Darcy flow is the presence of S(p) and k_r(S): saturation and relative permeability both depend on pressure, making the equation strongly nonlinear.

### Van Genuchten-Mualem Retention Model

The van Genuchten model relates effective saturation S_eff to the capillary pressure head:

```
S_eff(p) = (1 + (alpha * |p|)^n)^(-m)    for p < 0
S_eff = 1                                  for p >= 0
```

where:
- `alpha = 3.5e-4 /Pa` — inverse of air-entry pressure (relates to pore size)
- `n` related to `m` by: `m = 1 - 1/n`, and the input uses `van_genuchten_m = 0.5` (corresponding to n = 2)
- At p = -50 kPa: `S_eff = (1 + (3.5e-4 * 50000)^2)^(-0.5) = (1 + 17.5^2)^(-0.5) ~ 0.057`

This means the initial condition represents a very dry soil with only about 6% effective saturation.

### Corey Relative Permeability

The Corey model gives:

```
k_r(S_eff) = S_eff^n    with n = 2 (exponent)
```

At the initial dry state (S_eff ~ 0.057): k_r ~ 0.003 — the soil is nearly impermeable in this state. As the wetting front passes and S_eff rises toward 1, k_r rises toward 1 and flow becomes much easier. This strong nonlinearity is what creates a sharp wetting front rather than a smooth diffusing pressure wave.

### Wetting Front Velocity

Behind the wetting front (fully saturated): k_r = 1 and flow is governed by hydraulic conductivity K = k * rho * g / mu = 1e-11 * 1000 * 9.81 / 0.001 = 9.81e-5 m/s. The frontal velocity is approximately:

```
v_front ~ K / phi = 9.81e-5 / 0.35 ~ 2.8e-4 m/s
```

Over 1 hour (3600 s), the front travels roughly 2.8e-4 * 3600 ~ 1.0 m, reaching approximately the middle of the 2 m column. The actual front is sharper and advances somewhat slower due to capillary effects at the front.

### Domain and Boundary Conditions

```
y = 2 m  ----  TOP: saturated (p = 0, ponded water)
          |    |
          | Initially dry soil  |
          | p_init = -50 kPa   |
          | Sandy soil:        |
          |   k = 1e-11 m^2    |
          |   phi = 0.35       |
          |    |
y = 0 m   ----  BOTTOM: free drainage (natural zero-flux)

Gravity: (0, -9.81, 0) acts downward
Mesh: 2 x 100 elements (quasi-1D in y direction)
```

The top boundary is held at p = 0 (fully saturated, representing ponded surface water). The bottom boundary uses the natural zero-flux Neumann condition, representing free drainage or an infinitely deep column.

---

## Input File Walkthrough

The input file is `case61_unsaturated_flow.i`.

### `[PorousFlowUnsaturated]` Action

```
[PorousFlowUnsaturated]
  porepressure    = porepressure
  coupling_type   = Hydro
  gravity         = '0 -9.81 0'
  fp              = water
  relative_permeability_type = COREY
  relative_permeability_exponent = 2
  van_genuchten_alpha = 3.5e-4
  van_genuchten_m = 0.5
[]
```

`PorousFlowUnsaturated` is the companion Action to `PorousFlowBasicTHM` for unsaturated problems. It automatically includes:
- The capillary pressure (S-p relationship) materials
- The relative permeability model specified by `relative_permeability_type`
- The van Genuchten retention curve parameterised by `van_genuchten_alpha` and `van_genuchten_m`
- All kernels for the Richards mass balance

Setting `gravity = '0 -9.81 0'` is critical here: gravity drives downward infiltration. Without gravity, capillary suction alone would draw water into the dry soil, but the front would advance much more slowly and the capillary-rise physics would dominate over the gravitational-flow physics.

### `[Variables]`

```
[porepressure]
  initial_condition = -5e4   # -50 kPa
[]
```

The single primary variable is pore pressure. Negative pressure in unsaturated flow represents capillary suction — water is held by surface tension in the pore throats against atmospheric pressure. The initial -50 kPa corresponds to moderately dry sandy soil (about 6% effective saturation).

### `[Materials]`

The material block follows the same pattern as Cases 59 and 60. Note that `PorousFlowMatrixInternalEnergy` and `PorousFlowThermalConductivityIdeal` are included in this input file even though `coupling_type = Hydro`. This is a quirk of `PorousFlowUnsaturated`: certain material hierarchies require these even in isothermal mode. They do not affect the flow solution.

**`PorousFlowPorosity`** (`porosity_zero = 0.35`): representative of sandy soil.

**`PorousFlowPermeabilityConst`**: k = 1e-11 m^2. Saturated hydraulic conductivity K = k * rho * g / mu = 9.81e-5 m/s (approximately 0.01 cm/s, typical of a fine sand).

### `[BCs]`

```
[top_saturated]
  type     = DirichletBC
  variable = porepressure
  boundary = top
  value    = 0      # p = 0 Pa (atmospheric = saturated)
[]
```

Only the top boundary needs a BC. p = 0 at the surface represents ponded water with zero ponding depth — the pressure at the very top of the soil is at atmospheric (zero gauge). The bottom is left as natural zero-flux (free drainage).

### `[Executioner]`

```
type     = Transient
dt       = 60       # 60 s time steps
end_time = 3600     # 1 hour total
nl_max_its = 20
```

Sixty time steps of 60 s each. The `nl_max_its = 20` limit is set because Richards' equation has a very nonlinear pressure-saturation relationship, particularly in the dry regime where small pressure changes cause large saturation changes. Newton convergence is sometimes slow; the limit prevents runaway iterations.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case61-unsaturated-flow \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case61_unsaturated_flow.i 2>&1 | tail -30'
```

Output files:
- `case61_unsaturated_flow_out.e` — Exodus with pore pressure field at each time step
- `case61_unsaturated_flow_out.csv` — Time series of pressure at four heights

---

## Expected Results

The wetting front advances downward from the surface. Points above the front transition from dry (-50 kPa) to near-zero (wet) pressure; points below the front remain dry.

### Pressure at Monitoring Points over Time

| Time (s) | p_1p5m (Pa) | p_1p0m (Pa) | p_0p5m (Pa) | p_base (Pa) |
|----------|------------|------------|------------|------------|
| 0        | -50,000    | -50,000    | -50,000    | -50,000    |
| 600      | -6,165     | -40,619    | -49,728    | -47,603    |
| 1200     | -2,294     | -13,972    | -43,820    | -46,319    |
| 1800     | -1,348     | -5,679     | -25,979    | -41,700    |
| 3600     | -544       | -1,456     | -3,816     | -6,595     |

The pressure at 1.5 m height (p_1p5m) wets rapidly — within 600 s it has risen from -50 kPa to near zero because the wetting front passes through the upper part of the column early. At 1 hour, even the base (y = 0) is wetting with pressure at -6.6 kPa, indicating the front has nearly reached the bottom of the 2 m column.

### Wetting Front Progression

In the Exodus file, the pore pressure field shows a sharp wetting front advancing downward. The front is characteristically steep because:
1. The dry soil ahead of the front has very low relative permeability (k_r ~ S^2 << 1)
2. The wet soil behind the front is nearly saturated (k_r ~ 1)
3. Gravity and capillary forces both drive the front downward

The domain-average pressure `p_avg` rises from -43 kPa at t = 60 s toward zero as more of the column wets.

---

## Key Takeaways

- Richards' equation is the correct governing equation for variably saturated flow; it reduces to Darcy's law (fully saturated) when S = 1 everywhere, and to pure capillary diffusion when gravity is negligible.
- `PorousFlowUnsaturated` is the MOOSE Action for Richards' equation; it automatically incorporates the van Genuchten retention curve and the specified relative permeability model.
- The van Genuchten parameters alpha and m control the shape of the retention curve: larger alpha means water is released at lower suction (more sand-like); smaller alpha means tighter pores that hold water at higher suction (more clay-like).
- The Corey relative permeability `k_r = S_eff^n` captures the strong reduction in flow capacity in partially saturated soil; for n = 2 and S_eff = 0.1, the soil is 100 times less permeable than when fully saturated.
- The wetting front sharpness is a direct consequence of the nonlinear k_r(S) relationship: this is a physically real feature of unsaturated flow, not a numerical artifact.
- Gravity is essential for downward infiltration; without it, flow would only be driven by the capillary pressure gradient, which would smear the front and produce slower, more diffusive behaviour.
- Richards' equation is notoriously difficult to solve numerically in the very dry regime (near p = -50 kPa) because small pressure changes cause large saturation changes. The `nl_max_its` limit and NEWTON solve with MUMPS are essential for robust convergence.
