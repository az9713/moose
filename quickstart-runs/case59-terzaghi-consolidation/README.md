# Case 59: Terzaghi 1D Consolidation

## Overview

This case solves the most fundamental problem in soil mechanics: a saturated soil column subjected to a sudden surface load. When a load is applied instantaneously, the pore fluid carries the entire load as excess pore pressure. Over time, as the fluid drains through the permeable top boundary, the excess pressure dissipates and the effective stress in the solid skeleton increases. This process — Terzaghi consolidation — controls settlement rates of buildings, embankments, and offshore foundations.

The governing equation is a linear diffusion equation in pore pressure, with the consolidation coefficient c_v playing the role of diffusivity. MOOSE's `PorousFlowBasicTHM` Action with `coupling_type = Hydro` assembles the single-phase, isothermal mass balance equation automatically. The result matches the classical analytical Terzaghi series solution.

Key concepts demonstrated:

- `PorousFlowBasicTHM` with `coupling_type = Hydro` for isothermal single-phase flow
- `PorousFlowConstantBiotModulus` for fluid-storage compressibility
- Drained top boundary (`p = 0`) with natural zero-flux bottom (undrained)
- Transient pressure dissipation matching the Terzaghi analytical series

---

## The Physics

### Governing Equation

In 1D (vertical direction z), the consolidation equation is:

```
d(p)/dt = c_v * d^2(p)/dz^2
```

where:
- `p` — excess pore pressure [Pa]
- `c_v = k / (mu * S)` — consolidation coefficient [m^2/s]
- `k = 1e-15 m^2` — intrinsic permeability
- `mu = 0.001 Pa.s` — fluid dynamic viscosity
- `S ~ 2e-10 /Pa` — specific storage (from Biot modulus)

This is derived from Darcy's law (the Darcy flux is proportional to the pressure gradient) combined with the fluid mass balance (what flows in must be stored or flow out).

### Consolidation Coefficient and Characteristic Time

```
c_v = k / (mu * S) = 1e-15 / (0.001 * 2e-10) = 5e-3 m^2/s
```

The characteristic time for the column to consolidate is:

```
t_c = H^2 / c_v = (10)^2 / 5e-3 = 20000 s
```

Note: the actual input uses k = 1e-15 m^2 which is very low permeability (tight clay). The simulation runs for 5000 s = 0.25 * t_c, so only partial dissipation is expected.

### Analytical Solution (Terzaghi, 1925)

For a column of height H with initial uniform excess pressure p0, drained at the top (z = H, p = 0) and undrained at the bottom (z = 0, zero flux), the solution is:

```
p(z, t) = sum_{n=0}^{inf} (4 * p0 / ((2n+1) * pi)) *
           sin((2n+1) * pi * z / (2H)) *
           exp(-(2n+1)^2 * pi^2 * c_v * t / (4 * H^2))
```

The first term (n=0) dominates at long times: a quarter-sine spatial profile decaying exponentially with time constant t_c / (pi^2 / 4) = 4 * t_c / pi^2.

### Degree of Consolidation

The degree of consolidation U quantifies how much excess pressure has dissipated:

```
U(t) = 1 - p_avg(t) / p0
```

U = 0 means undrained (no dissipation); U = 1 means fully consolidated. At t = 5000 s, the simulation reaches roughly U ~ 56% (p_avg ~ 44 kPa from initial 100 kPa).

### Domain and Boundary Conditions

```
y = 10 m  ----  TOP: drained (p = 0, free drainage)
          |    |
          | Saturated soil column  |
          | H = 10 m              |
          | k = 1e-15 m^2         |
          | phi = 0.4             |
          |    |
y = 0 m   ----  BOTTOM: undrained (zero flux, natural BC)

Mesh: 2 x 50 elements (quasi-1D in y direction)
Initial condition: p = 100 kPa everywhere
```

The left and right boundaries also use the natural zero-flux condition, consistent with the 1D nature of the problem.

---

## Input File Walkthrough

The input file is `case59_terzaghi_consolidation.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim  = 2
nx   = 2
ny   = 50
xmin = 0
xmax = 0.5
ymin = 0
ymax = 10
```

A quasi-1D mesh: the physics varies only in y (the vertical direction). Two elements in x are the minimum to avoid degenerate geometry. Fifty elements in y resolve the smooth pressure profile. The column height is 10 m.

### `[GlobalParams]`

```
PorousFlowDictator = dictator
```

All PorousFlow objects communicate through the `PorousFlowDictator` UserObject, which is created automatically by the `PorousFlowBasicTHM` Action. Declaring it in `[GlobalParams]` means every kernel and material finds it without repeating the name in each block.

### `[Variables]`

```
[porepressure]
  initial_condition = 1e5   # 100 kPa
[]
```

A single scalar field variable. The initial condition of 100 kPa represents the sudden load applied at t = 0 — before any drainage, the entire applied stress is carried as excess pore pressure.

### `[PorousFlowBasicTHM]` Action

```
[PorousFlowBasicTHM]
  porepressure    = porepressure
  coupling_type   = Hydro
  gravity         = '0 0 0'
  fp              = simple_fluid
  multiply_by_density = true
[]
```

This Action assembles the complete single-phase flow system:
- The `PorousFlowDictator` UserObject
- The fluid mass balance kernel with Darcy flux
- The time derivative of fluid mass storage
- All internal material linking objects

`coupling_type = Hydro` activates only the pressure equation (no temperature). `gravity = '0 0 0'` removes the hydrostatic component, isolating the excess pore pressure problem. `multiply_by_density = true` ensures the mass balance is in mass flux (rho * q) rather than volume flux, giving conservative mass accounting.

### `[Materials]`

**`PorousFlowPorosity`** (`porosity_zero = 0.4`): constant porosity with `mechanical = false`, `thermal = false`, `fluid = false` since there is no coupling to displacement, temperature, or fluid pressure in the porosity model.

**`PorousFlowConstantBiotModulus`**: computes the storage coefficient S from the Biot coefficient (1.0), solid bulk compliance (1e-10 /Pa), and fluid bulk modulus (2e9 Pa). The resulting storage coefficient is approximately 2e-10 /Pa. This controls how much pressure change occurs per unit mass of fluid added or removed — a high storage coefficient means the pressure diffuses slowly.

**`PorousFlowConstantThermalExpansionCoefficient`**: required by the Action even for Hydro-only coupling; all expansion coefficients are set to zero since there is no temperature effect.

**`PorousFlowPermeabilityConst`**: isotropic permeability tensor k = 1e-15 m^2 (9-component row-major tensor `'1e-15 0 0  0 1e-15 0  0 0 1e-15'`). This is characteristic of tight clay — orders of magnitude lower than sand.

### `[BCs]`

```
[top_drained]
  type     = DirichletBC
  variable = porepressure
  boundary = top
  value    = 0
[]
```

Only one BC is needed: the drained top surface where excess pore pressure is held at zero (atmospheric). The bottom, left, and right boundaries use the natural zero-flux Neumann condition automatically — no block is required.

### `[Postprocessors]`

| Name | Point | Physical meaning |
|------|-------|-----------------|
| `p_base` | (0.25, 0, 0) | Pressure at column base — slowest to drain |
| `p_mid` | (0.25, 5, 0) | Pressure at mid-height |
| `p_avg` | domain average | Overall consolidation state (U = 1 - p_avg/p0) |

### `[Executioner]`

```
type       = Transient
solve_type = NEWTON
dt         = 100       # 100 s steps
end_time   = 5000      # 5000 s total
```

Fifty time steps of 100 s each. `solve_type = NEWTON` uses the full analytical Jacobian provided by PorousFlow via automatic differentiation, giving fast quadratic convergence per step.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case59-terzaghi-consolidation \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case59_terzaghi_consolidation.i 2>&1 | tail -30'
```

Output files:
- `case59_terzaghi_consolidation_out.e` — Exodus with pore pressure field at each time step
- `case59_terzaghi_consolidation_out.csv` — Time series of p_base, p_mid, p_avg

---

## Expected Results

The simulation runs 50 time steps. The pore pressure at the base (`p_base`) starts at 100 kPa and decays as the pressure wave propagates upward toward the drained boundary.

### Pressure Evolution

| Time (s) | p_base (Pa) | p_mid (Pa) | p_avg (Pa) | U (degree of consol.) |
|----------|------------|-----------|-----------|----------------------|
| 0        | 100,000    | 100,000   | 100,000   | 0.00                 |
| 1000     | 99,411     | 88,797    | 75,078    | 0.25                 |
| 2000     | 94,593     | 74,024    | 64,538    | 0.35                 |
| 3000     | 86,349     | 63,593    | 56,491    | 0.44                 |
| 5000     | 68,753     | 48,925    | 43,952    | 0.56                 |

At t = 5000 s (0.25 * t_c), the simulation is partway through consolidation. The base pressure at ~68.8 kPa confirms the slow drainage through the low-permeability clay.

### Spatial Profile

The pressure profile in ParaView shows the characteristic shape of the Terzaghi solution: pressure remains near p0 at the bottom (far from the drained boundary), decreasing more steeply near the top. At early times the profile is nearly flat except near y = 10 m; at later times the decay penetrates deeper into the column.

The degree of consolidation U ~ 0.56 at t = 5000 s is consistent with the analytical series solution evaluated at t/t_c = 0.25.

---

## Key Takeaways

- Terzaghi consolidation is a pressure diffusion problem: the governing equation d(p)/dt = c_v * laplacian(p) is mathematically identical to the heat equation, with c_v playing the role of thermal diffusivity.
- `PorousFlowBasicTHM` with `coupling_type = Hydro` assembles the full single-phase flow system from a compact block; no individual kernels or fluid materials need to be listed.
- `PorousFlowConstantBiotModulus` controls the storage coefficient S, which determines c_v and therefore the consolidation rate; larger storage (more compressible system) slows consolidation.
- The natural (zero-flux) Neumann boundary condition is the correct undrained BC — it appears automatically when no DirichletBC is applied to a boundary.
- Low permeability (k = 1e-15 m^2) means extremely slow drainage; the consolidation time scale t_c = H^2 / c_v can span years for real clay layers, explaining why building settlements take decades.
- The degree of consolidation U = 1 - p_avg/p0 is the primary engineering quantity: U = 0.5 means half the settlement has occurred; U = 1 means full consolidation.
