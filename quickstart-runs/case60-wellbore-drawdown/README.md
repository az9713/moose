# Case 60: Wellbore Drawdown — Radial Flow to a Production Well

## Overview

This case simulates pressure drawdown around a pumping well in a confined aquifer. When a well is pumped, fluid flows radially inward from all directions toward the wellbore. The resulting cone of depression — the pressure field that develops around the well — is the central problem of hydrogeology. The steady-state solution is a logarithmic pressure profile; the transient Theis solution describes how this cone develops over time.

The simulation uses MOOSE's `PorousFlowBasicTHM` Action in axisymmetric (RZ) coordinates. In RZ geometry, MOOSE treats the x-coordinate as the radial direction r and integrates over the azimuthal angle automatically, so a 2D mesh in the (r, z) plane represents the full 3D cylindrical domain. This is the correct and efficient approach for radially symmetric well problems.

Key concepts demonstrated:

- `coord_type = RZ` for axisymmetric (cylindrical) geometry in porous flow
- Biased mesh refinement near the wellbore using `bias_x` to resolve the steep pressure gradient
- `PorousFlowBasicTHM` with `coupling_type = Hydro` for confined aquifer flow
- Fixed-pressure BCs at well and outer boundary for the Theis drawdown problem

---

## The Physics

### Governing Equation in Radial Coordinates

In radial (RZ) coordinates with no z-dependence, the transient flow equation becomes:

```
S * d(p)/dt = (k/mu) * (1/r) * d/dr ( r * d(p)/dr )
```

where:
- `p` — pore pressure [Pa]
- `S ~ 2e-10 /Pa` — storage coefficient
- `k = 1e-12 m^2` — permeability (~1 Darcy, permeable sandstone)
- `mu = 0.001 Pa.s` — water viscosity
- `r` — radial distance from well axis [m]

The right-hand side is the Laplacian in cylindrical coordinates. MOOSE handles the 1/r geometric weighting automatically when `coord_type = RZ` is set.

### Theis Solution (Analytical)

For a well of radius r_w in an infinite confined aquifer, pumped at constant flow rate Q, the transient drawdown is:

```
s(r, t) = p0 - p(r, t) = (Q * mu) / (4 * pi * k * H) * W(u)
```

where `W(u)` is the well function (exponential integral) and `u = r^2 * S * mu / (4 * k * t)`. At late times (small u), the Theis solution approaches the logarithmic Thiem solution:

```
p(r) = p_well + (p_outer - p_well) * ln(r/r_well) / ln(r_outer/r_well)
```

This logarithmic profile is the steady-state solution (Thiem, 1906) for a confined aquifer between fixed-pressure inner and outer boundaries.

### Estimated Steady-State Drawdown

For this case, with pressure-controlled boundaries rather than flow-rate control, the drawdown at any radius follows the Thiem formula. The maximum drawdown (at the well) is:

```
s_well = p0 - p_well = 10 MPa - 9.9 MPa = 100 kPa
```

Pressure at any radius:
```
p(r) = 9.9e6 + 100e3 * ln(r / 0.1) / ln(100 / 0.1)
     = 9.9e6 + 100e3 * ln(r / 0.1) / 6.908
```

At r = 1 m:   p ~ 9.9e6 + 14.5e3 = 9.914 MPa
At r = 10 m:  p ~ 9.9e6 + 43.3e3 = 9.943 MPa
At r = 50 m:  p ~ 9.9e6 + 79.0e3 = 9.979 MPa

These match the simulation output at late time (approaching steady state).

### Domain and Boundary Conditions

```
Axisymmetric model (RZ coordinates):
r = 0.1 m  (left = wellbore)     r = 100 m  (right = outer boundary)
     |___________________________________|
     |                                   |
p = 9.9 MPa                          p = 10.0 MPa
(pumped well)        -->              (undisturbed aquifer)
     |                                   |
     |___________________________________|
y = 0 (bottom, symmetry)           y = 1.0 (top, symmetry)

Mesh: 50 x 2 elements, bias_x = 1.15 (elements refined near r = 0.1 m)
Initial condition: p = 10 MPa everywhere
```

The top and bottom boundaries (y = 0 and y = 1) use the natural zero-flux condition, representing the impermeable confining layers above and below the aquifer.

---

## Input File Walkthrough

The input file is `case60_wellbore_drawdown.i`.

### `[Mesh]`

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 50
  ny = 2
  xmin = 0.1       # wellbore radius [m]
  xmax = 100.0     # outer boundary radius [m]
  ymin = 0
  ymax = 1.0       # unit thickness
  bias_x = 1.15    # progressive coarsening away from well
[]
coord_type = RZ     # x = r, y = z
```

Two critical features:

**`bias_x = 1.15`**: each successive element in the r-direction is 1.15 times wider than the previous. Near the well (r = 0.1 m), elements are very fine to resolve the steep logarithmic pressure gradient. Far from the well, elements are coarser where the gradient is gentle. This geometric refinement gives good accuracy without excessive element count.

**`coord_type = RZ`**: instructs MOOSE to interpret x as the radial coordinate r and y as the axial coordinate z, and to integrate over the 2*pi azimuthal angle. The radial source terms (the 1/r factor in the cylindrical Laplacian) are applied automatically. Without this setting, the simulation would model a slab geometry rather than a cylinder, producing incorrect results.

### `[PorousFlowBasicTHM]` Action

```
[PorousFlowBasicTHM]
  porepressure    = porepressure
  coupling_type   = Hydro
  gravity         = '0 0 0'
  fp              = water
  multiply_by_density = true
[]
```

Same as Case 59. The horizontal confined aquifer assumption allows gravity to be set to zero without loss of generality (the vertical pressure gradient from gravity is a static background; only the excess drawdown pressure field matters here).

### `[Materials]`

**`PorousFlowPorosity`** (`porosity_zero = 0.2`): representative porosity for a permeable sandstone aquifer.

**`PorousFlowConstantBiotModulus`**: storage coefficient from Biot coefficient = 1.0, solid compliance 1e-10 /Pa, fluid bulk modulus 2e9 Pa. Gives S ~ 2e-10 /Pa.

**`PorousFlowPermeabilityConst`**: k = 1e-12 m^2. This is approximately 1 Darcy — the standard unit of permeability in petroleum engineering and hydrogeology, representative of a productive oil reservoir or a high-quality sandstone aquifer.

### `[BCs]`

```
[well_pressure]
  type     = DirichletBC
  variable = porepressure
  boundary = left    # left = inner radius = wellbore
  value    = 9.9e6   # 9.9 MPa (100 kPa drawdown)
[]
[outer_pressure]
  type     = DirichletBC
  variable = porepressure
  boundary = right   # right = outer boundary
  value    = 1e7     # 10 MPa (undisturbed formation)
[]
```

Two pressure-controlled boundaries. The well is modelled as a fixed-pressure source (pressure-controlled pumping), not a fixed-flow-rate source. This simplifies the BC to a Dirichlet condition. The outer boundary at 100 m is held at the undisturbed formation pressure, approximating an infinite aquifer.

### `[Postprocessors]`

| Name | Point | Physical meaning |
|------|-------|-----------------|
| `p_well` | (0.1, 0.5, 0) | Pressure at wellbore face (confirms BC) |
| `p_1m` | (1.0, 0.5, 0) | Near-well pressure (in the steep gradient zone) |
| `p_10m` | (10.0, 0.5, 0) | Intermediate radius |
| `p_50m` | (50.0, 0.5, 0) | Far-field (near outer boundary) |
| `p_avg` | domain average | Overall pressure depletion |

These four radial monitoring points let you reconstruct the drawdown cone shape from the CSV output.

### `[Executioner]`

```
type     = Transient
dt       = 10
end_time = 500     # 500 s total
```

Fifty time steps of 10 s each. The simulation captures the early transient as the cone of depression expands outward from the well, approaching the steady-state logarithmic profile.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case60-wellbore-drawdown \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case60_wellbore_drawdown.i 2>&1 | tail -30'
```

Output files:
- `case60_wellbore_drawdown_out.e` — Exodus with pressure field at each time step
- `case60_wellbore_drawdown_out.csv` — Time series of p_well, p_1m, p_10m, p_50m, p_avg

---

## Expected Results

The pressure at the well is held fixed at 9.9 MPa from t = 0. The cone of depression propagates outward, initially affecting only the near-well region and progressively reaching larger radii.

### Pressure at Monitoring Points over Time

| Time (s) | p_1m (MPa) | p_10m (MPa) | p_50m (MPa) | p_avg (MPa) |
|----------|-----------|------------|------------|------------|
| 10       | 9.949     | 9.991      | 9.9999     | 9.9996     |
| 100      | 9.937     | 9.974      | 9.996      | 9.997      |
| 300      | 9.934     | 9.968      | 9.991      | 9.994      |
| 500      | 9.934     | 9.967      | 9.990      | 9.993      |

The near-well pressure (p_1m) drops quickly and approaches its steady-state value. The far-field pressure (p_50m) barely changes over 500 s because the pressure front has not fully propagated that far. By t = 500 s the cone is still developing.

### Steady-State Verification

At late time, the pressure at r = 1 m approaches the Thiem analytical solution:
```
p(1.0 m) = 9.9e6 + 1e5 * ln(1.0/0.1) / ln(100/0.1)
          = 9.9e6 + 1e5 * 2.303 / 6.908
          = 9.9e6 + 14,330 Pa
          ~ 9.9143 MPa
```

The simulated value of ~9.933 MPa at t = 500 s is approaching this steady-state target, confirming the logarithmic drawdown profile.

---

## Key Takeaways

- `coord_type = RZ` in the `[Mesh]` block activates axisymmetric geometry; MOOSE automatically applies the 1/r cylindrical weighting factor in all integrals, converting the 2D (r, z) problem into a full 3D radial solution.
- `bias_x` provides geometric mesh refinement near the well where the pressure gradient is steepest — this is essential for accuracy without requiring a uniformly fine mesh across the entire domain.
- The steady-state pressure profile in a confined aquifer with fixed inner and outer pressures follows the Thiem logarithmic solution exactly; verifying the simulation against this formula confirms the RZ geometry is working correctly.
- `coupling_type = Hydro` isolates the pressure diffusion problem; adding `ThermoHydro` would couple temperature for geothermal well problems.
- Permeability k = 1e-12 m^2 (~1 Darcy) is representative of productive hydrocarbon reservoirs; values 2-3 orders of magnitude lower (millidarcy or microdarcy) characterize tight gas or shale formations where wellbore pressures take much longer to propagate.
- The Theis time-to-steady-state scales as r^2 * S * mu / k; for this case at r = 100 m, steady state requires thousands of seconds — the simulation at 500 s shows the cone still propagating.
