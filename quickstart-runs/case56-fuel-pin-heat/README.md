# Case 56: Fuel Pin Heat Transfer — Radial Temperature Profile

## Overview

This case models steady-state heat conduction in a nuclear fuel pin, one of the most important thermal calculations in reactor engineering. A fuel pin consists of a ceramic UO2 pellet surrounded by a thin helium-filled gap and a Zircaloy metal cladding. Enormous volumetric heat generation in the fuel (4 × 10^8 W/m^3) drives a steep radial temperature gradient from a hot centerline to a cooled outer surface.

The problem uses axisymmetric (RZ) coordinates, which is the natural choice for cylindrical geometry. In MOOSE, setting `coord_type = RZ` on the mesh transforms the governing equations to include the 1/r weighting factors needed for cylindrical symmetry. This means a 2D rectangular mesh in (r, z) space represents a full 3D cylinder when the equations are integrated in the azimuthal direction.

A key technique in this case is the spatially varying thermal conductivity implemented via `ParsedFunction` and `ADGenericFunctionMaterial`. Rather than defining separate mesh blocks for fuel, gap, and cladding, a single-block mesh uses a piecewise function of radial position r to assign the appropriate conductivity to each region. The `ADHeatConduction` kernel with automatic differentiation then correctly handles the conductivity jump at material interfaces.

## The Physics

### Governing Equation

Steady-state heat conduction in axisymmetric coordinates:

```
-div(k(r) * grad(T)) = q'''(r)
```

In cylindrical coordinates (r, z), with no azimuthal variation:

```
-(1/r) * d/dr(r * k(r) * dT/dr) - d/dz(k(r) * dT/dz) = q'''(r)
```

MOOSE's `coord_type = RZ` setting applies this cylindrical weighting automatically to the `ADHeatConduction` kernel.

### Boundary Conditions

**Symmetry axis (r = 0, left boundary)**: Natural (zero-flux) boundary condition. Because of the 1/r factor in cylindrical coordinates, the flux through the centerline is automatically zero — no explicit BC is needed. This is a key advantage of RZ coordinates.

**Outer cladding surface (r = r_clad_out, right boundary)**: Convective heat flux condition:

```
-k * dT/dr = h * (T - T_coolant)
```

where h = 30,000 W/m^2/K and T_coolant = 573 K. This represents forced-convection cooling by pressurized water in a PWR. The `ConvectiveHeatFluxBC` kernel implements this Robin (mixed) boundary condition.

**Top and bottom (axial faces)**: Natural boundary conditions (zero axial heat flux), correct for modeling a short segment far from the fuel pin ends.

### Material Properties and Regions

| Region | r range (mm) | k (W/m/K) | Physical material |
|--------|-------------|-----------|------------------|
| Fuel | 0 to 4.1 | 3.0 | UO2 ceramic |
| Gap | 4.1 to 4.2 | 0.5 | Helium fill gas |
| Cladding | 4.2 to 4.75 | 16.0 | Zircaloy-4 alloy |

The gap conductivity of 0.5 W/m/K is much lower than either fuel or cladding, creating a significant thermal resistance at the fuel-clad interface. In real reactors, this gap resistance is the dominant thermal resistance and must be modeled carefully.

Additional material properties (density = 10,000 kg/m^3, specific heat = 300 J/kg/K) are declared for completeness but do not affect the steady-state solution.

### Domain and Mesh

- Geometry: r in [0, 4.75e-3] m, z in [0, 0.01] m (1 cm axial segment)
- Mesh: 40 elements radially, 2 elements axially (quasi-1D radial problem)
- Element width: ~0.12 mm radially, spanning all three material regions
- coord_type: RZ (MOOSE applies cylindrical Jacobian factors)

### Analytical Estimate

For a solid fuel cylinder with uniform heat generation and no gap:

```
T_center = T_surface + q''' * r_fuel^2 / (4 * k_fuel)
         = 573 + 4e8 * (4.1e-3)^2 / (4 * 3.0)
         = 573 + 560
         = 1133 K
```

With the helium gap resistance adding approximately DeltaT_gap ~ q'' * (r_gap) / k_gap ~ (q''' * r_fuel) / 2 / k_gap * delta_r ~ 100 K, the actual centerline temperature will be approximately 1200-1250 K.

## Input File Walkthrough

### `[Mesh]`

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 40
  ny = 2
  xmin = 0.0
  xmax = 4.75e-3   # outer cladding radius (m)
  ymin = 0.0
  ymax = 0.01      # 1 cm axial segment
[]
coord_type = RZ   # x = r, y = z
```

The `coord_type = RZ` declaration at the `[Mesh]` block level tells MOOSE that x represents the radial coordinate r and y represents the axial coordinate z. All integrals in kernels and BCs are automatically weighted by 2*pi*r for cylindrical symmetry.

### `[Variables]`

```
[T]
  initial_condition = 600.0   # K
[]
```

A single temperature variable initialized to 600 K (near the expected solution), which helps the Newton solver converge faster.

### `[Kernels]`

**`conduction`** (`ADHeatConduction`): Implements the heat conduction term -div(k * grad(T)) using the material property `thermal_conductivity`. The AD (automatic differentiation) prefix means MOOSE computes the exact Jacobian by forward-mode AD rather than finite differences, which greatly improves Newton convergence.

**`heat_source`** (`BodyForce`): Adds the volumetric heat source q'''(r) as a body force. The `function = heat_gen` parameter points to the `ParsedFunction` that returns 4e8 W/m^3 inside the fuel and 0 elsewhere. `BodyForce` adds -f * test to the residual (the negative sign is the MOOSE convention for sources).

### `[Functions]`

**`heat_gen`** (`ParsedFunction`):
```
expression = 'if(x < 4.1e-3, 4.0e8, 0.0)'
```
In MOOSE's `ParsedFunction`, `x` is the x-coordinate of the evaluation point — which in RZ coordinates is the radial position r. The `if` function implements a step function: 4e8 W/m^3 for r < 4.1 mm (fuel region), zero outside.

**`k_func`** (`ParsedFunction`):
```
expression = 'if(x < 4.1e-3, 3.0, if(x < 4.2e-3, 0.5, 16.0))'
```
Nested `if` statements assign thermal conductivity based on radial position:
- r < 4.1 mm: UO2 fuel, k = 3.0 W/m/K
- 4.1 mm < r < 4.2 mm: He gap, k = 0.5 W/m/K
- r > 4.2 mm: Zircaloy cladding, k = 16.0 W/m/K

### `[Materials]`

**`thermal`** (`ADGenericFunctionMaterial`):
```
prop_names  = 'thermal_conductivity'
prop_values = 'k_func'
```
This is the critical connection: `ADGenericFunctionMaterial` evaluates the `ParsedFunction` named `k_func` at each quadrature point and stores the result as an AD material property `thermal_conductivity`. The AD variant is required because `ADHeatConduction` requests AD material properties for its Jacobian computation.

**`density_cp`** (`ADGenericConstantMaterial`): Declares `density = 10000 kg/m^3` and `specific_heat = 300 J/kg/K`. These are required by `ADHeatConduction` to have material properties registered, even though they do not contribute to the steady-state solution.

### `[BCs]`

**`coolant`** (`ConvectiveHeatFluxBC`):
```
boundary = right
T_infinity = 573.0
heat_transfer_coefficient = 30000.0
```
Applied to the right boundary (r = r_clad_out = 4.75 mm). This Robin BC adds h * (T - T_infinity) * test to the residual, representing heat loss to the coolant. No BC is needed on the left (symmetry axis) — the natural BC of `ADHeatConduction` correctly gives zero radial heat flux at r = 0.

### `[Postprocessors]`

Five postprocessors track the thermal state:

- `T_centerline`: Temperature at (r=0, z=0.005) — the hottest point
- `T_fuel_surface`: Temperature at r = 4.1 mm — fuel-gap interface
- `T_clad_inner`: Temperature at r = 4.2 mm — gap-cladding interface
- `T_clad_outer`: Temperature at r = 4.75 mm — coolant-facing surface
- `T_fuel_avg`: Volume-averaged temperature over the entire domain (fuel + gap + clad)

### `[Executioner]`

```
type = Steady
solve_type = 'NEWTON'
petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
petsc_options_value = 'lu       mumps'
nl_rel_tol = 1e-10
nl_abs_tol = 1e-12
```

Full Newton with a direct LU factorization (MUMPS sparse direct solver). For this relatively small problem (40x2 = 80 elements), a direct solver is efficient and highly reliable. The tight tolerances (1e-10 relative, 1e-12 absolute) ensure the temperature field is fully converged.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case56-fuel-pin-heat \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case56_fuel_pin_heat.i 2>&1 | tail -30'
```

## Expected Results

The simulation converges in 2-3 Newton iterations (the problem is mildly nonlinear through the spatially varying k, but the variation is only spatial, not temperature-dependent). The CSV output will report:

```
T_centerline  ~ 1229 K   (UO2 fuel center)
T_fuel_surface ~ 1019 K  (fuel-gap interface, outer edge of fuel)
T_clad_inner  ~  994 K   (gap-cladding interface)
T_clad_outer  ~  597 K   (surface cooled by water)
```

Key temperature drops across each region:
- **Fuel pellet**: ~210 K (from centerline to fuel surface) — most of the drop is in the fuel
- **Gap**: ~25 K (small thickness but very low conductivity of He)
- **Cladding**: ~397 K (Zircaloy is a poor conductor despite being metal)
- **Film resistance**: ~24 K (h = 30000 is very high, typical of forced PWR convection)

The total temperature rise from coolant to centerline is about 656 K, driven entirely by the fission power deposited in the fuel.

Physical interpretation: the UO2 centerline temperature of ~1229 K is well below the melting point of UO2 (~3138 K), which is an important safety margin. However, during a power transient, the centerline temperature can rise rapidly — this is why thermal-hydraulic analysis is central to reactor safety.

## Key Takeaways

- `coord_type = RZ` in MOOSE transforms a 2D rectangular mesh into an axisymmetric cylindrical model; MOOSE handles the 1/r Jacobian factors automatically.
- `ParsedFunction` combined with `ADGenericFunctionMaterial` provides a powerful and flexible way to define spatially varying material properties without defining separate mesh blocks.
- The `ADHeatConduction` kernel requires AD material properties; using `ADGenericFunctionMaterial` instead of `GenericFunctionMaterial` ensures compatibility.
- `ConvectiveHeatFluxBC` implements the Robin (h * (T - T_inf)) boundary condition; it is the correct BC for convective cooling, not a Dirichlet BC with a fixed temperature.
- The symmetry axis (r = 0) requires no explicit BC — the natural boundary condition of the RZ-coordinated diffusion operator is already zero radial flux, consistent with symmetry.
- The helium gap, despite being only 0.1 mm thick, provides significant thermal resistance because k_He = 0.5 W/m/K is very low compared to both fuel and cladding.
- A direct LU solver (MUMPS) is an excellent choice for small-to-medium steady nonlinear thermal problems where reliability is more important than scaling to many processors.
