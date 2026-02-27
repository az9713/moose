# Case 58: Control Rod Worth — Eigenvalue Shift from Absorber

## Overview

This case demonstrates control rod worth: the reduction in reactivity (and therefore k_eff) caused by inserting a neutron-absorbing material into the reactor core. Control rods are the primary mechanism for adjusting reactor power and for shutting down the reactor safely. Their effectiveness depends on the magnitude of the absorption cross section, the rod's position relative to the flux peak, and the spatial flux depression they cause.

The model uses a 20 cm slab with a central 4 cm absorber region (x in [8, 12] cm) where the macroscopic absorption cross section is enhanced from 0.08 /cm to 0.13 /cm. The eigenvalue solver finds k_eff for this heterogeneous configuration. The result can be compared to the unrodded case (k_eff ~ 1.146) to compute the rod worth in units of delta-k or percent mille (pcm).

A key numerical technique introduced here is using `ADGenericFunctionMaterial` and `ADMatReaction` to implement a spatially varying absorption cross section without requiring separate mesh blocks. The `ParsedFunction` approach is more flexible than block-based materials when the region boundaries do not align with element boundaries, and it naturally produces a smooth material property field evaluated at quadrature points.

## The Physics

### Governing Equation

One-group neutron diffusion with spatially varying absorption:

```
-D * laplacian(phi) + Sigma_a(x) * phi = (1/k) * nu_Sigma_f * phi
```

where the absorption cross section is piecewise:

```
Sigma_a(x) = 0.13 /cm    for x in [8, 12] cm   (rod region)
Sigma_a(x) = 0.08 /cm    for x elsewhere        (base fuel)
```

The fission cross section nu_Sigma_f = 0.12 /cm is uniform throughout. Inserting the rod does not affect fission (the rod is purely absorbing), so the fission source is uniform while the loss term is spatially enhanced.

### Analytical Estimates

**Without control rod (uniform Sigma_a = 0.08):**

The geometric buckling for a 20 cm slab with vacuum BCs:
```
B^2 = (pi/L)^2 = (pi/20)^2 = 0.02467 /cm^2
```

The one-group k_eff formula:
```
k_eff = nu_Sigma_f / (Sigma_a + D * B^2)
      = 0.12 / (0.08 + 1.0 * 0.02467)
      = 0.12 / 0.1047
      = 1.146
```

**With control rod (Sigma_a enhanced in center):**

The rod suppresses flux in the high-importance central region. Because the central flux is the highest (fundamental mode peaks at center), the rod is inserted in the most effective location. The simple one-group analytical formula does not apply directly to the heterogeneous case; the MOOSE eigenvalue calculation gives k_eff ~ 1.019.

**Control rod worth:**

```
rho = (k_rod - k_no_rod) / k_rod / k_no_rod    (in delta-k units)
    = (1.019 - 1.146) / (1.019 * 1.146)
    ~ -0.109   (about -11 cents or -1090 pcm)

Or more simply:
Delta_rho = (1 - 1/k_rod) - (1 - 1/k_no_rod) = 1/k_no_rod - 1/k_rod
           = 1/1.146 - 1/1.019
           ~ 0.873 - 0.981
           ~ -0.108   (~-10800 pcm)
```

The rod worth of ~11,000 pcm means this rod would more than shut down the reactor if it were critical (k_eff = 1.0 is 0 pcm; shutdown requires negative reactivity).

### Sign Convention in `ADMatReaction`

The `ADMatReaction` kernel computes the residual contribution as:

```
R = -reaction_rate * phi * test
```

Note the negative sign. For absorption (a loss term), we want +Sigma_a * phi * test in the residual. Therefore, we must set `reaction_rate = -Sigma_a` (negative of the physical cross section). This is why the `ParsedFunction` for `sigma_a_func` returns negative values (-0.13 and -0.08).

This sign convention differs from `Reaction` (which has R = +rate * phi * test). The choice of which kernel to use depends on whether the material property is available as an AD material property (required for `ADMatReaction`) or as a plain scalar (suitable for `Reaction`).

### Flux Depression in the Rod Region

With the rod inserted, the flux profile is no longer a simple cosine. The enhanced absorption in [8, 12] cm acts as a local neutron sink, depressing the flux in that region. The shape is:

- Outside the rod: approximately cosine, but with a steeper curvature due to leakage toward the rod
- Inside the rod: a flatter, suppressed flux profile with a local minimum at the rod center
- The ratio phi_center/phi_quarter (rod center vs. quarter point) is significantly less than 1.0 for the rodded case vs. approximately 0.7 for the unrodded cosine (phi at x=5 vs. x=10 for a cosine on [0,20])

### Domain and Mesh

- Geometry: 1D slab, x in [0, 20] cm, modeled as 2D quasi-1D (100 x 2 elements)
- Rod region: x in [8, 12] cm, 20 elements (element width 0.2 cm, fine enough to resolve the flux depression)
- Both regions use the same mesh; material properties vary continuously at quadrature points

## Input File Walkthrough

### `[Mesh]`

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 100
  ny = 2
  xmin = 0
  xmax = ${L}    # 20.0 cm
  ymin = 0
  ymax = 0.5
[]
```

A uniform 100-element mesh. The rod region [8, 12] covers 20% of the domain and contains 20 elements. No mesh refinement at the rod boundaries is used; the sharp cross-section step is captured by the quadrature-point evaluation of the `ParsedFunction`.

### `[Variables]`

```
[phi]
  initial_condition = 1.0
[]
```

Single flux variable, initialized to 1.0 uniformly. The eigenvalue solver normalizes the mode shape, so the initial condition only affects the starting point for iteration.

### `[Functions]`

```
[sigma_a_func]
  type = ParsedFunction
  expression = 'if(x > 8.0 & x < 12.0, -0.13, -0.08)'
[]
```

Returns the negated absorption cross section (because `ADMatReaction` uses R = -reaction_rate * phi * test). The `&` operator in `ParsedFunction` expressions is the logical AND. Values:
- Rod region: -0.13 = -(0.08 + 0.05), total absorption including rod penalty
- Outside rod: -0.08, base fuel absorption only

### `[Kernels]`

**`diffusion`** (`MatDiffusion`): Standard diffusion term -D * laplacian(phi) using the constant material property D = 1.0.

**`absorption`** (`ADMatReaction`): Spatially varying absorption term. Uses the AD material property `sigma_a_ad` evaluated from `sigma_a_func`. The kernel computes R = -sigma_a_ad * phi * test. Since sigma_a_ad is negative (from the ParsedFunction), this gives R = +|Sigma_a| * phi * test — a positive residual contribution, meaning a loss term.

**`fission`** (`Reaction`, rate = -0.12, `extra_vector_tags = 'eigen'`): Uniform fission source -(1/k) * nu_Sigma_f * phi. Negative rate = source. Tagged eigen for the B-matrix.

### `[Materials]`

**`diffusion_coeff`** (`GenericConstantMaterial`): Constant D = 1.0 for the `MatDiffusion` kernel.

**`sigma_a_mat`** (`ADGenericFunctionMaterial`):
```
prop_names = 'sigma_a_ad'
prop_values = 'sigma_a_func'
```
Converts the `ParsedFunction` into an AD material property. The AD suffix on `ADGenericFunctionMaterial` is required because `ADMatReaction` requests AD properties for Jacobian computation via automatic differentiation.

### `[BCs]`

Four boundary conditions at x = 0 and x = L (same pattern as Cases 54 and 55):
- `DirichletBC`: phi = 0 for the A-matrix (absorption + diffusion)
- `EigenDirichletBC`: phi = 0 for the B-matrix (fission eigen source)

Both must be applied to enforce zero-flux vacuum boundaries in the complete generalized eigenproblem.

### `[VectorPostprocessors]`

```
[flux_profile]
  type = LineValueSampler
  variable = phi
  start_point = '0 0.25 0'
  end_point = '${L} 0.25 0'
  num_points = 101
  sort_by = x
[]
```

101 sample points along the slab centerline. The rod-induced flux depression is clearly visible in this profile as a dip around x = 10 cm relative to what a pure cosine would show.

### `[Postprocessors]`

- `phi_center` (x = 10): Flux at center of rod — shows depression
- `phi_quarter` (x = 5): Flux at quarter point — outside rod, elevated relative to center

The ratio phi_center/phi_quarter < 1 (inverted from the unrodded cosine where center > quarter) demonstrates the rod-induced flux peaking at off-center locations.

### `[Executioner]`

```
type = Eigenvalue
solve_type = PJFNK
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre    boomeramg'
```

Same as Cases 54 and 55. The heterogeneous cross-section profile increases the condition number of the system slightly compared to the uniform case, but `boomeramg` handles this well.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case58-control-rod-worth \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case58_control_rod_worth.i 2>&1 | tail -30'
```

## Expected Results

The eigenvalue solver converges to:

```
k_eff ~ 1.019    (with control rod fully inserted in central 4 cm)
```

Compare to the analytical unrodded value:
```
k_eff_no_rod ~ 1.146    (uniform Sigma_a = 0.08 /cm, 20 cm slab)
```

The flux profile from `flux_profile_0001.csv` shows a characteristic "top-hat depression" in the rod region:

```
x = 0:    phi = 0.000   (vacuum BC)
x = 5:    phi ~ 0.75    (quarter point, outside rod, elevated)
x = 8:    phi ~ 0.65    (just entering rod, flux starts dropping)
x = 10:   phi ~ 0.55    (rod center, flux minimum)
x = 12:   phi ~ 0.65    (exiting rod, flux recovers)
x = 15:   phi ~ 0.75    (quarter point on other side, symmetric)
x = 20:   phi = 0.000   (vacuum BC)
```

Values are SLEPc-normalized and may differ by an overall scale factor. The key observable is the ratio phi_center/phi_quarter, which should be noticeably less than 1.0 (flux peaked off-center due to the rod).

The reactivity worth of the rod is approximately:
```
Delta_rho = 1/k_no_rod - 1/k_rod = 1/1.146 - 1/1.019 ~ -0.107
```
or about -10,700 pcm (percent mille = 10^-5 delta-k). This is very large — a rod worth this much could shut down even a substantially supercritical reactor.

## Key Takeaways

- `ADMatReaction` has residual R = -reaction_rate * phi * test (note the negative sign); for a loss term like absorption, the material property must be provided as the negative of the physical cross section.
- `ADGenericFunctionMaterial` converts a `ParsedFunction` to an AD material property; the AD version is required when the consuming kernel (`ADMatReaction`) uses automatic differentiation.
- `ParsedFunction` with logical operators (`&`, `|`) enables piecewise material properties on a single mesh block, avoiding the need to define and track separate element subdomains.
- Both `DirichletBC` and `EigenDirichletBC` are needed for eigenvalue problems with Dirichlet constraints — one enforces zero flux in the A-matrix (diffusion + absorption), the other in the B-matrix (fission source).
- Control rod worth depends strongly on placement: a rod at the flux peak is maximally effective (highest importance), while a rod at the flux minimum contributes little. This is the principle of "rod worth" in reactor physics.
- The flux depression in the rod region redistributes power to peripheral fuel, which must be accounted for in thermal analysis to avoid hot spots outside the rod region.
- Comparing two k_eff values (rodded vs. unrodded) requires running two separate eigenvalue calculations; this input file computes only the rodded case. The unrodded case would use uniform Sigma_a = 0.08 /cm.
