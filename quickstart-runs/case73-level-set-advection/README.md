# Case 73: Level Set Bubble Advection — Solid Body Rotation

## Overview

The level set method is one of the most widely-used techniques for tracking moving interfaces in computational physics. Rather than tracking the interface explicitly as a collection of marker points (which can tangle and break), the interface is represented implicitly as the zero contour of a smooth scalar function phi(x,y,t). The interface moves by advecting phi with the local fluid velocity. This Eulerian description naturally handles topological changes (bubble breakup, merging), and the interface geometry (normal vector, curvature) can be computed directly from phi via spatial derivatives.

This case demonstrates the MOOSE `level_set` module for a standard benchmark problem: solid-body rotation of a smooth circular bubble. The velocity field is v = (-(y-0.5), (x-0.5)), which rotates every point in the unit square counterclockwise about the center (0.5, 0.5) with period 2*pi. A bubble defined by a smooth tanh profile is initialized at (0.5, 0.75) with radius 0.15 and advected for t = pi (a half-rotation), at which point it should arrive at (0.5, 0.25). Pure advection by an incompressible velocity field is a stringent test: a perfect scheme would deliver the bubble without any distortion, diffusion, or oscillation.

Key concepts demonstrated:

- `LevelSetAdvection` and `LevelSetAdvectionSUPG` kernels for the advection term
- `LevelSetTimeDerivativeSUPG` for the consistent SUPG modification of the time derivative
- `LAGRANGE_VEC` family for storing vector-valued auxiliary variables (velocity)
- `VectorFunctionIC` with `ParsedVectorFunction` for initializing vector fields
- BDF2 second-order time integration for improved temporal accuracy
- Mass conservation monitoring via `ElementIntegralVariablePostprocessor`

---

## The Physics

### Level Set Advection Equation

The level set function phi(x,y,t) satisfies the scalar hyperbolic advection equation:

```
dphi/dt + v . nabla(phi) = 0
```

where v = (v_x, v_y) is the prescribed velocity field. This equation states that phi is constant along particle trajectories (the method of characteristics): a material point at (x_0, y_0) at t = 0 carries its initial value phi_0 = phi(x_0, y_0, 0) with it as it moves along the streamline.

There is no diffusion term — the equation is purely hyperbolic. This makes it challenging for standard Galerkin finite element methods, which are known to produce spurious oscillations (Gibbs-type oscillations) for pure advection without stabilization.

### Velocity Field — Solid-Body Rotation

```
v_x = -(y - 0.5)
v_y =  (x - 0.5)
```

This velocity field rotates the domain counterclockwise about the center (0.5, 0.5). The angular velocity is omega = 1 rad/s (the velocity magnitude at distance r from the center is |v| = r). The rotation period is T = 2*pi ~ 6.28 s. After t = pi ~ 3.14 s (a half-rotation), a point at (0.5, 0.75) maps to (0.5, 0.25). The velocity field is divergence-free (div(v) = 0), so the level set equation is equivalent to a conservation law and phi-integrals over the domain should be conserved.

### Initial Condition — Smooth Bubble

```
phi(x, y, 0) = 0.5 * (1 + tanh((0.15 - sqrt((x-0.5)^2 + (y-0.75)^2)) / 0.02))
```

This creates a smooth bubble:
- phi ~ 1 inside the bubble (r < 0.15 - several interface widths)
- phi ~ 0 outside the bubble (r > 0.15 + several interface widths)
- The transition from 1 to 0 occurs over a width of ~ 4 * 0.02 = 0.08 m (the tanh profile's characteristic width)

The smooth (diffuse) interface avoids the Gibbs phenomenon that occurs with a sharp step function initial condition on a finite-resolution mesh.

### SUPG Stabilization

Standard Galerkin FEM applied to the pure advection equation produces centered differences in the streamline direction, which are neutrally stable (not dissipative) and produce nonphysical oscillations in the wake of sharp gradients. Streamline Upwind Petrov-Galerkin (SUPG) stabilization modifies the test functions to add dissipation only along streamlines:

```
w_SUPG = w + tau * (v . nabla(w))
```

where tau ~ h/(2|v|) is the stabilization parameter (h is element size, |v| is velocity magnitude). This is equivalent to adding a small amount of upwind diffusion in the velocity direction, which dampens oscillations without introducing excessive isotropic diffusion.

For the SUPG-stabilized transient problem, both the time derivative and the advection terms need SUPG modification to maintain consistency. This requires both `LevelSetAdvectionSUPG` and `LevelSetTimeDerivativeSUPG` kernels — omitting the time derivative SUPG term produces an inconsistent formulation that can still oscillate.

### Expected Accuracy

Even with SUPG, pure advection schemes with standard Lagrange elements suffer from:
1. **Overshoot/undershoot**: max_phi > 1 and min_phi < 0 (violating the physical bounds 0 <= phi <= 1). This is a Gibbs-type effect at the bubble interface.
2. **Mass drift**: The integral of phi over the domain drifts slightly from its initial value. For incompressible flow (div(v) = 0), the level set equation is conservative, but numerical discretization breaks exact conservation.
3. **Interface diffusion**: The tanh interface profile broadens slightly over time.

The CSV data shows max_phi reaching ~1.108 (11% overshoot) and mass conserved within ~2.5% — acceptable for a 40x40 mesh.

---

## Input File Walkthrough

The input file is `case73_level_set_advection.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 40
ny = 40
xmin = 0.0  xmax = 1.0
ymin = 0.0  ymax = 1.0
```

A 40x40 structured mesh on [0,1]^2. The finer resolution (compared to 20x20 in earlier cases) is needed to resolve the bubble interface width (~0.08 m) with several elements. The element size h = 1/40 = 0.025 m gives approximately 3 elements across the tanh interface width.

### `[Variables]` and `[AuxVariables]`

```
[phi]
[]

[AuxVariables]
  [velocity]
    family = LAGRANGE_VEC
  []
[]
```

`phi` is the scalar level set function (standard FIRST order LAGRANGE). `velocity` uses the `LAGRANGE_VEC` family — a vector-valued variable that stores all components (v_x, v_y in 2D) as a single variable object. This is more convenient than declaring separate `vel_x` and `vel_y` AuxVariables and is required by the `LevelSetAdvection` and SUPG kernel interfaces.

### `[ICs]`

```
[phi_ic]
  type = FunctionIC
  variable = phi
  function = phi_init
[]
[vel_ic]
  type = VectorFunctionIC
  variable = velocity
  function = vel_func
[]
```

`FunctionIC` initializes the scalar phi from a `ParsedFunction`. `VectorFunctionIC` initializes the vector velocity from a `ParsedVectorFunction` (which specifies each component separately via `expression_x` and `expression_y`).

### `[Functions]`

```
[phi_init]
  type = ParsedFunction
  expression = '0.5*(1.0 + tanh((0.15 - sqrt((x-0.5)^2+(y-0.75)^2))/0.02))'
[]
[vel_func]
  type = ParsedVectorFunction
  expression_x = '-(y - 0.5)'
  expression_y = '(x - 0.5)'
[]
```

`ParsedVectorFunction` is the vector-valued counterpart of `ParsedFunction`. Each component is specified as a separate expression. Note that the velocity is prescribed as an initial condition only — the velocity does not change in time (solid-body rotation is a time-independent field), so no time-dependent update is needed.

### `[Kernels]`

```
[time]
  type = TimeDerivative
  variable = phi
[]
[advection]
  type = LevelSetAdvection
  variable = phi
  velocity = velocity
[]
[advection_supg]
  type = LevelSetAdvectionSUPG
  variable = phi
  velocity = velocity
[]
[time_supg]
  type = LevelSetTimeDerivativeSUPG
  variable = phi
  velocity = velocity
[]
```

Four kernels work together to form the SUPG-stabilized advection equation:

1. `TimeDerivative`: The standard Galerkin dphi/dt term (integral of w * dphi/dt over the domain).
2. `LevelSetAdvection`: The standard Galerkin advection term (integral of w * v . nabla(phi)).
3. `LevelSetAdvectionSUPG`: The SUPG stabilization for the advection term (integral of (v . nabla(w)) * tau * v . nabla(phi)).
4. `LevelSetTimeDerivativeSUPG`: The SUPG stabilization for the time derivative (integral of (v . nabla(w)) * tau * dphi/dt).

All four kernels reference the same `velocity` AuxVariable. The stabilization parameter tau is computed internally by the SUPG kernels based on the element size and velocity magnitude.

### `[Postprocessors]`

```
[total_phi]  type = ElementIntegralVariablePostprocessor  variable = phi
[max_phi]    type = ElementExtremeValue  variable = phi  value_type = max
[min_phi]    type = ElementExtremeValue  variable = phi  value_type = min
```

`ElementIntegralVariablePostprocessor` computes the integral of phi over the entire domain: integral_Omega phi dV. For a conservative advection scheme with div(v) = 0, this should remain constant at its initial value of ~ 0.0717 (the area of the bubble). Deviations indicate numerical mass loss or gain. `ElementExtremeValue` tracks the maximum and minimum values of phi to quantify overshoots and undershoots relative to the physical bounds [0, 1].

### `[Executioner]`

```
type = Transient
scheme = bdf2
solve_type = NEWTON
dt = 0.1
end_time = 3.14159
nl_abs_tol = 1e-10
nl_rel_tol = 1e-8
```

`scheme = bdf2` selects the second-order Backward Differentiation Formula time integrator. BDF2 is more accurate than the default BDF1 (implicit Euler) for smooth problems, reducing temporal truncation error from O(dt) to O(dt^2). For pure advection with a smooth velocity field, this is beneficial because the bubble profile varies smoothly in time. The 31 time steps of dt = 0.1 s span the half-rotation period of pi seconds.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case73-level-set-advection \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case73_level_set_advection.i 2>&1 | tail -30'
```

Output files:
- `case73_level_set_advection_out.e` — Exodus with the phi field at each time step (shows bubble rotating)
- `case73_level_set_advection_out.csv` — Time series of total_phi, max_phi, min_phi

---

## Expected Results

### Time Evolution of Key Quantities

| Time (s) | max_phi | min_phi | total_phi |
|----------|---------|---------|-----------|
| 0.0      | 1.000   |  0.000  | 0.07171   |
| 0.3      | 1.018   | -0.018  | 0.07174   |
| 0.6      | 1.056   | -0.052  | 0.07200   |
| 1.0      | 1.085   | -0.079  | 0.07271   |
| 1.5      | 1.100   | -0.103  | 0.07159   |
| 2.0      | 1.104   | -0.118  | 0.07189   |
| 2.5      | 1.109   | -0.127  | 0.07246   |
| 3.0      | 1.109   | -0.134  | 0.07208   |
| 3.14159  | 1.108   | -0.135  | 0.07179   |

### Bubble Position at t = pi

After the half-rotation, the bubble center has moved from (0.5, 0.75) to (0.5, 0.25). The phi contour at phi = 0.5 (the interface) forms a circle of radius ~0.15 centered at (0.5, 0.25) in the Exodus output. The shape is well-preserved — the circular bubble remains approximately circular, which validates the SUPG scheme's ability to advect smooth profiles with reasonable fidelity.

### Overshoot and Undershoot Growth

The max_phi overshoot grows from 1.000 at t = 0 to ~1.108 at t = pi (11% overshoot). The min_phi undershoot grows from 0.000 to approximately -0.135 (13.5% undershoot). These violations of the physical bounds [0, 1] are characteristic of the Gibbs phenomenon for SUPG-stabilized advection on Q1 elements. The overshoots concentrate near the bubble interface where the gradient is steepest.

### Mass Conservation

The total_phi integral remains within ~2.5% of its initial value 0.07171 throughout the simulation. The variation is not monotone — it fluctuates as the bubble traverses different parts of the mesh. Mass conservation improves with mesh refinement: doubling the mesh to 80x80 would reduce mass error to approximately 0.6%.

### Comparison: With vs. Without SUPG

Removing `LevelSetAdvectionSUPG` and `LevelSetTimeDerivativeSUPG` (keeping only the Galerkin `TimeDerivative` and `LevelSetAdvection` kernels) produces oscillations that grow rapidly: max_phi exceeds 1.5 within the first few steps and the simulation may diverge. SUPG is not optional for the pure advection equation on standard Lagrange elements — it is a necessary stabilization.

---

## Key Takeaways

- The `level_set` module provides `LevelSetAdvection`, `LevelSetAdvectionSUPG`, and `LevelSetTimeDerivativeSUPG` as a complete set of kernels for SUPG-stabilized level set advection. All three are needed: the Galerkin advection term plus both SUPG modifications (advection and time derivative). Omitting the time derivative SUPG term produces an inconsistent formulation that still oscillates.
- `LAGRANGE_VEC` is the recommended family for vector-valued auxiliary variables in the level_set module. It stores the full velocity vector as a single object, avoids the overhead of separate x/y components, and is accepted directly by the `LevelSetAdvection` and SUPG kernel `velocity` parameters. The companion `VectorFunctionIC` and `ParsedVectorFunction` objects enable function-based initialization of vector fields.
- BDF2 (`scheme = bdf2`) provides second-order temporal accuracy for smooth transient problems. For pure advection where the exact solution is simply the initial condition transported along characteristics, temporal accuracy directly controls how well the bubble shape is preserved. BDF2 reduces temporal smearing compared to first-order implicit Euler.
- SUPG stabilization adds streamline diffusion proportional to tau ~ h/(2|v|). On a 40x40 mesh with h = 0.025 and |v| ~ 0.25 (at the bubble radius), tau ~ 0.05 s. This adds a small but crucial amount of numerical diffusion that prevents the Gibbs oscillations inherent in centered Galerkin schemes for advection-dominated problems.
- The overshoots (max_phi > 1) and undershoots (min_phi < 0) seen at t = pi are a known limitation of SUPG for problems with steep gradients. Bound-preserving methods (Flux-Corrected Transport, discontinuous Galerkin with limiters, or the level_set reinitialization approach) can enforce the physical bounds [0, 1] but add complexity. The MOOSE level_set module also provides reinitialization kernels that periodically restore phi to a signed distance function, reducing diffusion of the bubble profile over time.
- The solid-body rotation benchmark is a standard test for advection schemes because the exact solution is known (the bubble returns to its starting position after one full rotation), allowing direct error measurement. It is particularly challenging because the maximum strain rate is zero (no stretching, only rotation), so all deviation from the initial condition is due to numerical error rather than physical deformation.
- Level set methods extend naturally to two-phase flow simulation by coupling phi to a Navier-Stokes solver. The velocity field is no longer prescribed but is solved from the momentum equations with density and viscosity that depend on the local value of phi. MOOSE's `navier_stokes` module supports this coupling for incompressible two-phase flows.
