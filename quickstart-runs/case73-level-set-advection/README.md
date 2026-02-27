# Case 73: Level Set Bubble Advection — SUPG Stabilization

## Overview

This case demonstrates the **level_set module** for interface tracking using the
level set method with **SUPG (Streamline Upwind Petrov-Galerkin) stabilization**.
A smooth circular bubble defined by a level set function is advected by a
solid-body rotation velocity field — a standard benchmark for advection schemes.

The level set equation is a pure advection PDE, which is notoriously difficult to
solve with standard Galerkin FEM due to oscillations. SUPG stabilization adds
streamline diffusion to suppress these oscillations while maintaining accuracy
along characteristics.

**New concepts introduced**: `LevelSetAdvection`, `LevelSetAdvectionSUPG`,
`LevelSetTimeDerivativeSUPG`, `LAGRANGE_VEC` vector variables,
`ParsedVectorFunction`, `VectorFunctionIC`, mass conservation monitoring.

---

## The Physics

**Governing equation**: Level set advection (pure transport)

    dphi/dt + v . nabla(phi) = 0    in [0, 1]^2

**Velocity field**: Solid-body rotation about (0.5, 0.5)

    v_x = -(y - 0.5)
    v_y =  (x - 0.5)

This velocity field rotates the entire domain counterclockwise with period 2*pi.
Every point traces a circle about the center.

**Initial condition**: Smooth bubble centered at (0.5, 0.75) with radius 0.15

    phi(x, y, 0) = 0.5 * (1 + tanh((0.15 - r) / 0.02))

where r = sqrt((x-0.5)^2 + (y-0.75)^2). The tanh profile creates a smooth
transition from phi = 1 (inside) to phi = 0 (outside) over a width of ~0.02.

**Expected behavior**: The bubble rotates counterclockwise. After a half-rotation
(t = pi ~ 3.14), the bubble should be at (0.5, 0.25). With perfect advection,
the bubble shape and peak value would be preserved exactly.

**Conservation**: The integral of phi over the domain should remain constant
(the level set equation is conservative).

---

## Input File Walkthrough

### [Mesh]
40x40 generated mesh on [0, 1]^2. Finer than many previous cases because
advection problems require good resolution to minimize numerical diffusion.

### [Variables] / [AuxVariables]
- `phi`: The level set function (primary variable, LAGRANGE)
- `velocity`: The advection velocity field (auxiliary variable, `LAGRANGE_VEC` —
  a vector-valued nodal variable storing both components)

### [ICs]
- `phi_ic`: Smooth tanh bubble profile via `FunctionIC`
- `vel_ic`: Solid-body rotation via `VectorFunctionIC` with a `ParsedVectorFunction`

### [Functions]
- `phi_init`: Smooth bubble `0.5*(1 + tanh((0.15 - r)/0.02))`
- `vel_func`: Rotation field `v_x = -(y-0.5)`, `v_y = (x-0.5)`

### [Kernels]
Four kernels working together:

1. **`TimeDerivative`**: Standard dphi/dt term
2. **`LevelSetAdvection`**: Standard Galerkin advection v . nabla(phi)
3. **`LevelSetAdvectionSUPG`**: SUPG stabilization for the advection term — adds
   streamline diffusion proportional to the element size and velocity magnitude
4. **`LevelSetTimeDerivativeSUPG`**: SUPG stabilization for the time derivative
   term — ensures consistency of the SUPG formulation

All four kernels take `velocity` as a parameter. The SUPG kernels use the velocity
to compute the streamline direction and stabilization parameter tau.

### [Postprocessors]
- `total_phi`: Integral of phi (mass conservation check — should be constant)
- `max_phi`: Peak value (ideally stays at 1.0; numerical diffusion/oscillation
  causes deviations)
- `min_phi`: Minimum value (ideally stays at 0.0; oscillations cause undershoots)

### [Executioner]
BDF2 time integration (second-order accurate), NEWTON solver with LU
preconditioner. dt = 0.1 gives ~31 steps for the half-rotation.

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

---

## Expected Results

The simulation converges in ~32 time steps (31 steps of dt = 0.1, plus one smaller
step to reach t = pi exactly). Each time step converges in 1 nonlinear iteration.

**Final values** (t = pi):
- total_phi ~ 0.0718 (initial ~ 0.0715, conserved within ~2.5%)
- max_phi ~ 1.108 (11% overshoot from the ideal 1.0)
- min_phi ~ -0.135 (13.5% undershoot from the ideal 0.0)

**Physical interpretation**:
- The bubble completes a half-rotation from (0.5, 0.75) to (0.5, 0.25)
- SUPG stabilization prevents catastrophic oscillations but some overshoot/undershoot
  remains — this is inherent to the level set method without reinitialization
- Mass is conserved within ~2.5%, which is acceptable for this mesh resolution
- Finer meshes or reinitialization strategies would improve both mass conservation
  and extrema preservation

**Comparison without SUPG**: Without the SUPG kernels (removing `LevelSetAdvectionSUPG`
and `LevelSetTimeDerivativeSUPG`), the simulation would show severe oscillations,
with max_phi >> 1 and min_phi << 0, potentially causing solver divergence.

---

## Key Takeaways

- **Level set method**: Represents interfaces implicitly as the zero contour of a
  scalar field phi. The interface evolves by solving the advection equation.
- **SUPG stabilization is essential**: Pure Galerkin discretization of advection-
  dominated PDEs produces spurious oscillations. SUPG adds just enough streamline
  diffusion to stabilize without excessive smearing.
- **Three SUPG kernels**: The level_set module provides `LevelSetAdvection` (Galerkin),
  `LevelSetAdvectionSUPG` (SUPG for advection), and `LevelSetTimeDerivativeSUPG`
  (SUPG for time derivative). All three are needed for a consistent formulation.
- **Vector auxiliary variables**: `LAGRANGE_VEC` family stores multi-component
  vector fields as a single variable, used with `VectorFunctionIC` and
  `ParsedVectorFunction`.
- **Mass conservation monitoring**: Level set methods are not inherently mass-
  conserving. Tracking `total_phi` over time quantifies the numerical mass
  loss/gain. Reinitialization (available in the module) can improve this.
- These skills extend directly to two-phase flow interface tracking when combined
  with the `navier_stokes` module.
