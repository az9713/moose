# ============================================================
# Case 15: Lid-Driven Cavity — Incompressible Navier-Stokes
# Classic CFD benchmark: Re = rho*U*L/mu = 1*1*1/0.01 = 100
# Steady state, 2D square cavity, top wall moves at U=1
# ============================================================
#
# Domain: unit square [0,1] x [0,1]
# Boundary conditions:
#   top    -- moving lid: vel_x = 1, vel_y = 0
#   left   -- no-slip:   vel = 0
#   right  -- no-slip:   vel = 0
#   bottom -- no-slip:   vel = 0
# Fluid properties:
#   rho = 1    (density)
#   mu  = 0.01 (dynamic viscosity)
#   Re  = rho * U * L / mu = 1 * 1 * 1 / 0.01 = 100
# ============================================================

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx   = 30
    ny   = 30
  []
[]

# The NSFVAction automatically creates all FV variables, kernels,
# and boundary conditions needed for incompressible Navier-Stokes.
# It uses cell-centered finite volume (CCFV) discretization.
[Modules]
  [NavierStokesFV]
    compressibility          = 'incompressible'
    porous_medium_treatment  = false
    add_energy_equation      = false

    # Material property names for density and dynamic viscosity.
    # Actual values are defined in the [FunctorMaterials] block below.
    density          = 'rho'
    dynamic_viscosity = 'mu'

    # Initial conditions for the velocity and pressure fields.
    # A quiescent start is standard for lid-driven cavity.
    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    # No-slip walls: left, right, and bottom boundaries.
    wall_boundaries   = 'left right bottom'
    momentum_wall_types = 'noslip noslip noslip'

    # The lid (top boundary) is treated as a fixed-velocity inlet.
    # This is the standard MOOSE FV approach for a moving wall:
    # prescribe (vel_x, vel_y) = (1, 0) via momentum_inlet_function.
    inlet_boundaries        = 'top'
    momentum_inlet_types    = 'fixed-velocity'
    momentum_inlet_functors = '1 0'

    # Advection interpolation scheme.
    # 'average' (central difference) is stable at Re=100 on a 30x30 mesh.
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'average'

    # Pin the domain-average pressure to zero to fix the pressure level.
    # Incompressible flow determines pressure only up to an additive constant;
    # pinning removes this null space and ensures a unique solution.
    pin_pressure       = true
    pinned_pressure_type  = average
    pinned_pressure_value = 0
  []
[]

# ADGenericFunctorMaterial declares functor (cell-averaged) material
# properties compatible with the FV Navier-Stokes action.
[FunctorMaterials]
  [fluid_properties]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho mu'
    prop_values = '1   0.01'    # Re = rho*U*L/mu = 1*1*1/0.01 = 100
  []
[]

[Postprocessors]
  # Maximum velocity magnitude anywhere in the domain.
  # At Re=100 the peak occurs inside the primary vortex, slightly
  # below the centerline toward the bottom-right of the cavity.
  [max_vel_x]
    type       = ElementExtremeValue
    variable   = vel_x
    value_type = max
  []

  [min_vel_x]
    type       = ElementExtremeValue
    variable   = vel_x
    value_type = min
  []

  [max_vel_y]
    type       = ElementExtremeValue
    variable   = vel_y
    value_type = max
  []

  [min_vel_y]
    type       = ElementExtremeValue
    variable   = vel_y
    value_type = min
  []

  # Domain-average pressure (should remain near zero due to pinning).
  [avg_pressure]
    type     = ElementAverageValue
    variable = pressure
  []
[]

[Executioner]
  type = Steady

  # NEWTON is preferred over PJFNK for the saddle-point structure of
  # incompressible N-S. LU factorization (via MUMPS or SuperLU_dist)
  # provides a robust direct solve for a 30x30 FV system.
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Tight tolerances are achievable at Re=100 with Newton + LU.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50

  l_tol     = 1e-6
  l_max_its = 200
[]

[Outputs]
  exodus = true
  csv    = true
[]
