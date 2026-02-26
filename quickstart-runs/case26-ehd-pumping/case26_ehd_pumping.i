# ============================================================
# Case 26: EHD Pumping — Coulomb Force Drives Fluid Flow
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 9 §9.11
#
# A prescribed electric body force (Coulomb force) drives
# recirculating flow in a closed cavity. The force mimics
# the effect of space charge in an applied electric field:
#   f_x = A · (1-x) · sin(π·y)
#
# This is the EHD analog of natural convection — instead of
# buoyancy, an electric body force drives the fluid.
#
# Domain: unit square [0,1]², 25×25 mesh
# BCs: no-slip on all walls, pressure pinned to zero average
# Fluid: ρ = 1, μ = 0.01  (Re_eff ~ 50-100)
# Force amplitude: A = 500
# Steady-state solution.
# ============================================================

force_amp = 500.0   # Coulomb force amplitude

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx   = 25
    ny   = 25
  []
[]

# The NSFVAction constructs all FV variables (vel_x, vel_y, pressure),
# all momentum/continuity kernels, Rhie-Chow interpolation, and all wall
# BCs from the single block below.  The EHD body force is then appended
# as an explicit INSFVBodyForce kernel in [FVKernels].
[Modules]
  [NavierStokesFV]
    compressibility         = 'incompressible'
    porous_medium_treatment = false
    add_energy_equation     = false

    # Material property names for density and dynamic viscosity.
    # Actual values are defined in the [FunctorMaterials] block below.
    density           = 'rho'
    dynamic_viscosity = 'mu'

    # Initial conditions: quiescent fluid and zero pressure.
    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    # All four walls are no-slip for velocity.
    wall_boundaries     = 'left right top bottom'
    momentum_wall_types = 'noslip noslip noslip noslip'

    # Advection scheme: upwind for momentum gives better stability when
    # the EHD body force creates strong interior flow gradients.
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'

    # Pin the domain-average pressure to zero.  Incompressible flow in a
    # closed cavity with no outlets determines pressure only up to an
    # additive constant; pinning removes the null space.
    pin_pressure          = true
    pinned_pressure_type  = average
    pinned_pressure_value = 0
  []
[]

# Additional FV kernel for the EHD (Coulomb) body force.
# INSFVBodyForce adds a prescribed body force to the momentum equation.
# The force is evaluated as a functor — here the functor is the
# 'coulomb_force_x' property defined in [FunctorMaterials].
# The rhie_chow_user_object reference is required by INSFVBodyForce when
# using the NSFVAction; the action creates this object with the fixed name
# 'ins_rhie_chow_interpolator'.
[FVKernels]
  [ehd_force_x]
    type                  = INSFVBodyForce
    variable              = vel_x
    functor               = coulomb_force_x
    momentum_component    = 'x'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  []
[]

[FunctorMaterials]
  # Fluid properties: dimensionless, rho = 1, mu = 0.01.
  # With the EHD body force amplitude A = 500, the effective velocity scale
  # is U ~ A * L^2 / mu = 500 * 1 / 0.01 = 50000 (Stokes estimate), giving
  # a very large effective Reynolds number.  In practice inertia limits this,
  # but the problem is strongly convection-driven.
  [fluid_props]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho mu'
    prop_values = '1.0 0.01'
  []

  # Prescribed Coulomb body force: f_x = A * (1 - x) * sin(π·y).
  #
  # Physical interpretation: the force models space charge of density ρ_e
  # acted on by an applied electric field E.  The factor (1 - x) means the
  # force is strongest near the left wall (where charge injection occurs)
  # and decays to zero at the right wall (collecting electrode), consistent
  # with Melcher Ch. 9's picture of unipolar injection across a gap.
  # The sin(π·y) variation drives a single vortex cell: the force is
  # directed rightward in the upper half (sin > 0) and leftward in the lower
  # half (sin < 0), which would normally create two counter-rotating cells;
  # however with the domain going from y=0 to y=1 and sin(π·y) ≥ 0
  # throughout, the force is purely rightward everywhere and the return flow
  # runs along the walls, producing one large recirculation loop.
  #
  # ADParsedFunctorMaterial evaluates an algebraic expression at each cell
  # centre using the cell's spatial coordinates (x, y, z, t are always
  # available as built-in symbols).  The result is an AD-compatible functor
  # property that INSFVBodyForce can consume directly.
  [ehd_force]
    type                  = ParsedFunctorMaterial
    property_name         = coulomb_force_x
    expression            = '${force_amp} * (1 - x) * sin(3.14159265358979 * y)'
    disable_fpoptimizer   = true
    enable_jit            = false
  []
[]

[Postprocessors]
  # Peak x-velocity: positive values confirm the EHD force is driving
  # rightward flow in the interior.
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []

  # Minimum x-velocity (most negative): indicates the magnitude of the
  # return flow along the right and lateral walls.
  [min_vel_x]
    type       = ADElementExtremeFunctorValue
    functor    = vel_x
    value_type = min
  []

  # Peak y-velocity: measures the turning flow at the ends of the cavity.
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []

  # Domain-average pressure: should remain near zero due to pressure pinning.
  [avg_pressure]
    type     = ElementAverageValue
    variable = pressure
  []
[]

[Executioner]
  type = Steady

  # Full Newton with exact AD Jacobians.  The saddle-point velocity-pressure
  # system benefits from full Newton over PJFNK because the coupling between
  # momentum and continuity is captured exactly in the Jacobian.
  solve_type = 'NEWTON'

  # Direct LU factorization: reliable for the 25x25 FV system (~1875 cells,
  # ~5625 DOF for vel_x, vel_y, pressure).  The NONZERO shift prevents zero-
  # pivot failures in the pressure rows of the saddle-point block.
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Tight nonlinear tolerances.  EHD-driven flow at moderate body-force
  # amplitude typically converges in 10-30 Newton iterations from rest.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
  l_tol      = 1e-6
  l_max_its  = 200
[]

[Outputs]
  exodus = true   # Full vel_x, vel_y, pressure fields for ParaView
  csv    = true   # Postprocessor values at convergence
[]
