# ============================================================
# Case 39: Blasius Boundary Layer — Flat Plate Laminar Flow
# Rieutord, Fluid Dynamics (Springer, 2015), Ch. 4, Sec. 4.3
#
# Laminar boundary layer growing over a semi-infinite flat plate.
# Uniform free-stream flow U=1 enters from the left; viscous
# effects near the no-slip bottom wall create a growing shear
# layer whose thickness scales as:
#
#   delta_99(x) = 5.0 * sqrt(mu * x / (rho * U))
#
# The velocity profile u/U as a function of eta = y / delta_99
# collapses onto the universal Blasius similarity solution
# f'(eta), where f satisfies:
#
#   2 f''' + f f'' = 0,   f(0)=0, f'(0)=0, f'(inf)=1
#
# Key results from the Blasius solution:
#   f''(0) = 0.332   (determines wall shear stress)
#   Skin friction:  C_f(x)  = 0.664 / sqrt(Re_x)
#   Displacement thickness:  delta_1 = 1.72 * sqrt(mu*x/(rho*U))
#   Momentum thickness:      delta_2 = 0.664 * sqrt(mu*x/(rho*U))
#
# Parameters (chosen for Re_L = 400, well inside laminar regime):
#   rho = 1.0,  U = 1.0,  mu = 0.005
#   Re_L = rho * U * L / mu = 1 * 1 * 2 / 0.005 = 400
#
# Boundary layer thickness at x=1: delta_99 = 5*sqrt(0.005*1) = 0.354
# Boundary layer thickness at x=2: delta_99 = 5*sqrt(0.005*2) = 0.500
#
# Domain: [0, 2] x [0, 1.5], 60x30 cells with bias_y=0.5 (wall clustering)
# BCs:
#   left   -- fixed-velocity inlet: (u,v) = (1, 0)
#   bottom -- no-slip flat plate: (u,v) = (0, 0)
#   top    -- slip (symmetry): models free-stream far field
#   right  -- fixed-pressure outlet: p = 0
# ============================================================

mu_val  = 0.005   # dynamic viscosity  [Pa s]
rho_val = 1.0     # density            [kg/m^3]
U_inf   = 1.0     # free-stream speed  [m/s]

[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1.5
    nx    = 60      # streamwise resolution — captures developing velocity profile
    ny    = 30      # wall-normal resolution with clustering near bottom wall
    bias_y = 0.5   # bias < 1 clusters cells toward the bottom (flat plate side)
  []
[]

# The NSFVAction constructs all FV variables (vel_x, vel_y, pressure),
# kernels, and boundary conditions for incompressible Navier-Stokes.
[Modules]
  [NavierStokesFV]
    compressibility         = 'incompressible'
    porous_medium_treatment = false
    add_energy_equation     = false

    # Material property names — values defined in [FunctorMaterials] below
    density           = 'rho'
    dynamic_viscosity = 'mu'

    # Initial condition: uniform free-stream throughout the domain.
    # Newton converges from this starting point for Re_L = 400.
    initial_velocity = '${U_inf} 0 0'
    initial_pressure = 0.0

    # Inlet: uniform flow U=1 from the left boundary
    inlet_boundaries        = 'left'
    momentum_inlet_types    = 'fixed-velocity'
    momentum_inlet_functors = '${U_inf} 0'

    # Outlet: zero pressure on the right boundary (Dirichlet p=0)
    # This fixes the pressure level, so pressure pinning is not needed.
    outlet_boundaries       = 'right'
    momentum_outlet_types   = 'fixed-pressure'
    pressure_functors       = '0'

    # Bottom: no-slip flat plate — the surface that generates the boundary layer
    # Top:    slip (symmetry) — approximates the undisturbed free stream at y=1.5
    wall_boundaries       = 'bottom top'
    momentum_wall_types   = 'noslip slip'

    # Upwind advection: more robust for the developing boundary layer where
    # the local Re in cells near the inlet can be moderately high.
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
  []
[]

[FunctorMaterials]
  [fluid_props]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu'
    prop_values = '${rho_val}  ${mu_val}'
    # Re_L = rho * U * L / mu = 1 * 1 * 2 / 0.005 = 400
    # Blasius solution valid for Re_x < ~5e5 (laminar regime)
  []
[]

[Postprocessors]
  # Free-stream velocity (should remain near 1.0 away from the plate)
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []

  # Maximum wall-normal velocity (induced by boundary layer displacement)
  # Blasius theory predicts v/U ~ 0.86 * sqrt(mu/(rho*U*x)) at large y
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []

  # Domain-average streamwise velocity — decreases from 1.0 as the
  # boundary layer displaces mass toward the free stream
  [avg_vel_x]
    type     = ElementAverageValue
    variable = vel_x
  []

  # Sample u-velocity at x=1 (Re_x=200), midway through the domain.
  # Blasius: u(x=1, y=delta_99) / U should be ~0.99; at y=delta_99/2 ~0.63.
  # delta_99 at x=1 = 5*sqrt(0.005*1/1) = 0.354
  [u_at_x1_y035]
    type     = PointValue
    variable = vel_x
    point    = '1.0 0.35 0'   # approximately y = delta_99 at x=1
  []

  [u_at_x1_y010]
    type     = PointValue
    variable = vel_x
    point    = '1.0 0.10 0'   # inside the boundary layer at x=1
  []

  # Sample u-velocity at x=2 (Re_x=400), at the plate exit.
  # delta_99 at x=2 = 5*sqrt(0.005*2/1) = 0.500
  [u_at_x2_y050]
    type     = PointValue
    variable = vel_x
    point    = '2.0 0.50 0'   # approximately y = delta_99 at x=2
  []

  [u_at_x2_y020]
    type     = PointValue
    variable = vel_x
    point    = '2.0 0.20 0'   # inside the boundary layer at x=2
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'

  # Direct LU factorization: robust for the incompressible saddle-point
  # system at Re_L = 400. The NONZERO shift prevents factorization failure
  # from zeros on the pressure diagonal.
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30

  # Automatic scaling balances velocity and pressure degrees of freedom,
  # which improves Newton convergence for advection-dominated flows.
  automatic_scaling = true
[]

[Outputs]
  exodus = true   # vel_x, vel_y, pressure fields for ParaView visualization
  csv    = true   # postprocessor values for comparison with Blasius theory
[]
