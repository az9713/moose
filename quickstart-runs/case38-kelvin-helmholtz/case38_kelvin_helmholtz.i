# ============================================================
# Case 38: Kelvin-Helmholtz Instability — Shear Layer Rollup
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 6, Sec 6.3.1
#
# Two counterflowing streams separated by a shear layer become
# unstable to Kelvin-Helmholtz billows. The velocity profile
# is a tanh shear layer, and a passive scalar (temperature
# used as dye) marks the two streams.
#
# A small sinusoidal perturbation in vel_y seeds the instability.
# The shear layer rolls up into vortex structures that grow and
# eventually merge.
#
# Velocity profile: vel_x = tanh((y-0.5)/0.05)
#   Upper stream (y>0.5): vel_x → +1
#   Lower stream (y<0.5): vel_x → -1
#
# Passive scalar: T = 0.5*(1 + tanh((y-0.5)/0.05))
# Perturbation:   vel_y = 0.01*sin(2*pi*x)
#
# Parameters: rho=1, mu=0.001 (Re~1000), kappa=0.001
# KH growth rate: sigma_max ~ 0.2*dU/delta ~ 0.2*2/0.05 = 8
#
# Domain: [0, 2] x [0, 1], 80x20 = 1600 elements
#   accommodates 2 full wavelengths along x (lambda = 1 each)
# BCs: inlet/outlet (left/right) to approximate periodic,
#       top/bottom: slip walls (no normal flow, no friction)
# ============================================================

mu_val = 0.001
k_val  = 0.001

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 1
    nx   = 80
    ny   = 20
  []
[]

# NavierStokesFV action sets up all FV momentum, pressure, and energy kernels.
# The energy equation is repurposed as a passive scalar transport: cp=1, k=kappa,
# no buoyancy. Temperature T tracks the dye concentration of each stream.
[Modules]
  [NavierStokesFV]
    compressibility         = 'incompressible'
    porous_medium_treatment = false
    add_energy_equation     = true

    # Functor material property names (defined in [FunctorMaterials] below)
    density              = 'rho'
    dynamic_viscosity    = 'mu'
    thermal_conductivity = 'k'
    specific_heat        = 'cp'

    # Initial conditions: use function names for the tanh shear-layer profile.
    # The NSFV action creates FVFunctionICs from these function references.
    initial_velocity    = 'vel_x_init_func vel_y_init_func 0'
    initial_pressure    = 0.0
    initial_temperature = T_init_func

    # Left inlet: prescribe the tanh velocity profile and the corresponding
    # passive-scalar value. This continuously injects the shear layer and keeps
    # the left edge consistent with the interior initial condition.
    inlet_boundaries        = 'left'
    momentum_inlet_types    = 'fixed-velocity'
    momentum_inlet_functors = 'vel_x_inlet 0'
    energy_inlet_types      = 'fixed-temperature'
    energy_inlet_functors   = 'T_inlet'

    # Right outlet: fixed pressure = 0 lets fluid exit freely.
    outlet_boundaries      = 'right'
    momentum_outlet_types  = 'fixed-pressure'
    pressure_functors      = '0'

    # Top and bottom: slip walls (symmetry planes).
    # No normal velocity (enforced by MOOSE), no tangential friction.
    # Adiabatic (zero heat flux) so the scalar does not diffuse out the top/bottom.
    wall_boundaries        = 'top bottom'
    momentum_wall_types    = 'slip slip'
    energy_wall_types      = 'heatflux heatflux'
    energy_wall_functors   = '0 0'

    # Upwind for momentum and energy suppresses spurious oscillations near the
    # steep tanh shear layer. Average for mass conservation (continuity).
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
    energy_advection_interpolation   = 'upwind'
  []
[]

[Functions]
  # Inlet velocity: tanh shear layer — right-going above y=0.5, left-going below.
  # The shear layer half-thickness is delta = 0.05 (resolved by 20 cells over
  # a height of 1, giving dy = 0.05 per cell — just at the resolution limit).
  [vel_x_inlet]
    type       = ParsedFunction
    expression = 'tanh((y-0.5)/0.05)'
  []

  # Inlet passive scalar: 1 in upper stream (y>0.5), 0 in lower stream (y<0.5).
  # The smooth tanh transition matches the initial condition to avoid a sharp
  # discontinuity that would cause oscillations at the inlet face.
  [T_inlet]
    type       = ParsedFunction
    expression = '0.5*(1 + tanh((y-0.5)/0.05))'
  []

  # Initial velocity profile: same tanh everywhere (not just at inlet).
  [vel_x_init_func]
    type       = ParsedFunction
    expression = 'tanh((y-0.5)/0.05)'
  []

  # Initial perturbation: sinusoidal in x with amplitude 0.01.
  # Two full wavelengths fit in [0, 2], seeding two KH billows.
  # The amplitude 0.01 << dU = 2 keeps the perturbation in the linear
  # growth regime initially; nonlinear rollup develops around t ~ 0.5-1.
  [vel_y_init_func]
    type       = ParsedFunction
    expression = '0.01*sin(2*pi*x)'
  []

  # Initial temperature: same tanh scalar profile as the inlet.
  [T_init_func]
    type       = ParsedFunction
    expression = '0.5*(1 + tanh((y-0.5)/0.05))'
  []
[]

[FunctorMaterials]
  [const]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu        k        cp'
    prop_values = '1.0  ${mu_val} ${k_val} 1.0'
    # Re = rho * dU * delta / mu = 1 * 2 * 0.05 / 0.001 = 100 (shear-layer Re)
    # Full-domain Re = rho * dU * L / mu = 1 * 2 * 1 / 0.001 = 2000
    # Pr = mu*cp/k = 0.001*1/0.001 = 1  (scalar diffuses at same rate as momentum)
  []
[]

[Postprocessors]
  # Velocity extremes track the evolving shear layer.
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []
  [min_vel_x]
    type       = ADElementExtremeFunctorValue
    functor    = vel_x
    value_type = min
  []

  # max_vel_y and min_vel_y bracket the perturbation amplitude.
  # During the linear growth phase (t < ~0.5), max_vel_y grows exponentially
  # at rate sigma_max ~ 8. By t ~ 0.3 it should reach ~ 0.01 * exp(8*0.3) ~ 0.1.
  # Saturation into the nonlinear regime is visible when growth slows.
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []
  [min_vel_y]
    type       = ADElementExtremeFunctorValue
    functor    = vel_y
    value_type = min
  []

  # Average temperature remains near 0.5 throughout: the scalar is merely
  # rearranged by the vortices, not created or destroyed.
  [avg_T]
    type     = ElementAverageValue
    variable = T_fluid
  []
[]

[Executioner]
  type       = Transient
  solve_type = 'NEWTON'

  # Direct LU factorization handles the coupled saddle-point system robustly.
  # At Re~1000 and on a coarse mesh, the momentum equations are not stiff enough
  # to require a specialized block preconditioner for this element count.
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20

  dt       = 0.01
  end_time = 2.0

  # Automatic scaling balances the very different scales of pressure (~1),
  # velocity (~1), and temperature (~0.5-1). Essential for Newton convergence
  # with the coupled momentum+energy system.
  automatic_scaling = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 0.01
    growth_factor      = 1.2
    cutback_factor     = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true   # full 2D field at every output step for ParaView animation
  csv    = true   # postprocessor time series for instability growth tracking
[]
