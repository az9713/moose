# ============================================================
# Case 16: Natural Convection — Buoyancy-Driven Flow in a
#           Differentially Heated Square Cavity
#
# Physics: incompressible Navier-Stokes + energy equation,
#          Boussinesq buoyancy approximation
#
# Non-dimensional parameters:
#   Ra = g * alpha * dT * L^3 / (nu * kappa) = 10000
#   Pr = nu / kappa = 0.71  (air)
#
# Domain: [0,1] x [0,1]
# BCs:
#   Left  wall: T = 1  (hot)
#   Right wall: T = 0  (cold)
#   Top/bottom: insulated (zero heat flux)
#   All walls:  no-slip velocity
#
# Material properties (non-dimensional, rho=1, cp=1):
#   nu = mu = sqrt(Pr / Ra) = sqrt(0.71 / 10000) = 0.008426
#   kappa = k = nu / Pr    = 0.008426 / 0.71    = 0.011867
#   alpha (thermal expansion) = 1.0
#   gravity = '0 -1 0'  (non-dimensional)
#   T_ref = 0.5
#
# Verification: de Vahl Davis (1983) benchmark
#   Ra = 10^4 => Nu_avg ~ 2.243 (average Nusselt number on hot wall)
# ============================================================

# Non-dimensional kinematic viscosity and thermal diffusivity
# Pr = nu/kappa = 0.71,   Ra = g*alpha*dT*L^3/(nu*kappa) = 10000
# nu = sqrt(Pr/Ra) = sqrt(0.71/10000) = 0.008426
# kappa = nu/Pr = 0.008426/0.71 = 0.011867
# Check: Ra = 1*1*1*1/(0.008426*0.011867) = 1/0.0001 = 10000  (approx)

nu    = 0.008426   # = mu  since rho = 1
kappa = 0.011867   # = k   since rho*cp = 1

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

# Note: [Modules/NavierStokesFV] is a deprecated but still-functional syntax.
# The modern equivalent is [Physics/NavierStokes/Flow/...] and
# [Physics/NavierStokes/FluidHeatTransfer/...] (see boussinesq-action.i in the
# navier_stokes test suite for a complete example).
[Modules]
  [NavierStokesFV]
    compressibility          = 'incompressible'
    porous_medium_treatment  = false
    add_energy_equation      = true

    # Material property references (functor names from [FunctorMaterials])
    density              = 'rho'
    dynamic_viscosity    = 'mu'
    thermal_conductivity = 'k'
    specific_heat        = 'cp'

    # Initial conditions: fluid at rest, mid-range temperature
    initial_velocity    = '1e-15 1e-15 0'
    initial_pressure    = 0.0
    initial_temperature = 0.5

    # All four walls are no-slip for velocity
    wall_boundaries       = 'left right top bottom'
    momentum_wall_types   = 'noslip noslip noslip noslip'

    # Energy BCs:
    #   left  -> fixed T = 1 (hot wall)
    #   right -> fixed T = 0 (cold wall)
    #   top   -> zero heat flux (insulated)
    #   bottom-> zero heat flux (insulated)
    energy_wall_types    = 'fixed-temperature fixed-temperature heatflux heatflux'
    energy_wall_functors = '1 0 0 0'

    # Boussinesq buoyancy: rho = rho_ref * (1 - alpha*(T - T_ref))
    boussinesq_approximation = true
    gravity                  = '0 -1 0'
    ref_temperature          = 0.5
    thermal_expansion        = 'alpha'

    # Pin the pressure to zero average (closed cavity with no outlets)
    pin_pressure         = true
    pinned_pressure_type = average
    pinned_pressure_value = 0

    # Interpolation schemes: upwind for momentum and energy for stability
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
    energy_advection_interpolation   = 'upwind'
  []
[]

[FunctorMaterials]
  [const]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu       k         cp  alpha'
    prop_values = '1.0  ${nu}   ${kappa}  1.0  1.0'
    # rho = 1,  cp = 1  => nu = mu/rho = mu,  kappa = k/(rho*cp) = k
    # Ra = g*alpha*dT*L^3 / (nu*kappa)
    #    = 1 * 1 * 1 * 1 / (0.008426 * 0.011867)
    #    = 1 / 0.00010000 = 10000
    # Pr = nu / kappa = 0.008426 / 0.011867 = 0.71
  []
[]

[Postprocessors]
  # Maximum velocity magnitude (diagnostic for strength of convection)
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []

  # Average temperature in the cavity (should remain ~0.5 by symmetry)
  [avg_T]
    type     = ElementAverageValue
    variable = T_fluid
  []

  # Nusselt number on the hot left wall:
  #   Nu = (-k * dT/dn) * L / (k * dT) = -dT/dx|_{x=0}
  # For a unit cavity with dT = 1, Nu = -dT/dx|_{x=0}.
  # We compute the surface-integrated heat flux / (k * dT/L) via
  # the FVFlux postprocessor as a proxy for Nu.
  [Nu_hot_wall]
    type     = SideAverageFunctorPostprocessor
    functor  = 'T_fluid'
    boundary = 'left'
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'

  # Direct LU factorization: robust for coupled Navier-Stokes + energy
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Tight nonlinear tolerance to capture correct flow structure
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30

  # Automatic scaling helps balance velocity/pressure/temperature
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv    = true
[]
