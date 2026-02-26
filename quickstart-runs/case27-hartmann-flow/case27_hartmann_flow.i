# ============================================================
# Case 27: MHD Hartmann Flow — Magnetic Braking of Channel Flow
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 9 §9.9-9.10
#
# Pressure-driven channel flow with a uniform transverse magnetic
# field B₀. The Lorentz force j × B opposes the flow, flattening
# the velocity profile from parabolic (Poiseuille) toward the
# Hartmann profile with thin boundary layers.
#
# Steady momentum balance:
#   0 = -dp/dx + μ·∂²v_x/∂y² - σB₀²·v_x
#
# The Lorentz drag -σB₀²·v is mathematically identical to Darcy
# friction in a porous medium: -μ/(ρ·K)·v.  We exploit this by
# using porous_medium_treatment with K = μ/(ρ·σB₀²) = 1/Ha².
#
# Hartmann number:
#   Ha = B₀·L·√(σ/μ) = 5
#
# Analytical Hartmann profile (y ∈ [0,1], walls at y=0 and y=1):
#   v(y) = v_max · [1 - cosh(Ha·(y-0.5)) / cosh(Ha/2)]
#
# Non-dimensional parameters:
#   ρ = 1, μ = 1, Ha = 5 → K = 1/Ha² = 0.04
#   Inlet velocity: v_in = 1.0 (flow enters from left)
#
# Domain: [0,5] × [0,1], 50×25 channel mesh
# BCs: no-slip top/bottom walls, fixed velocity inlet, pressure outlet
# ============================================================

Ha2 = 25.0   # Hartmann number squared = 1/K_perm

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 1
    nx   = 50
    ny   = 25
  []
[]

# NavierStokesFV with porous_medium_treatment = true adds Darcy friction:
#   f_drag = -μ/(ρ·K) · v = -(1/(1·0.04)) · v = -25 · v
# This is exactly the Lorentz drag -σB₀²·v for Ha² = 25.
# The pressure gradient driving the flow is established naturally by
# the inlet velocity and outlet pressure boundary conditions.
[Modules]
  [NavierStokesFV]
    compressibility         = 'incompressible'
    porous_medium_treatment = true
    add_energy_equation     = false

    density           = 'rho'
    dynamic_viscosity = 'mu'

    initial_velocity = '1 0 0'
    initial_pressure = 0.0

    # Top and bottom walls: no-slip (Hartmann walls)
    wall_boundaries     = 'top bottom'
    momentum_wall_types = 'noslip noslip'

    # Left inlet: fixed velocity = 1.0 in x-direction
    inlet_boundaries       = 'left'
    momentum_inlet_types   = 'fixed-velocity'
    momentum_inlet_functors = '1 0'

    # Right outlet: zero pressure
    outlet_boundaries      = 'right'
    momentum_outlet_types  = 'fixed-pressure'
    pressure_functors      = '0'

    # Darcy friction for the Lorentz drag
    friction_types  = 'darcy'
    friction_coeffs = 'Darcy_coefficient'

    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
  []
[]

[FunctorMaterials]
  [fluid_props]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho mu'
    prop_values = '1.0 1.0'
  []
  # Darcy coefficient vector: 1/K for each direction.
  # K = 0.04 → 1/K = 25 in all directions (isotropic drag).
  [darcy_coeff]
    type        = ADGenericVectorFunctorMaterial
    prop_names  = 'Darcy_coefficient'
    prop_values = '${Ha2} ${Ha2} ${Ha2}'
  []
  # Porosity = 1.0 (not a real porous medium, just using the drag term)
  [porosity]
    type        = ADGenericFunctorMaterial
    prop_names  = 'porosity'
    prop_values = '1.0'
  []
[]

[Postprocessors]
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = superficial_vel_x
  []
  [avg_vel_x]
    type     = ElementAverageValue
    variable = superficial_vel_x
  []
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = superficial_vel_y
  []
  [avg_pressure]
    type     = ElementAverageValue
    variable = pressure
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
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
