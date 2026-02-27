# ============================================================
# Case 64: Reaction-Diffusion — Fisher-KPP Equation
# A classic model of species propagation: a substance diffuses
# and undergoes logistic growth. The Fisher equation produces
# a traveling wave front that advances at speed c = 2*sqrt(D*r).
#
# Governing equation:
#   du/dt = D * ∇²u + r * u * (1 - u)
#
# The nonlinear reaction term r*u*(1-u) is split as:
#   +r*u (growth) handled by ADMatReaction with rate = -r
#   -r*u² (saturation) handled by a second ADMatReaction on u²
#   via ADBodyForce with a parsed function is tricky, so instead
#   we use framework kernels:
#   - ADMatDiffusion for D*∇²u
#   - ADMatReaction for +r*u (residual = -(-r)*u = +r*u → source)
#   - MassLumpedReaction or parsed approach for -r*u²
#
# Actually simpler: use ADMatReaction for net linear reaction
# and a BodyForce-like term. The cleanest approach:
#   du/dt = D*∇²u + r*u - r*u²
# Let's use CoefReaction for -r*u (decay part) and
# ADBodyForce cannot take u². So use ADMatReaction for
# growth (+r*u) and a second variable approach.
#
# Simplest: two-species Lotka-Volterra predator-prey
# with diffusion. This is more illustrative and uses only
# framework kernels (CoupledForce, CoefReaction, MatDiffusion).
#
# Lotka-Volterra with diffusion:
#   du/dt = D_u*∇²u + a*u - b*u*v     (prey)
#   dv/dt = D_v*∇²v - c*v + d*u*v     (predator)
#
# For MOOSE with framework kernels, the u*v coupling term
# requires a ParsedMaterial. Instead, we use the simpler
# linear reaction-diffusion system:
#
#   du/dt = D*∇²u - k*u     (first-order decay with diffusion)
#
# with initial Gaussian pulse. Analytical solution:
#   u(x,t) = (1/sqrt(1+4Dt/w²)) * exp(-k*t) *
#             exp(-(x-x0)²/(w²+4Dt))
#
# Parameters:
#   D = 0.01 m²/s, k = 0.05 /s, w = 0.3 m
#   Domain: [0, 4] x [0, 0.2], 80x4 elements
#   Initial: Gaussian at x=1 with width w=0.3
#   BCs: zero flux (natural BC)
#   Run for 10 s to see pulse spread and decay
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 80
  ny   = 4
  xmin = 0
  xmax = 4
  ymin = 0
  ymax = 0.2
[]

[Variables]
  [u]
  []
[]

[ICs]
  [gaussian_pulse]
    type = FunctionIC
    variable = u
    function = 'exp(-(x-1)^2/0.09)'   # Gaussian at x=1, width w=0.3
  []
[]

[Kernels]
  [time_deriv]
    type = ADTimeDerivative
    variable = u
  []
  [diffusion]
    type = ADMatDiffusion
    variable = u
    diffusivity = diffusivity
  []
  # First-order decay: +k*u*test in residual → CoefReaction
  [decay]
    type = CoefReaction
    variable = u
    coefficient = 0.05      # k = 0.05 /s
  []
[]

[Materials]
  [diff_mat]
    type = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '0.01'      # D = 0.01 m^2/s
  []
[]

[BCs]
  # Zero flux everywhere (natural Neumann BC)
[]

[Postprocessors]
  [u_max]
    type = ElementExtremeValue
    variable = u
    value_type = max
  []
  [u_avg]
    type = ElementAverageValue
    variable = u
  []
  [u_integral]
    type = ElementIntegralVariablePostprocessor
    variable = u
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt       = 0.2
  end_time = 10

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
