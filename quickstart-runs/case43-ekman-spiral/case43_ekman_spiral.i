# ============================================================
# Case 43: Ekman Spiral — Rotating Boundary Layer
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 8, Sec 8.4
#
# Steady viscous flow near a wall in a rotating frame. The
# Coriolis force couples the x and y velocity components,
# producing the famous Ekman spiral — the velocity vector
# rotates through the boundary layer.
#
# Governing equations (1D in z, steady):
#   nu * d²vx/dz² + 2*Omega*vy = 0
#   nu * d²vy/dz² - 2*Omega*vx = -2*Omega*U_g
#
# Rewritten for MOOSE (move coupling to RHS as source terms):
#   nu * d²vx/dz² = -2*Omega*vy           (CoupledForce)
#   nu * d²vy/dz² = +2*Omega*vx - 2*Omega*U_g  (CoupledForce + BodyForce)
#
# Analytical solution:
#   vx(z) = U_g * (1 - exp(-z/delta_E) * cos(z/delta_E))
#   vy(z) = U_g * exp(-z/delta_E) * sin(z/delta_E)
# where delta_E = sqrt(nu/Omega) is the Ekman layer thickness.
#
# Parameters: nu=0.01, Omega=1.0, U_g=1.0
#   delta_E = sqrt(0.01/1.0) = 0.1
#
# Validation:
#   vx(delta_E) = 1*(1 - e^{-1}*cos(1)) = 1 - 0.1988 = 0.801
#   vy(delta_E) = 1*e^{-1}*sin(1) = 0.310
#   max(vy) ~ 0.322 at z ~ 0.0785 (= pi/4 * delta_E)
#
# Domain: [0, 0.7] x [0, 0.01], 200x1 elements (quasi-1D in z)
# Note: x-coordinate represents the vertical z-direction
# ============================================================

nu = 0.01
Omega = 1.0
U_g = 1.0

# Coriolis coupling coefficient: 2*Omega = 2.0
# Sign derivation for CoupledForce:
#   CoupledForce::computeQpResidual returns -coef * v[qp]
#   So the kernel adds -coef * v to the residual for variable u.
#   ADMatDiffusion (strong form): contributes -nu * d²u/dz² to the LHS residual.
#
# For the vx equation, the full residual is:
#   R_vx = -nu * d²vx/dz² + (-coef_vx) * vy = 0
# We want:  -nu * d²vx/dz² - 2*Omega*vy = 0
#   => -coef_vx = -2*Omega  =>  coef_vx = 2*Omega  (positive)
#
# For the vy equation, the full residual is:
#   R_vy = -nu * d²vy/dz² + (-coef_vy) * vx + (-body_val) = 0
# We want:  -nu * d²vy/dz² + 2*Omega*vx - 2*Omega*U_g = 0
#   => -coef_vy = +2*Omega  =>  coef_vy = -2*Omega  (negative)
#   => -body_val = -2*Omega*U_g  =>  body_val = 2*Omega*U_g  (positive)

coriolis = '${fparse 2.0 * Omega}'
body_force_val = '${fparse 2.0 * Omega * U_g}'

[Mesh]
  type = GeneratedMesh
  dim = 2
  # x represents z (vertical coordinate through the Ekman boundary layer)
  # y is a dummy dimension for the quasi-1D geometry
  xmin = 0
  xmax = 0.7
  ymin = 0
  ymax = 0.01
  nx = 200
  ny = 1
[]

[Variables]
  [vx]
  []
  [vy]
  []
[]

[Kernels]
  # --- Equation for vx: nu*d²vx/dz² + 2*Omega*vy = 0 ---
  [vx_diff]
    type = ADMatDiffusion
    variable = vx
    diffusivity = nu_mat
  []
  # Coriolis coupling from vy: coef = 2*Omega (see sign derivation above)
  [vx_coriolis]
    type = CoupledForce
    variable = vx
    v = vy
    coef = ${coriolis}
  []

  # --- Equation for vy: nu*d²vy/dz² - 2*Omega*vx + 2*Omega*U_g = 0 ---
  [vy_diff]
    type = ADMatDiffusion
    variable = vy
    diffusivity = nu_mat
  []
  # Coriolis coupling from vx: coef = -2*Omega (see sign derivation above)
  [vy_coriolis]
    type = CoupledForce
    variable = vy
    v = vx
    coef = '${fparse -2.0 * Omega}'
  []
  # Geostrophic pressure gradient forcing: body_force_val = 2*Omega*U_g
  [vy_pressure]
    type = BodyForce
    variable = vy
    value = ${body_force_val}
  []
[]

[BCs]
  # Wall at z=0 (x=0): no-slip, vx=0, vy=0
  [vx_wall]
    type = DirichletBC
    variable = vx
    boundary = left
    value = 0
  []
  [vy_wall]
    type = DirichletBC
    variable = vy
    boundary = left
    value = 0
  []
  # Far field at z=0.7 (x=0.7): geostrophic flow vx=U_g, vy=0
  [vx_farfield]
    type = DirichletBC
    variable = vx
    boundary = right
    value = ${U_g}
  []
  [vy_farfield]
    type = DirichletBC
    variable = vy
    boundary = right
    value = 0
  []
[]

[Materials]
  [diffusivity]
    type = ADGenericConstantMaterial
    prop_names = 'nu_mat'
    prop_values = '${nu}'
  []
[]

[Postprocessors]
  # vx at z = delta_E = 0.1 (analytical: 0.801)
  [vx_at_delta]
    type = PointValue
    variable = vx
    point = '0.1 0.005 0'
  []
  # vy at z = delta_E = 0.1 (analytical: 0.310)
  [vy_at_delta]
    type = PointValue
    variable = vy
    point = '0.1 0.005 0'
  []
  # Maximum vy over domain (analytical: ~0.322 at z ~ 0.0785 = pi/4 * delta_E)
  [max_vy]
    type = ElementExtremeValue
    variable = vy
    value_type = max
  []
  # Average vx — tracks approach to geostrophic value U_g=1
  [avg_vx]
    type = ElementAverageValue
    variable = vx
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 20
[]

[Outputs]
  exodus = true
  csv = true
[]
