# ============================================================
# Case 63: Gravity Dam — Self-Weight Loading & Foundation Stress
# A gravity dam resting on a deformable foundation under
# self-weight demonstrates stress distribution in geotechnical
# structures. The dam transmits its weight to the foundation,
# producing a stress concentration at the toe and heel.
#
# This is a plane-strain solid mechanics problem with:
#   - Gravity body force (self-weight)
#   - Two material blocks (dam: concrete, foundation: rock)
#   - Stress analysis showing the load transfer pattern
#
# Parameters:
#   Dam: H=10 m, base=8 m, E=25 GPa, nu=0.2, rho=2400 kg/m^3
#   Foundation: 20 m wide x 10 m deep, E=10 GPa, nu=0.3, rho=2600 kg/m^3
#
# The dam shape is approximated as a rectangular block (simple
# geometry from GeneratedMesh). In reality, dams are trapezoidal
# but the rectangular shape suffices for educational purposes.
#
# Domain: [0, 20] x [-10, 10] m
#   y > 0: dam (centered, x in [6, 14])
#   y < 0: foundation (full width)
# Mesh: 20x20 elements
#
# BCs: bottom fixed, sides roller, top free
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
    xmin = 0
    xmax = 20
    ymin = -10
    ymax = 10
  []
  # Define subdomains: dam (y>0, 6<x<14) vs foundation
  [dam_block]
    type = SubdomainBoundingBoxGenerator
    input = gmg
    block_id = 1
    block_name = dam
    bottom_left = '6 0 0'
    top_right = '14 10 0'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_yy stress_xx vonmises_stress'
  []
[]

[Kernels]
  # Gravity body force in y-direction
  [gravity_y]
    type = Gravity
    variable = disp_y
    value = -9.81     # m/s^2 downward
  []
[]

[Materials]
  # Dam material (concrete)
  [elasticity_dam]
    type = ComputeIsotropicElasticityTensor
    block = dam
    youngs_modulus = 2.5e10    # 25 GPa
    poissons_ratio = 0.2
  []
  [density_dam]
    type = GenericConstantMaterial
    block = dam
    prop_names = 'density'
    prop_values = '2400'       # kg/m^3
  []
  [stress_dam]
    type = ComputeLinearElasticStress
    block = dam
  []

  # Foundation material (rock)
  [elasticity_found]
    type = ComputeIsotropicElasticityTensor
    block = 0
    youngs_modulus = 1e10      # 10 GPa
    poissons_ratio = 0.3
  []
  [density_found]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'density'
    prop_values = '2600'       # kg/m^3
  []
  [stress_found]
    type = ComputeLinearElasticStress
    block = 0
  []
[]

[BCs]
  # Bottom: fully fixed
  [fix_bottom_x]
    type     = DirichletBC
    variable = disp_x
    boundary = bottom
    value    = 0
  []
  [fix_bottom_y]
    type     = DirichletBC
    variable = disp_y
    boundary = bottom
    value    = 0
  []
  # Left and right: roller (no horizontal displacement)
  [fix_left_x]
    type     = DirichletBC
    variable = disp_x
    boundary = left
    value    = 0
  []
  [fix_right_x]
    type     = DirichletBC
    variable = disp_x
    boundary = right
    value    = 0
  []
[]

[Postprocessors]
  # Max vertical stress (should be at dam base / foundation interface)
  [max_stress_yy]
    type = ElementExtremeValue
    variable = stress_yy
  []
  # Max von Mises stress
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
  # Max vertical displacement (settlement)
  [max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
  []
  # Stress at dam centerline base (x=10, y=0)
  [stress_yy_center]
    type = PointValue
    variable = stress_yy
    point = '10 0 0'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv    = true
[]
