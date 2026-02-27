# ============================================================
# Case 50: Finite Strain — Large Deformation Compression
# A 2D rubber-like block is compressed 30% in the y-direction.
# Uses FINITE strain formulation with hyperelastic stress
# (ComputeFiniteStrainElasticStress).  The nonlinear geometric
# effects (rotation of material frame) distinguish this from
# small-strain linear elasticity.
#
# Domain: [0,1] x [0,1], 10x10 elements
# Material: Rubber-like — E = 10 MPa, nu = 0.45 (near-incompressible)
# Loading: displacement-controlled compression, disp_y = -0.3 on top
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Functions]
  [top_compress]
    type = ParsedFunction
    expression = '-0.30 * t'   # 30% compression at t=1
  []
[]

# Finite-strain QuasiStatic action
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE            # finite (large) strain formulation
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_xy vonmises_stress'
  []
[]

[BCs]
  # Bottom: fully fixed
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []

  # Top: prescribed vertical compression, free in x
  [compress_top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = top_compress
  []

  # Left: symmetry (no horizontal displacement)
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
[]

[Materials]
  # Isotropic elasticity tensor
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 10.0        # MPa — soft rubber
    poissons_ratio = 0.45        # near-incompressible
  []
  # Finite strain elastic stress (St. Venant-Kirchhoff type)
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[Postprocessors]
  # Reaction force on top (average stress_yy * area)
  [avg_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
  [max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    value_type = max
  []
  [top_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 30
  l_max_its  = 50

  # Small steps for convergence of the nonlinear geometric problem
  dt       = 0.02
  end_time = 1.0
[]

[Outputs]
  exodus = true
  csv    = true
[]
