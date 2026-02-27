# ============================================================
# Case 70: Contact Mechanics — Mortar Frictionless Contact
# Two elastic blocks. Left block is pushed into the right
# block by a prescribed displacement. Mortar-based
# frictionless contact prevents interpenetration.
#
# Left block: [-1, 0] x [-0.5, 0.5]
# Right block: [0, 1] x [-0.6, 0.6]
# Initial gap: 0.01 (left block shifted left by gap)
#
# BCs:
#   Left boundary (lb_left): disp_x = 0.1*t, disp_y = 0
#   Right boundary (rb_right): disp_x = 0, disp_y = 0
#
# Timeline:
#   t < 0.1: gap closing (no contact force)
#   t = 0.1: contact initiates
#   t > 0.1: contact pressure builds
#   t = 0.5: total displacement = 0.05, compression = 0.04
#
# Materials: E = 1e6 Pa, nu = 0.3 (both blocks)
# Mesh: 10x10 left, 10x12 right (MeshCollectionGenerator)
# ============================================================

gap = 0.01

[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1.0
    xmax = 0.0
    ymin = -0.5
    ymax = 0.5
    nx = 10
    ny = 10
    elem_type = QUAD4
    boundary_name_prefix = lb
  []
  [left_block_id]
    type = SubdomainIDGenerator
    input = left_block
    subdomain_id = 1
  []
  [right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.0
    xmax = 1.0
    ymin = -0.6
    ymax = 0.6
    nx = 10
    ny = 12
    elem_type = QUAD4
    boundary_name_prefix = rb
    boundary_id_offset = 10
  []
  [right_block_id]
    type = SubdomainIDGenerator
    input = right_block
    subdomain_id = 2
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'left_block_id right_block_id'
  []
  [block_rename]
    type = RenameBlockGenerator
    input = combined
    old_block = '1 2'
    new_block = 'left_block right_block'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    block = 'left_block right_block'
  []
[]

[AuxVariables]
  [vonmises]
    order = CONSTANT
    family = MONOMIAL
    block = 'left_block right_block'
  []
[]

[AuxKernels]
  [vonmises_kernel]
    type = RankTwoScalarAux
    variable = vonmises
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = 'timestep_end'
    block = 'left_block right_block'
  []
[]

[Functions]
  [horizontal_movement]
    type = ParsedFunction
    expression = '0.1*t'
  []
[]

[BCs]
  # Push left block rightward
  [push_x]
    type = FunctionDirichletBC
    preset = true
    variable = disp_x
    boundary = lb_left
    function = horizontal_movement
  []
  # Fix left block vertically
  [fix_y_left]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = lb_left
    value = 0.0
  []
  # Fix right block completely
  [fix_x_right]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = rb_right
    value = 0.0
  []
  [fix_y_right]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = rb_right
    value = 0.0
  []
[]

[Contact]
  [leftright]
    secondary = lb_right
    primary = rb_left
    model = frictionless
    formulation = mortar
  []
[]

[ICs]
  # Shift left block to create initial gap
  [disp_x_left]
    type = ConstantIC
    block = left_block
    variable = disp_x
    value = -${gap}
  []
[]

[Materials]
  [elasticity_left]
    type = ComputeIsotropicElasticityTensor
    block = left_block
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
  []
  [stress_left]
    type = ComputeFiniteStrainElasticStress
    block = left_block
  []
  [elasticity_right]
    type = ComputeIsotropicElasticityTensor
    block = right_block
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
  []
  [stress_right]
    type = ComputeFiniteStrainElasticStress
    block = right_block
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [max_disp_x]
    type = NodalExtremeValue
    variable = disp_x
    block = 'left_block right_block'
    execute_on = 'initial timestep_end'
  []
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises
    block = 'left_block right_block'
    execute_on = 'initial timestep_end'
  []
  [avg_vonmises]
    type = ElementAverageValue
    variable = vonmises
    block = 'left_block right_block'
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       1e-5          NONZERO               1e-15'
  dt = 0.05
  end_time = 0.5
  l_tol = 1e-4
  l_max_its = 100
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 50
[]

[Outputs]
  exodus = true
  csv = true
[]
