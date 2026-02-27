# ============================================================
# Case 47 Forward App: Steady Heat Conduction with Unknown Source
# Solves -k*laplacian(T) = q on [0,2]x[0,2]
# T = 200 (bottom), T = 100 (top), insulated sides
# q is received from the optimizer via ParsedOptimizationFunction
# ============================================================

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmax = 2
    ymax = 2
  []
[]

[Variables]
  [temperature]
  []
[]

[Kernels]
  [heat_conduction]
    type = MatDiffusion
    variable = temperature
    diffusivity = thermal_conductivity
  []
  [heat_source]
    type = BodyForce
    function = volumetric_heat_func
    variable = temperature
  []
[]

[BCs]
  [left]
    type = NeumannBC
    variable = temperature
    boundary = left
    value = 0
  []
  [right]
    type = NeumannBC
    variable = temperature
    boundary = right
    value = 0
  []
  [bottom]
    type = DirichletBC
    variable = temperature
    boundary = bottom
    value = 200
  []
  [top]
    type = DirichletBC
    variable = temperature
    boundary = top
    value = 100
  []
[]

[Functions]
  [volumetric_heat_func]
    type = ParsedOptimizationFunction
    expression = q
    param_symbol_names = 'q'
    param_vector_name = 'params/q'
  []
[]

[Materials]
  [steel]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 5
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Reporters]
  [measure_data]
    type = OptimizationData
    objective_name = objective_value
    variable = temperature
  []
  [params]
    type = ConstantReporter
    real_vector_names = 'q'
    real_vector_values = '0'
  []
[]

[Outputs]
  console = false
[]
