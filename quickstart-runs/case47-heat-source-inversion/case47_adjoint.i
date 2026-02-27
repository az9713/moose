# ============================================================
# Case 47 Adjoint App: Sensitivity Equation for Heat Source
# Solves the adjoint of the forward heat equation:
#   -k*laplacian(lambda) = 0     in the domain
# with DiracKernel point sources at measurement locations
# weighted by (T_sim - T_meas), i.e., the misfit.
#
# The gradient of the objective is then:
#   dJ/dq = integral(lambda * d(forward_PDE)/dq) dV
# For a constant body force q, this simplifies to:
#   dJ/dq = integral(lambda) dV
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
  [adjoint_T]
  []
[]

[Kernels]
  [heat_conduction]
    type = MatDiffusion
    variable = adjoint_T
    diffusivity = thermal_conductivity
  []
[]

[DiracKernels]
  [pt]
    type = ReporterPointSource
    variable = adjoint_T
    x_coord_name = misfit/measurement_xcoord
    y_coord_name = misfit/measurement_ycoord
    z_coord_name = misfit/measurement_zcoord
    value_name = misfit/misfit_values
  []
[]

[Reporters]
  [misfit]
    type = OptimizationData
  []
  [params]
    type = ConstantReporter
    real_vector_names = 'q'
    real_vector_values = '0'
  []
[]

[BCs]
  [left]
    type = NeumannBC
    variable = adjoint_T
    boundary = left
    value = 0
  []
  [right]
    type = NeumannBC
    variable = adjoint_T
    boundary = right
    value = 0
  []
  [bottom]
    type = DirichletBC
    variable = adjoint_T
    boundary = bottom
    value = 0
  []
  [top]
    type = DirichletBC
    variable = adjoint_T
    boundary = top
    value = 0
  []
[]

[Materials]
  [steel]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 5
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

[VectorPostprocessors]
  [gradient_vpp]
    type = ElementOptimizationSourceFunctionInnerProduct
    variable = adjoint_T
    function = volumetric_heat_func
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  console = false
[]
