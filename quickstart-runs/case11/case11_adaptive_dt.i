# ============================================================
# Case 11: Adaptive Time Stepping with IterationAdaptiveDT
# Transient heat equation with a step-change source at t=0
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []
  [diffusion]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []
  [source]
    type     = BodyForce
    variable = T
    # A step-function heat source that switches on at t=0.
    # We use a large constant value to create rapid initial transient.
    value    = 100.0
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  [props]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '1.0'
  []
[]

[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []
  # TimestepSize reports the current dt — useful for verifying adaptation.
  [dt_pp]
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # IterationAdaptiveDT adjusts dt based on Newton iteration count.
  [TimeStepper]
    type = IterationAdaptiveDT

    # Starting timestep.
    dt = 0.001

    # Target number of Newton iterations per timestep.
    # If actual < this, grow dt. If actual > this (+ window), shrink dt.
    optimal_iterations = 5

    # Multiply dt by this factor when fewer than optimal_iterations.
    growth_factor = 2.0

    # Multiply dt by this factor when too many iterations.
    cutback_factor = 0.5

    # The iteration window around optimal_iterations that is acceptable
    # without triggering growth or cutback.
    iteration_window = 2
  []

  start_time = 0.0
  end_time   = 1.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
