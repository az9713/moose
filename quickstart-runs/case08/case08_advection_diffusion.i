# ============================================================
# Case 8: Transient Advection-Diffusion
# dc/dt + div(v*c) - div(D*grad c) = 0
# v = (0.5, 0, 0),  D = 0.01
# ============================================================

vx = 0.5    # x-velocity component  (HIT top-level scalar variable)
D  = 0.01   # diffusion coefficient

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 60
  ny   = 30
  xmin = 0
  xmax = 3     # longer domain so the blob has room to travel
  ymin = 0
  ymax = 1
[]

[Variables]
  [c]
  []
[]

[ICs]
  # Start with a Gaussian blob centered at (0.3, 0.5).
  [blob]
    type     = FunctionIC
    variable = c
    function = blob_fn
  []
[]

[Functions]
  [blob_fn]
    type       = ParsedFunction
    expression = 'exp(-((x-0.3)^2 + (y-0.5)^2) / 0.01)'
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = c
  []

  # ConservativeAdvection implements  -int( grad(phi_i) . v*c ) dV
  # which is the conservative (divergence) form of advection.
  # The velocity is provided as a material property (vector type).
  [advection]
    type              = ConservativeAdvection
    variable          = c
    velocity_material = velocity_vec   # references a vector material property
  []

  [diffusion]
    type        = MatDiffusion
    variable    = c
    diffusivity = D_coeff
  []
[]

[BCs]
  # On the right boundary let the concentration flow out naturally
  # using a zero-flux Neumann BC (automatic, no entry needed).
  # On the left, let the advective flux carry concentration out
  # using a DirichletBC that fixes c=0 (outflow if v is into domain):
  [left_wall]
    type     = DirichletBC
    variable = c
    boundary = left
    value    = 0
  []
[]

[Materials]
  # GenericConstantVectorMaterial declares a RealVectorValue material property.
  # prop_values lists all components: (vx, vy, vz) for one property.
  [velocity_mat]
    type        = GenericConstantVectorMaterial
    prop_names  = 'velocity_vec'
    prop_values = '${vx} 0 0'   # v = (0.5, 0, 0)
  []

  [diffusion_mat]
    type        = GenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D}'
  []
[]

[Postprocessors]
  # Track total concentration (should be approximately conserved
  # until concentration reaches the outflow boundary).
  [total_c]
    type     = ElementIntegralVariablePostprocessor
    variable = c
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'ilu'

  dt         = 0.02
  end_time   = 2.0

  nl_rel_tol = 1e-6
[]

[Outputs]
  exodus = true
  csv    = true
[]
