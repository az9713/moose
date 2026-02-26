# ============================================================
# Case 18: Cahn-Hilliard Spinodal Decomposition
# Split form of the Cahn-Hilliard equation for phase separation
# in a binary mixture on a periodic domain.
#
# Equations (split form):
#   dc/dt  = div( M * grad(w) )          [concentration evolution]
#   w      = dF/dc - kappa * laplacian(c) [chemical potential]
#
# Free energy:  F(c) = c^2 * (1-c)^2    (double-well)
# Mobility:     M = 1.0
# Gradient energy coefficient: kappa = 1.0
#
# Initial condition: c = 0.5 + small random noise
# Boundary conditions: periodic in x and y
# Domain: [0,25] x [0,25],  40x40 mesh
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 40
  ny   = 40
  xmin = 0
  xmax = 25
  ymin = 0
  ymax = 25
[]

[Variables]
  # Concentration / order parameter.
  # c = 0 and c = 1 are the two equilibrium phases.
  [c]
    order  = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = RandomIC
      min  = 0.44
      max  = 0.56
    []
  []

  # Chemical potential (split variable).
  # Splitting the 4th-order CH equation into two coupled 2nd-order
  # equations avoids the need for C1-continuous (Hermite) elements.
  [w]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # -------------------------------------------------------
  # Equation 1 (solved on variable w):
  #   dc/dt - div( M * grad(w) ) = 0
  # -------------------------------------------------------

  # CoupledTimeDerivative: adds  (dc/dt, phi_w) on variable w.
  # This couples the time rate of c into the w equation.
  [c_dot]
    type     = CoupledTimeDerivative
    variable = w
    v        = c
  []

  # SplitCHWRes: adds  ( M * grad(w), grad(phi_w) ) on variable w.
  # This is the mobility-weighted diffusion of the chemical potential.
  [w_res]
    type     = SplitCHWRes
    variable = w
    mob_name = M
  []

  # -------------------------------------------------------
  # Equation 2 (solved on variable c):
  #   w - dF/dc + kappa * laplacian(c) = 0
  # -------------------------------------------------------

  # SplitCHParsed: adds the chemical potential equation.
  # It reads dF/dc from the DerivativeParsedMaterial and the
  # gradient energy term using kappa_c.
  [c_res]
    type       = SplitCHParsed
    variable   = c
    f_name     = F
    kappa_name = kappa_c
    w          = w
  []
[]

[BCs]
  # Periodic boundary conditions eliminate boundary effects and let
  # the simulation model bulk behaviour of an infinite medium.
  [Periodic]
    [all]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  # Double-well free energy: F(c) = c^2 * (1-c)^2
  # Minima at c=0 and c=1 drive phase separation.
  # DerivativeParsedMaterial automatically computes dF/dc and d2F/dc2
  # which are required by SplitCHParsed and the Jacobian assembly.
  [free_energy]
    type             = DerivativeParsedMaterial
    property_name    = F
    coupled_variables = 'c'
    expression       = 'c^2*(1-c)^2'
    derivative_order = 2
    disable_fpoptimizer = true
    enable_jit          = false   # JIT requires mpicxx; disable for Docker portability
  []

  # Constant material properties:
  #   kappa_c = 1.0  (gradient energy penalty; controls interface width)
  #   M       = 1.0  (atomic mobility; controls kinetics)
  [const]
    type        = GenericConstantMaterial
    prop_names  = 'kappa_c M'
    prop_values = '1.0     1.0'
  []
[]

[Postprocessors]
  # Volume-averaged concentration — should remain ~0.5 throughout
  # (mass conservation check: total solute is conserved by the CH equation).
  [avg_c]
    type       = ElementAverageValue
    variable   = c
    execute_on = 'initial timestep_end'
  []

  # Integral of the bulk free energy F(c) over the domain.
  # As spinodal decomposition proceeds the bulk free energy decreases
  # because more material sits near the equilibrium values c=0 and c=1.
  [bulk_energy]
    type          = ElementIntegralMaterialProperty
    mat_prop      = F
    execute_on    = 'initial timestep_end'
  []

  # Minimum and maximum concentration: tracks how far the phases have
  # separated. Initially both stay near 0.5; they should approach 0 and 1
  # as domains coarsen.
  [min_c]
    type       = ElementExtremeValue
    variable   = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [max_c]
    type       = ElementExtremeValue
    variable   = c
    value_type = max
    execute_on = 'initial timestep_end'
  []

  # Current time step size (diagnostic).
  [dt_size]
    type = TimestepSize
  []
[]

[Executioner]
  type       = Transient
  solve_type = 'NEWTON'
  scheme     = bdf2

  # AMG preconditioner well-suited for the elliptic sub-problems
  # arising in the split CH system.
  # LU direct solver is the most robust choice for the split CH system
  # on a small 2-D mesh. AMG can struggle with the saddle-point structure.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  l_max_its  = 30
  nl_max_its = 30
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11

  start_time = 0.0
  end_time   = 100.0

  # Adaptive time stepping: grow the step when Newton converges easily,
  # cut it when convergence is slow. Early spinodal kinetics are fast;
  # late-stage coarsening is slow — adaptive dt captures both efficiently.
  [TimeStepper]
    type           = IterationAdaptiveDT
    dt             = 0.1
    growth_factor  = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv    = true

  # Keep the two most recent checkpoints so the run can be restarted.
  [checkpoint]
    type      = Checkpoint
    num_files = 2
  []
[]
