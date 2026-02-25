# ============================================================
# Case 1: Steady-State 1-D Diffusion
# Solves -d^2u/dx^2 = 0, u(0)=0, u(1)=1
# Exact solution: u(x) = x
# ============================================================

[Mesh]
  # GeneratedMesh builds a structured mesh without an external file.
  type = GeneratedMesh
  dim  = 1        # 1-D problem: a line segment
  nx   = 20       # 20 uniform elements along x
  xmin = 0        # left endpoint
  xmax = 1        # right endpoint
[]

[Variables]
  # u is the primary unknown.  The default FE family is LAGRANGE,
  # order FIRST (linear elements), which is appropriate here.
  [u]
  []
[]

[Kernels]
  # The Diffusion kernel contributes the weak-form integral
  #   int( grad(phi_i) . grad(u) ) dV
  # to the residual.  It implements the -div(grad(u)) operator.
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]

[BCs]
  # DirichletBC pins the value of u at a boundary.
  # Named boundaries on a GeneratedMesh in 1-D are 'left' and 'right'.
  [pin_left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [pin_right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
[]

[Executioner]
  # Steady means we solve F(u)=0 once, with no time loop.
  type = Steady

  # PJFNK: Preconditioned Jacobian-Free Newton-Krylov.
  # This is the workhorse nonlinear solver in MOOSE.
  solve_type = 'PJFNK'

  # BoomerAMG is an algebraic multigrid preconditioner from PETSc/Hypre.
  # It works well for elliptic problems.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  # Write an Exodus II file (case01_diffusion_1d_out.e).
  # Exodus is the preferred format for visualization in ParaView or VisIt.
  exodus = true

  # Also write a plain-text CSV summary (useful for quick checks).
  csv = true
[]
