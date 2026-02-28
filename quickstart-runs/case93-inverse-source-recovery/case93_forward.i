# ============================================================
# Case 93 — Forward App: Screened Poisson Equation
#
# Solves the forward PDE for the current source amplitude p₁:
#
#   −∇²u + α u = q(x,y; p₁)    on Ω = [0,1]²
#   u = 0                        on ∂Ω
#
# where α = 1.0 and q(x,y; p₁) = p₁ × sin(πx) × sin(πy).
#
# After solving, the OptimizationData reporter:
#   (1) Evaluates u at the 9 sensor locations {(xᵢ, yᵢ)}
#   (2) Computes the misfit:  eᵢ = u(xᵢ) − u_obs,ᵢ
#   (3) Computes the objective: J = (1/2) Σᵢ eᵢ²
#   (4) Returns {eᵢ} and J to the main optimizer app
#
# Called by the main app (case93_inverse_source_recovery.i) at
# execute_on = FORWARD in every L-BFGS iteration.
#
# ============================================================
# KERNEL SIGN CONVENTIONS
# ============================================================
#
# MOOSE kernel residuals (all contributions ADDED to residual R,
# solved as R = 0):
#
#   Diffusion:     R += +∫ ∇u · ∇v dΩ    (strong: −∇²u)
#   CoefReaction:  R += +α ∫ u · v dΩ    (strong: +α u)
#   BodyForce:     R += −∫ f · v dΩ      (strong: −f)
#
# Combined residual R = 0:
#   ∫(∇u·∇v + αuv − fv) dΩ = 0
#   ⟺  −∇²u + αu − f = 0
#   ⟺  −∇²u + αu = f   ✓
#
# where f = q = p₁ × sin(πx) × sin(πy) (the parameterised source).
#
# Note on CoefReaction sign:
#   CoefReaction adds +coefficient × u × test to the residual.
#   This gives the strong-form term +α u (POSITIVE), which is
#   the screening term on the LHS of the screened Poisson equation.
#   This is DIFFERENT from ADMatReaction / MatReaction, which add
#   NEGATIVE rate × u × test (strong form: −rate × u).
# ============================================================

[Mesh]
  # 2D unit square [0,1] × [0,1] with 30×30 bilinear Q1 elements.
  # Element size: h = 1/30 ≈ 0.0333.
  #
  # Resolution check for the screened Poisson problem:
  #   Characteristic length:   L_c = 1/√α = 1/√1 = 1 (full domain)
  #   Source spatial scale:    λ_source = 1 (sin(πx) has period 2)
  #   Elements per period:     30 × 2 = 60 elements per sin period
  #
  # 30×30 is more than adequate for smooth sin-mode solutions.
  # The FEM solution converges as O(h²) for linear elements on smooth
  # functions; at h = 1/30, the relative error in u is O(1/900) ≈ 0.1%.
  # This error level ensures that the optimization convergence is
  # limited by the gradient tolerance (1e-8), not by discretisation.
  [mesh]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 30
  []
[]

[Variables]
  # u: solution of the screened Poisson equation.
  # First-order (linear) Lagrange elements on the 30×30 mesh.
  # C⁰ continuity is enforced at all interior nodes and boundaries.
  # Dirichlet BCs (u = 0) on all four walls set the solution to zero
  # at the boundary nodes, consistent with the Green's function structure.
  [u]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # ------------------------------------------------------------------
  # Diffusion term: contributes −∇²u to the strong-form equation.
  # Standard Laplacian diffusion with coefficient = 1 (no material prop).
  # Weak form: +∫ ∇u · ∇v dΩ (integration by parts of ∇²u with zero BCs).
  # ------------------------------------------------------------------
  [diffusion]
    type     = Diffusion
    variable = u
  []

  # ------------------------------------------------------------------
  # Reaction (screening) term: contributes +α u to the strong form.
  # CoefReaction adds +coefficient × ∫ u v dΩ to the residual (POSITIVE).
  # With coefficient = 1.0 (α = 1): adds the term +u to the equation.
  # Combined with Diffusion: −∇²u + u = f (screened Poisson) ✓
  #
  # Note: Do NOT use ADMatReaction or MatReaction here — those add
  # NEGATIVE rate × u (strong: −rate×u), which would give −∇²u − u = f,
  # i.e., the Helmholtz equation with NEGATIVE k², not screened Poisson.
  # ------------------------------------------------------------------
  [screening]
    type        = CoefReaction
    variable    = u
    coefficient = 1.0
  []

  # ------------------------------------------------------------------
  # Source term: the parameterised source q(x,y; p₁) = p₁ sin(πx) sin(πy).
  # BodyForce adds −∫ f v dΩ to the residual (NEGATIVE sign).
  # With residual R = ∫∇u·∇v + uv dΩ − ∫fv dΩ = 0:
  #   strong form: −∇²u + u = f   (where f = p₁ sin(πx) sin(πy)) ✓
  #
  # The function 'source_fn' is a ParsedOptimizationFunction that reads
  # the current parameter p₁ from the reporter 'src_values/values',
  # which is updated by the main optimizer before each forward solve.
  # ------------------------------------------------------------------
  [source]
    type     = BodyForce
    variable = u
    function = source_fn
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # ParsedOptimizationFunction for the source q(x,y; p₁)
  #
  # This function links the optimizer's parameter vector to the PDE source.
  # At each Newton iteration within the forward solve, MOOSE evaluates:
  #   source_fn(x,y) = p₁ × sin(π × x) × sin(π × y)
  # using the current p₁ from 'src_values/values'.
  #
  # The expression 'p0 * sin(pi * x) * sin(pi * y)' uses the symbol p0,
  # which is bound to the first component of the reporter vector 'src_values/values'.
  # This reporter is updated by the main app Transfer [to_forward_params]
  # before each FORWARD execute_on call.
  #
  # For the gradient computation in the adjoint app, we need the same
  # functional form (the "basis function" for the parameter). The adjoint
  # app uses the analogous ParsedFunction (not ParsedOptimizationFunction)
  # to evaluate the basis shape dq/dp₁ = sin(πx)sin(πy).
  # ------------------------------------------------------------------
  [source_fn]
    type               = ParsedOptimizationFunction
    expression         = 'p0 * sin(pi * x) * sin(pi * y)'
    param_symbol_names = 'p0'
    param_vector_name  = 'src_values/values'
  []
[]

[BCs]
  # Dirichlet BCs: u = 0 on all four boundaries.
  # This is consistent with the Green's function for the screened Poisson
  # equation on a finite domain (homogeneous Dirichlet = "grounded" walls).
  # The zero BCs ensure that the analytical solution C sin(πx) sin(πy)
  # automatically satisfies the boundary conditions (sin vanishes at 0,1).
  [all_walls]
    type     = DirichletBC
    variable = u
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Reporters]
  # ------------------------------------------------------------------
  # ConstantReporter 'src_values': receives the optimizer parameter p₁.
  #
  # The main app transfers OptimizationReporter/source_amplitude to
  # src_values/values before each forward solve. The ParsedOptimizationFunction
  # 'source_fn' reads from this reporter at each quadrature point evaluation.
  #
  # Initial value '1.0' matches the initial_condition of the OptimizationReporter,
  # ensuring the first forward solve uses the correct starting guess.
  # ------------------------------------------------------------------
  [src_values]
    type                  = ConstantReporter
    real_vector_names     = 'values'
    real_vector_values    = '1.0'
  []

  # ------------------------------------------------------------------
  # OptimizationData 'measure_data': evaluation, misfit, and objective.
  #
  # This reporter:
  #   (1) Stores measurement locations (xᵢ, yᵢ, zᵢ) and observations u_obs,ᵢ
  #       received from the main app via [to_forward_data] transfer.
  #   (2) After the forward solve, evaluates u(xᵢ, yᵢ) by interpolation.
  #   (3) Computes misfits: eᵢ = u(xᵢ) − u_obs,ᵢ
  #   (4) Computes the objective: objective_value = (1/2) Σᵢ eᵢ²
  #   (5) Exports misfit_values and objective_value for transfer back to main.
  #
  # The 'variable' parameter tells the reporter which FEM variable to
  # evaluate at the measurement points (u in this case).
  # ------------------------------------------------------------------
  [measure_data]
    type             = OptimizationData
    objective_name   = objective_value
    variable         = u
  []
[]

[Executioner]
  # Steady-state Newton solve for the linear screened Poisson system.
  # The system is linear (u appears at most linearly in all kernels),
  # so Newton converges in exactly ONE iteration (the linear residual
  # equals the full nonlinear residual for a linear PDE).
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 30×30 mesh with 1 variable gives
  # 961 DOFs (31×31 nodes minus Dirichlet BCs ≈ 841 free DOFs).
  # This is tiny for direct LU; solve time is negligible.
  # Direct LU is preferred over iterative solvers because:
  #   (1) The system is small enough for exact factorisation.
  #   (2) The screened Poisson operator (−∇² + I) is SPD → CG would
  #       also work, but LU avoids any convergence concerns.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-14
[]

[Outputs]
  # Output the forward solution field for diagnostics:
  #   - Exodus: visualize u(x,y) to verify the sin-mode shape
  #   - CSV: misfit and objective values at each optimizer iteration
  # Both outputs are suppressed in the sub-app console to reduce noise.
  exodus = true
  csv    = true
[]
