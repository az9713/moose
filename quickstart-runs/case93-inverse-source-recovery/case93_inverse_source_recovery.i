# ============================================================
# Case 93: Inverse Source Recovery — PDE-Constrained Estimation
# (Main optimizer file)
#
# Vogel, "Computational Methods for Inverse Problems"
#   (SIAM Frontiers in Applied Mathematics, 2002), Ch. 1–3
# Gunzburger, "Perspectives in Flow Control and Optimization"
#   (SIAM, 2003), Ch. 1–2
# Biegler et al. (eds.), "Large-Scale PDE-Constrained Optimization"
#   (Springer, 2003), Ch. 1
#
# Recover an unknown source distribution q(x,y) in a 2D screened
# Poisson equation from sparse interior measurements of the field u.
# Uses MOOSE's adjoint-based PDE-constrained optimization module.
#
# ============================================================
# MATHEMATICAL FORMULATION
# ============================================================
#
# FORWARD PROBLEM (screened Poisson equation):
#
#   −∇²u + α u = q(x,y)    on Ω = [0,1]²
#   u = 0                   on ∂Ω
#
# where α = 1.0 is the reaction (screening) coefficient and
# q(x,y) is the unknown source distribution to be recovered.
#
# This is the screened (or modified) Helmholtz / Yukawa equation.
# It governs:
#   - Steady-state diffusion with first-order decay
#   - The pressure correction equation in incompressible flow
#   - Green's function problems in potential theory
#   - Inverse design in diffusion-based manufacturing
#
# INVERSE PROBLEM:
#
# Given measurements u_obs at N sensor locations {(xᵢ, yᵢ)}:
#   Find q*(x,y) that minimises the data misfit functional:
#
#   J(q) = (1/2) Σᵢ [u(xᵢ, yᵢ; q) − u_obs,ᵢ]²
#
# subject to the forward PDE constraint (u solves the screened
# Poisson equation for the given q).
#
# ============================================================
# ADJOINT METHOD
# ============================================================
#
# The adjoint method computes the gradient dJ/dq analytically
# without finite differences, using a single additional PDE solve.
#
# Adjoint equation (same operator, different RHS):
#
#   −∇²λ + α λ = Σᵢ (u(xᵢ) − u_obs,ᵢ) δ(x − xᵢ)   on Ω
#   λ = 0                                               on ∂Ω
#
# The adjoint λ is the sensitivity of J to a unit perturbation
# of the source at each interior point. The Dirac delta sources
# on the RHS are the "measurement residuals" (misfits).
#
# Gradient of the objective functional:
#
#   dJ/dq = ∫ λ(x,y) × (∂f/∂q)(x,y) dΩ
#
# where f(x,y) = q × sin(πx)sin(πy) is the parameterised source
# shape and the inner product ∫ λ × sin(πx)sin(πy) dΩ gives the
# scalar gradient component for the amplitude parameter p₁.
#
# ============================================================
# PARAMETRISATION OF q
# ============================================================
#
# The unknown source is parameterised as:
#
#   q(x,y; p₁) = p₁ × sin(πx) × sin(πy)
#
# This single-parameter family is chosen so that the true source
# q_true(x,y) = 10 × sin(πx) × sin(πy) is exactly representable
# (p₁_true = 10). The parameterisation is a truncated modal expansion;
# the full problem would use many basis functions (pixels, splines, etc.)
# but a single mode suffices to demonstrate the algorithm.
#
# ============================================================
# ANALYTICAL SOLUTION
# ============================================================
#
# For q(x,y) = p₁ × sin(πx) × sin(πy) with u = 0 on ∂Ω:
#
# The exact solution is:
#   u(x,y) = p₁ / (2π² + α) × sin(πx) × sin(πy)
#
# Proof: substituting u = C × sin(πx) × sin(πy) into
#   −∇²u + u = p₁ sin(πx) sin(πy):
#   C(π² + π² + 1) = p₁  →  C = p₁/(2π² + 1)
#
# Numerically (α = 1, p₁ = 10):
#   C = 10 / (2π² + 1) = 10 / (19.739 + 1) = 10 / 20.739 = 0.48218...
#
# Note: the problem specification uses 2π²+1 ≈ 20.739, giving
# C ≈ 0.48218. The measurement values below use this exact value.
#
# Observation values at the 9 sensor locations (3×3 grid):
#
#   sin(π/4) = sin(3π/4) = √2/2 ≈ 0.70711
#   sin(π/2) = 1
#
#   u(0.25, 0.25) = 0.48218 × 0.70711 × 0.70711 = 0.48218 × 0.5 = 0.24109
#   u(0.50, 0.25) = 0.48218 × 1       × 0.70711 = 0.48218 × 0.70711 = 0.34099
#   u(0.75, 0.25) = 0.48218 × 0.70711 × 0.70711 = 0.24109
#   u(0.25, 0.50) = 0.48218 × 0.70711 × 1       = 0.34099
#   u(0.50, 0.50) = 0.48218 × 1       × 1       = 0.48218
#   u(0.75, 0.50) = 0.48218 × 0.70711 × 1       = 0.34099
#   u(0.25, 0.75) = 0.48218 × 0.70711 × 0.70711 = 0.24109
#   u(0.50, 0.75) = 0.48218 × 1       × 0.70711 = 0.34099
#   u(0.75, 0.75) = 0.48218 × 0.70711 × 0.70711 = 0.24109
#
# ============================================================
# OPTIMIZATION ALGORITHM
# ============================================================
#
# Algorithm: L-BFGS (Limited-memory BFGS, Broyden-Fletcher-
#   Goldfarb-Shanno quasi-Newton method via PETSc TAO)
# TAO solver: taoblmvm (L-BFGS with Moré–Thuente line search)
#
# For this linear PDE with a single parameter, the objective
# functional J(p₁) is a quadratic:
#   J(p₁) = (1/2) × K × (p₁ − p₁_true)²
#   dJ/dp₁ = K × (p₁ − p₁_true)
#
# where K = Σᵢ [C sin(πxᵢ) sin(πyᵢ)]². The L-BFGS method
# finds the exact minimum in ONE step for a quadratic objective.
#
# Convergence criterion: gradient tolerance |dJ/dp₁| < 1e-8.
# Starting guess: p₁⁰ = 1.0 (far from the true value of 10).
# Expected: ONE optimization iteration to reach p₁ ≈ 10.0.
#
# ============================================================
# THREE-FILE ARCHITECTURE
# ============================================================
#
# This main file orchestrates the optimization:
#   case93_inverse_source_recovery.i  — this file (TAO optimizer)
#   case93_forward.i                  — forward PDE solve
#   case93_adjoint.i                  — adjoint sensitivity solve
#
# Each TAO iteration:
#   1. FORWARD solve: given p₁, find u satisfying −∇²u + u = p₁ sin(πx)sin(πy)
#   2. Forward computes misfit: eᵢ = u(xᵢ) − u_obs,ᵢ and J = (1/2)Σeᵢ²
#   3. ADJOINT solve: given {eᵢ}, find λ satisfying −∇²λ + λ = Σᵢ eᵢ δ(x−xᵢ)
#   4. Adjoint computes gradient: dJ/dp₁ = ∫ λ sin(πx) sin(πy) dΩ
#   5. TAO updates p₁ using L-BFGS and the computed gradient.
# ============================================================

[Optimization]
  # The [Optimization] block activates the MOOSE Optimization module.
  # It registers the TAO-based executioner and the OptimizationReporter
  # system. No parameters are needed in this block; all configuration
  # is in [OptimizationReporter] and [Executioner].
[]

[Mesh]
  # The main optimizer mesh is NOT used for physics — it is required
  # syntactically by MOOSE but the actual PDE solves happen on the
  # sub-app meshes in case93_forward.i and case93_adjoint.i.
  # A minimal 1×1 mesh is sufficient here to satisfy the framework.
  [main_mesh]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 2
    ny   = 2
  []
[]

[OptimizationReporter]
  # GeneralOptimization reporter manages:
  #   - The parameter vector (p₁ = source amplitude)
  #   - The measurement data (9 sensor observations)
  #   - The objective value J and gradient dJ/dp₁
  #   - Communication between optimizer and sub-apps
  type = GeneralOptimization

  # The objective function value is collected from the forward app
  # (OptimizationData reporter computes (1/2)Σeᵢ²).
  objective_name = objective_value

  # Parameter names and sizes.
  # 'source_amplitude' is the label used throughout all Transfer blocks.
  # num_values = 1: only the amplitude p₁ is optimized (single-parameter
  # modal expansion of the source).
  parameter_names = 'source_amplitude'
  num_values      = '1'

  # Starting guess: p₁⁰ = 1.0 (true value is 10.0, so we start 10× too small).
  # L-BFGS should reach the optimum in 1 iteration for this linear problem.
  initial_condition = '1.0'
[]

[Reporters]
  # ------------------------------------------------------------------
  # OptimizationData 'main': stores measurement locations, observations,
  # and receives misfit values from the forward app.
  #
  # This reporter auto-parses measurement_points into:
  #   main/measurement_xcoord, main/measurement_ycoord, main/measurement_zcoord
  # and stores main/measurement_values, main/misfit_values.
  # These are transferred to the forward and adjoint sub-apps.
  # ------------------------------------------------------------------
  [main]
    type = OptimizationData
    measurement_points = '0.25 0.25 0
                           0.50 0.25 0
                           0.75 0.25 0
                           0.25 0.50 0
                           0.50 0.50 0
                           0.75 0.50 0
                           0.25 0.75 0
                           0.50 0.75 0
                           0.75 0.75 0'
    measurement_values = '0.24109 0.34099 0.24109
                           0.34099 0.48218 0.34099
                           0.24109 0.34099 0.24109'
  []
[]

[Executioner]
  # Optimize executioner: drives the TAO solver to minimise J(p₁).
  type = Optimize

  # taoblmvm: Limited-memory BFGS with Moré–Thuente line search.
  # BLMVM is the PETSc TAO implementation of L-BFGS (Liu & Nocedal 1989).
  # It approximates the inverse Hessian from the last m=5 gradient vectors
  # (memory parameter m is set via -tao_lmm_vectors if needed).
  # For a strictly convex quadratic J, BLMVM finds the exact minimum
  # in ONE step (since the true Hessian is constant and L-BFGS converges
  # exactly after m steps for an m-dimensional quadratic).
  tao_solver = taoblmvm

  # TAO solver tolerances:
  #   -tao_gatol 1e-8: gradient absolute tolerance. Stop when |∇J| < 1e-8.
  #   For the quadratic J = (1/2)K(p₁−10)² with K~0.1, |∇J|=K|p₁−10|.
  #   This requires |p₁−10| < 1e-8/K ≈ 1e-7, i.e., recovery to 8 digits.
  #   -tao_ls_type unit: use a unit (trivial) step for line search.
  #   For a quadratic objective, the exact Newton step is optimal so
  #   no line search is needed. Setting unit avoids wasted function evals.
  petsc_options_iname = '-tao_gatol -tao_ls_type'
  petsc_options_value = '1e-8       unit'

  # Print TAO convergence information to console at each iteration.
  # Shows: iteration number, objective J, gradient norm |∇J|, step size.
  verbose = true
[]

[MultiApps]
  # Forward sub-app: solves the screened Poisson PDE for given p₁.
  # FullSolveMultiApp runs the sub-app to steady-state convergence.
  # execute_on = FORWARD: run during the forward pass of each TAO iteration.
  [forward]
    type        = FullSolveMultiApp
    input_files = case93_forward.i
    execute_on  = FORWARD
  []

  # Adjoint sub-app: solves the adjoint (sensitivity) equation.
  # execute_on = ADJOINT: run during the gradient computation pass.
  [adjoint]
    type        = FullSolveMultiApp
    input_files = case93_adjoint.i
    execute_on  = ADJOINT
  []
[]

[Transfers]
  # ==================================================================
  # MAIN → FORWARD transfers: send current parameter to the forward app
  # ==================================================================

  # Transfer the current optimizer parameter p₁ to the forward app.
  # The forward app receives it in the ConstantReporter 'src_values/values'
  # and uses it in the ParsedOptimizationFunction that defines the source.
  [to_forward_params]
    type           = MultiAppReporterTransfer
    to_multi_app   = forward
    from_reporters = 'OptimizationReporter/source_amplitude'
    to_reporters   = 'src_values/values'
  []

  # Transfer measurement locations and observed values to the forward app.
  # The forward app needs these to compute the misfit eᵢ = u(xᵢ) − u_obs,ᵢ
  # and the objective value J = (1/2)Σeᵢ².
  [to_forward_data]
    type           = MultiAppReporterTransfer
    to_multi_app   = forward
    from_reporters = 'main/measurement_xcoord
                      main/measurement_ycoord
                      main/measurement_zcoord
                      main/measurement_time
                      main/measurement_values'
    to_reporters   = 'measure_data/measurement_xcoord
                      measure_data/measurement_ycoord
                      measure_data/measurement_zcoord
                      measure_data/measurement_time
                      measure_data/measurement_values'
  []

  # ==================================================================
  # FORWARD → MAIN transfers: collect objective and misfit
  # ==================================================================

  # Collect objective value J and misfit vector {eᵢ} from forward app.
  # The OptimizationReporter needs J to pass to TAO, and the misfit
  # {eᵢ} is needed by the adjoint app (it forms the adjoint RHS).
  [from_forward]
    type           = MultiAppReporterTransfer
    from_multi_app = forward
    from_reporters = 'measure_data/misfit_values measure_data/objective_value'
    to_reporters   = 'main/misfit_values OptimizationReporter/objective_value'
  []

  # ==================================================================
  # MAIN → ADJOINT transfers: send sensor locations and misfit
  # ==================================================================

  # Transfer sensor coordinates and misfit {eᵢ} to the adjoint app.
  # The adjoint RHS is: Σᵢ eᵢ δ(x − xᵢ), implemented as ReporterPointSource.
  # The adjoint needs the misfit (not the observations) as the source weights.
  [to_adjoint]
    type           = MultiAppReporterTransfer
    to_multi_app   = adjoint
    from_reporters = 'main/measurement_xcoord
                      main/measurement_ycoord
                      main/measurement_zcoord
                      main/misfit_values'
    to_reporters   = 'misfit_data/measurement_xcoord
                      misfit_data/measurement_ycoord
                      misfit_data/measurement_zcoord
                      misfit_data/misfit_values'
  []

  # Transfer current parameter values to the adjoint app for the
  # ParsedOptimizationFunction used in the gradient inner product.
  [to_adjoint_params]
    type           = MultiAppReporterTransfer
    to_multi_app   = adjoint
    from_reporters = 'OptimizationReporter/source_amplitude'
    to_reporters   = 'params/values'
  []

  # ==================================================================
  # ADJOINT → MAIN transfers: collect gradient
  # ==================================================================

  # Collect the gradient dJ/dp₁ = ∫ λ sin(πx) sin(πy) dΩ from the adjoint app.
  # This inner product is computed by ElementOptimizationSourceFunctionInnerProduct.
  # TAO uses this gradient to update p₁ in the L-BFGS step.
  [from_adjoint]
    type           = MultiAppReporterTransfer
    from_multi_app = adjoint
    from_reporters = 'gradient/inner_product'
    to_reporters   = 'OptimizationReporter/grad_source_amplitude'
  []
[]

[Outputs]
  # CSV output for the main app:
  # Records the optimization history: iteration, p₁, J, |∇J|.
  # This allows post-processing of the convergence curve.
  csv = true
[]
