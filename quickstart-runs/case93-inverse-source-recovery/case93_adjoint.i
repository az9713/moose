# ============================================================
# Case 93 — Adjoint App: Sensitivity Equation
#
# Solves the adjoint of the screened Poisson equation to compute
# the gradient dJ/dp₁ of the objective functional with respect
# to the source amplitude parameter p₁.
#
# Adjoint PDE (same operator as forward, adjoint RHS):
#
#   −∇²λ + α λ = Σᵢ eᵢ δ(x − xᵢ)    on Ω = [0,1]²
#   λ = 0                              on ∂Ω
#
# where eᵢ = u(xᵢ) − u_obs,ᵢ are the misfits from the forward solve,
# α = 1.0, and δ(x − xᵢ) are Dirac point sources at sensor locations.
#
# After solving, the gradient is computed as the inner product:
#
#   dJ/dp₁ = ∫ λ(x,y) × (dq/dp₁)(x,y) dΩ
#           = ∫ λ(x,y) × sin(πx) sin(πy) dΩ
#
# This inner product is computed by ElementOptimizationSourceFunctionInnerProduct.
#
# ============================================================
# WHY THE ADJOINT EQUATION HAS THE SAME OPERATOR
# ============================================================
#
# For a self-adjoint operator L = −∇² + αI (with Dirichlet BCs),
# the adjoint operator L* = L (self-adjoint). Therefore the adjoint
# equation is:
#
#   L* λ = eᵢ-weighted Dirac sum
#   ⟺  L λ = same RHS   (same PDE, different source)
#
# This means the forward and adjoint apps share the SAME kernel
# structure (Diffusion + CoefReaction) and the SAME boundary
# conditions (Dirichlet zero). The only difference is the source:
#   Forward: BodyForce with smooth ParsedOptimizationFunction
#   Adjoint: ReporterPointSource (Dirac deltas at sensor locations)
#
# ============================================================
# GRADIENT FORMULA DERIVATION
# ============================================================
#
# Lagrangian: L(u, λ, p₁) = J(u) + ∫ λ[−∇²u + αu − q(p₁)] dΩ
#
# Stationarity with respect to u (adjoint equation):
#   ∂L/∂u = 0:  −∇²λ + αλ = Σᵢ eᵢ δ(x − xᵢ)
#
# Stationarity with respect to λ (forward PDE):
#   ∂L/∂λ = 0:  −∇²u + αu = q(p₁)
#
# Gradient:
#   dJ/dp₁ = ∂L/∂p₁ = −∫ λ × ∂q/∂p₁ dΩ
#
# But note: BodyForce adds −∫ q v to the residual, meaning the
# Lagrangian sign for the source term is +∫ λ × ∂q/∂p₁ dΩ.
# The ElementOptimizationSourceFunctionInnerProduct computes
# ∫ λ × (basis function) dΩ, which is dJ/dp₁ (with consistent sign).
#
# For q = p₁ × sin(πx) sin(πy):
#   ∂q/∂p₁ = sin(πx) sin(πy)   (basis function 'grad_basis_fn')
#   dJ/dp₁ = ∫ λ × sin(πx) sin(πy) dΩ
#
# ============================================================
# DIRAC DELTA SOURCES (ReporterPointSource)
# ============================================================
#
# The misfit-weighted Dirac sources on the adjoint RHS are:
#
#   f_adj(x) = Σᵢ eᵢ δ(x − xᵢ)
#
# MOOSE implements this as a DiracKernel (ReporterPointSource) that
# adds a point load at each measurement location weighted by eᵢ.
# The Dirac DiracKernel adds to the residual at the element containing xᵢ,
# distributing the point source to the nearby nodes using element shape
# functions. This is the finite-element regularisation of the delta.
#
# The misfit values {eᵢ} and sensor coordinates {(xᵢ,yᵢ)} are
# received from the main app via [to_adjoint] Transfer and stored
# in the ConstantReporter 'misfit_data'.
# ============================================================

[Mesh]
  # Identical mesh to the forward app: 30×30 unit square.
  # The forward and adjoint problems must be solved on the SAME mesh
  # because the gradient formula ∫ λ × (∂q/∂p₁) dΩ is a volume integral
  # over the same domain. Using different meshes would introduce
  # interpolation errors in the gradient computation.
  [mesh]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 30
  []
[]

[Variables]
  # λ (adjoint variable, here named u_adj):
  # Solution of the adjoint screened Poisson equation.
  # Same function space as the forward variable u (FIRST Lagrange).
  # Self-adjointness of −∇² + I means L* = L, so the same FEM
  # discretisation is correct for both forward and adjoint.
  [u_adj]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # ------------------------------------------------------------------
  # Diffusion term: same as forward. Contributes −∇²λ (strong form).
  # ------------------------------------------------------------------
  [diffusion_adj]
    type     = Diffusion
    variable = u_adj
  []

  # ------------------------------------------------------------------
  # Screening (reaction) term: same coefficient α = 1.0 as forward.
  # Contributes +λ to the strong form: −∇²λ + λ = f_adj.
  # CoefReaction adds +coefficient × ∫ u_adj v dΩ (POSITIVE sign).
  # Combined with Diffusion: strong form is −∇²λ + λ = f_adj ✓
  # ------------------------------------------------------------------
  [screening_adj]
    type        = CoefReaction
    variable    = u_adj
    coefficient = 1.0
  []
[]

[DiracKernels]
  # ------------------------------------------------------------------
  # ReporterPointSource: adjoint RHS = Σᵢ eᵢ δ(x − xᵢ)
  #
  # This DiracKernel places weighted point sources at the 9 sensor
  # locations. At each location xᵢ, the weight eᵢ = u(xᵢ) − u_obs,ᵢ
  # is the misfit from the forward solve, received from the main app.
  #
  # The DiracKernel adds −∫ eᵢ δ(x−xᵢ) v dΩ = −eᵢ v(xᵢ) to the residual
  # (BodyForce sign convention). With residual = 0:
  #   −∇²λ + λ − Σᵢ eᵢ δ(x−xᵢ) = 0
  #   ⟺  −∇²λ + λ = Σᵢ eᵢ δ(x−xᵢ)   ✓
  #
  # The x,y,z coordinates and misfit values are read from the
  # ConstantReporter 'misfit_data', which is filled by [to_adjoint]
  # Transfer from the main app before each ADJOINT execute_on call.
  # ------------------------------------------------------------------
  [misfit_source]
    type         = ReporterPointSource
    variable     = u_adj
    x_coord_name = 'misfit_data/measurement_xcoord'
    y_coord_name = 'misfit_data/measurement_ycoord'
    z_coord_name = 'misfit_data/measurement_zcoord'
    value_name   = 'misfit_data/misfit_values'
  []
[]

[BCs]
  # Homogeneous Dirichlet BCs: λ = 0 on all four walls.
  # The adjoint BCs match the forward BCs for a self-adjoint operator.
  # Physical interpretation: the adjoint function λ represents the
  # sensitivity of J to a point source perturbation at (x,y). Since
  # the forward solution u = 0 on ∂Ω, the adjoint is also zero there
  # (boundary sensors would contribute nothing additional anyway, since
  # all 9 sensors are interior points).
  [all_walls]
    type     = DirichletBC
    variable = u_adj
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # grad_basis_fn: the basis function shape ∂q/∂p₁ = sin(πx) sin(πy)
  #
  # This is the derivative of the source q = p₁ × sin(πx) sin(πy)
  # with respect to the parameter p₁. It is independent of p₁ and
  # needs to be evaluated only once per mesh (not per iteration).
  #
  # The gradient ∫ λ × grad_basis_fn dΩ is computed by the
  # ElementOptimizationSourceFunctionInnerProduct reporter.
  #
  # For a multi-parameter problem with q = Σⱼ pⱼ φⱼ(x,y), one would
  # need one grad_basis_fn per parameter (∂q/∂pⱼ = φⱼ(x,y)), and
  # one ElementOptimizationSourceFunctionInnerProduct per parameter.
  # ------------------------------------------------------------------
  [grad_basis_fn]
    type               = ParsedOptimizationFunction
    expression         = 'p0 * sin(pi * x) * sin(pi * y)'
    param_symbol_names = 'p0'
    param_vector_name  = 'params/values'
  []
[]

[Reporters]
  # ------------------------------------------------------------------
  # ConstantReporter 'misfit_data': receives misfit and sensor locations.
  #
  # The main app transfers:
  #   - Sensor coordinates (xᵢ, yᵢ, zᵢ) from OptimizationReporter
  #   - Misfit values eᵢ from the forward app (via OptimizationReporter)
  # These are stored here and consumed by ReporterPointSource DiracKernel.
  #
  # Initial values are placeholder zeros (9 each for coords, 9 for misfits).
  # The main app always transfers updated values before the adjoint solves.
  # ------------------------------------------------------------------
  [misfit_data]
    type                  = ConstantReporter
    real_vector_names     = 'measurement_xcoord
                             measurement_ycoord
                             measurement_zcoord
                             misfit_values'
    real_vector_values    = '0.25 0.50 0.75 0.25 0.50 0.75 0.25 0.50 0.75;
                             0.25 0.25 0.25 0.50 0.50 0.50 0.75 0.75 0.75;
                             0 0 0 0 0 0 0 0 0;
                             0 0 0 0 0 0 0 0 0'
  []

  # Parameter reporter for the ParsedOptimizationFunction.
  # Receives the current parameter value from the main optimizer app.
  [params]
    type               = ConstantReporter
    real_vector_names  = 'values'
    real_vector_values = '1.0'
  []
[]

[VectorPostprocessors]
  # ------------------------------------------------------------------
  # ElementOptimizationSourceFunctionInnerProduct 'gradient':
  #
  # Computes the gradient of J with respect to the source amplitude p₁:
  #
  #   dJ/dp₁ = ∫ u_adj(x,y) × sin(πx) sin(πy) dΩ
  #
  # ElementOptimizationSourceFunctionInnerProduct is a VectorPostprocessor
  # that returns a vector of inner products, one per parameter. For our
  # single parameter (p₁), it returns the vector [dJ/dp₁].
  #
  # The main app reads this via [from_adjoint] Transfer referencing
  # 'gradient/inner_product' and stores it in
  # OptimizationReporter/grad_source_amplitude for the TAO L-BFGS step.
  #
  # Expected gradient at the starting guess p₁ = 1.0:
  #   u₁(x,y) = (1/(2π²+1)) sin(πx) sin(πy) ≈ 0.04822 sin(πx) sin(πy)
  #   misfit: eᵢ = u₁(xᵢ) − u_obs,ᵢ ≈ (0.04822 − 0.48218) × sin products
  #                                    ≈ −0.43396 × sin products (negative)
  #   adjoint RHS: Σᵢ eᵢ δ(x−xᵢ)  (negative point loads at sensors)
  #   adjoint solution: λ(x,y) has negative values near sensors
  #   gradient: dJ/dp₁ = ∫ λ sin(πx)sin(πy) dΩ < 0
  #
  # Negative gradient means increasing p₁ decreases J → the optimizer
  # correctly increases p₁ from 1.0 toward the true value 10.0. ✓
  # ------------------------------------------------------------------
  [gradient]
    type     = ElementOptimizationSourceFunctionInnerProduct
    variable = u_adj
    function = grad_basis_fn
  []
[]

[Executioner]
  # Steady-state Newton solve for the linear adjoint system.
  # Same solver setup as the forward app. The adjoint system has
  # the same SPD structure as the forward (self-adjoint operator),
  # so LU is fast and reliable.
  type       = Steady
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-14
[]

[Outputs]
  # Output the adjoint field for diagnostics:
  #   - Exodus: visualize λ(x,y) to verify the sensitivity structure.
  #     The adjoint should peak near the sensor locations and decay away.
  #     Near a sensor at (xᵢ, yᵢ), λ looks like a Green's function of
  #     the screened Poisson operator — a rounded "bump" with width ~ 1/√α.
  #   - CSV: inner_product (dJ/dp₁) at each optimizer iteration.
  #     The gradient should be large and negative in iteration 1,
  #     then near zero after the optimizer finds p₁ = 10.
  exodus = true
  csv    = true
[]
