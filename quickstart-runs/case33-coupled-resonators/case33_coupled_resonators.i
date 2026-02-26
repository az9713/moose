# ============================================================
# Case 33: Coupled Resonator Beating — Energy Exchange
# Haus, Electromagnetic Noise and Quantum Optical Measurements
#   (Springer, 2000), Ch. 3
#
# Two weakly coupled resonator amplitude envelopes u and v
# satisfy the Haus coupled-mode equations:
#
#   ∂u/∂t = D·∇²u − γ·u + κ·v
#   ∂v/∂t = D·∇²v − γ·v − κ·u
#
# The diffusion term D·∇²u broadens the initial Gaussian packet.
# The decay term −γ·u drains total energy as exp(−γ·t).
# The coupling term +κ·v transfers power back and forth between
# the two resonators at the beating frequency.
#
# Analytical beating (spatially uniform, no diffusion):
#   u(t) = cos(κ·t)·exp(−γ·t)   (starts full, drains to v)
#   v(t) = sin(κ·t)·exp(−γ·t)   (starts empty, fills from u)
#
# With diffusion the spatial Gaussian spreads but the temporal
# beating is preserved in the domain-average sense:
#   avg_u and avg_v oscillate 90° out of phase with each other.
#   Beating period:    T = π/κ ≈ 1.047 s  (half-period for u→v→u)
#   Total energy decay: E(t) ∝ exp(−γ·t)
#
# Domain: [0,1]², 30×30 mesh
# IC: u = Gaussian at (0.5, 0.5), v = 0
# BC: Neumann (zero-flux) on all walls — natural BC, no explicit entry needed
# Parameters: D = 0.01, γ = 0.5, κ = 3.0
# ============================================================

# HIT top-level variables — change these to re-scale the physics.
# All three are referenced with ${} below; the parser will error
# on any unused top-level variable, so they must all appear.
D_coeff = 0.01   # diffusivity D [m²/s, normalised]: spreads the spatial envelope
gamma   = 0.5    # decay rate γ [1/s]: e-folding time τ_decay = 1/γ = 2 s
kappa   = 3.0    # coupling rate κ [1/s]: beating period T = π/κ ≈ 1.047 s

[Mesh]
  type = GeneratedMesh
  dim  = 2
  # 30×30 QUAD4 elements give h = 1/30 ≈ 0.033.
  # The initial Gaussian has width σ = √0.01 = 0.1, so ~3 elements per σ.
  nx   = 30
  ny   = 30
[]

[Variables]
  # u — amplitude envelope of the first resonator.
  # Initialised to a Gaussian centred at (0.5, 0.5) in [ICs] below.
  [u]
  []

  # v — amplitude envelope of the second resonator.
  # Starts at zero: all energy is in u at t = 0.
  [v]
  []
[]

[ICs]
  # Gaussian initial condition for u: peak value 1 at the domain centre.
  # Width parameter 0.02 gives a spatial 1/e radius of √0.02 ≈ 0.14,
  # so the blob is well-contained within the unit square.
  [u_ic]
    type     = FunctionIC
    variable = u
    function = 'exp(-((x-0.5)^2+(y-0.5)^2)/0.02)'
  []

  # v starts identically zero; it fills via the κ·u coupling term.
  [v_ic]
    type     = ConstantIC
    variable = v
    value    = 0.0
  []
[]

[Kernels]
  # ------------------------------------------------------------------
  # Equation for u:   ∂u/∂t = D·∇²u − γ·u + κ·v
  # ------------------------------------------------------------------

  # ∂u/∂t — temporal rate of change of the u envelope.
  # ADTimeDerivative provides automatic-differentiation Jacobian entries.
  [u_time]
    type     = ADTimeDerivative
    variable = u
  []

  # D·∇²u — diffusive spreading of the spatial envelope.
  # ADMatDiffusion reads the diffusivity from the material property
  # named 'D_mat', which is defined in [Materials] below.
  [u_diff]
    type        = ADMatDiffusion
    variable    = u
    diffusivity = D_mat
  []

  # −γ·u — linear decay: each resonator loses energy at rate γ.
  # ADReaction adds +∫ rate·u·φ_i dV to the residual, which represents
  # the weak form of +γ·u.  Because MOOSE residuals are written as
  # F = 0, the ∂u/∂t term already carries a −1 sign in the weak form,
  # so ADReaction with rate = γ correctly implements −γ·u in the PDE.
  [u_decay]
    type     = ADReaction
    variable = u
    rate     = ${gamma}
  []

  # +κ·v — energy coupling from v into u (the beating driver).
  # CoupledForce adds −∫ coef·v·φ_i dV to the residual, which in
  # MOOSE's residual convention corresponds to +κ·v on the left-hand
  # side.  The negative sign in the weak residual is absorbed by the
  # Newton update direction, so coef = kappa is the correct sign here.
  [u_coupling]
    type     = CoupledForce
    variable = u
    v        = v          # the MOOSE variable v defined in [Variables]
    coef     = ${kappa}
  []

  # ------------------------------------------------------------------
  # Equation for v:   ∂v/∂t = D·∇²v − γ·v + κ·u
  # (mirror of the u equation with the coupling direction reversed)
  # ------------------------------------------------------------------

  # ∂v/∂t — temporal accumulation for the v envelope.
  [v_time]
    type     = ADTimeDerivative
    variable = v
  []

  # D·∇²v — diffusive spreading; same diffusivity as u (identical resonators).
  [v_diff]
    type        = ADMatDiffusion
    variable    = v
    diffusivity = D_mat
  []

  # −γ·v — decay at the same rate γ as u.
  [v_decay]
    type     = ADReaction
    variable = v
    rate     = ${gamma}
  []

  # −κ·u — antisymmetric coupling from u into v.
  # The coupling MUST be antisymmetric (+κ in the u equation, −κ in the v
  # equation) to produce oscillatory energy exchange (beating).  Symmetric
  # coupling (+κ in both) gives real eigenvalues −γ±κ and exponential
  # growth/decay instead of oscillation.
  # Eigenvalues of the coupling matrix [-γ, +κ; -κ, -γ] are λ = -γ ± jκ,
  # giving decaying oscillations at the beat frequency 2κ.
  [v_coupling]
    type     = CoupledForce
    variable = v
    v        = u          # couples u into the v equation
    coef     = ${fparse -kappa}
  []
[]

[Materials]
  # Constant diffusivity D shared by both resonator envelopes.
  # ADGenericConstantMaterial produces the AD-compatible property 'D_mat'
  # that ADMatDiffusion requires for exact Jacobian assembly.
  [diffusivity]
    type        = ADGenericConstantMaterial
    prop_names  = 'D_mat'
    prop_values = '${D_coeff}'
  []
[]

[Postprocessors]
  # Domain-average of u.  Should start near 1/A·∫exp(−r²/0.02) dA ≈ 0.063
  # (integral of the Gaussian over [0,1]²) and then oscillate while decaying.
  # avg_u and avg_v beat 90° out of phase with period T = π/κ ≈ 1.047 s.
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []

  # Domain-average of v.  Starts at zero, rises to a maximum near T/2 ≈ 0.52 s
  # (when u has transferred half its energy to v), then falls back.
  [avg_v]
    type     = ElementAverageValue
    variable = v
  []

  # Total integrated energy in u over the domain.
  # For the pure beating case (no diffusion, γ = 0) this would be
  #   ∫u dA = cos²(κ·t)·∫u₀ dA
  # With decay: E_u(t) ≈ cos²(κ·t)·exp(−γ·t)·∫u₀ dA.
  # Monitoring this confirms the beating and exponential envelope.
  [total_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = u
  []
[]

[Executioner]
  type = Transient

  # NEWTON with exact AD Jacobians: fast quadratic convergence per step.
  # Full SMP preconditioning (below) captures the u–v off-diagonal coupling.
  solve_type          = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # dt = 0.02 s gives ~52 steps per beating period T ≈ 1.047 s.
  # end_time = 3.0 s covers ~2.9 beating periods and 1.5 decay e-folds
  # (γ·t_end = 1.5 → energy reduces to exp(−1.5) ≈ 22% of initial).
  dt       = 0.02
  end_time = 3.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

# Full single-matrix preconditioning assembles the complete Jacobian
# including the off-diagonal blocks linking u and v.  This is essential
# for tight convergence of the coupled resonator system.
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true   # full u, v spatial fields at every timestep for ParaView
  csv    = true   # avg_u, avg_v, total_energy vs. time for beating verification
[]
