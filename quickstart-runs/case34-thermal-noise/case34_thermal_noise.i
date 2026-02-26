# ============================================================
# Case 34: Thermal Noise Relaxation — Fluctuation-Dissipation
# Haus, Electromagnetic Noise and Quantum Optical Measurements
# (Springer, 2000), Ch. 5 (Nyquist theorem)
#
# Heat diffusion from random initial conditions — spectral relaxation.
#   ∂T/∂t = D · ∇²T
#
# Each spatial Fourier mode k decays independently:
#   T̂_k(t) = T̂_k(0) · exp(−D · |k|² · t)
# Short-wavelength (large |k|) modes decay rapidly;
# long-wavelength modes persist.
#
# Snapshot timeline on [0,1]²:
#   t = 0.0  — speckle pattern (random IC, min=0.3, max=0.7)
#   t ≈ 0.1  — fine-scale features erased (high-k modes decayed)
#   t ≈ 0.5  — only fundamental mode (k = 2π) survives
#   t = 2.0  — domain uniform at T ≈ 0.5 (equilibrium)
#
# Fundamental mode decay rate (L=1 square, lowest k = π):
#   γ₁ = D · π² · (1² + 1²) = 0.1 · 2π² ≈ 1.97 s⁻¹
# so the max–min temperature range decays as ~exp(−1.97 t) at late times.
#
# Domain: [0,1]², 30×30 mesh
# IC: T ~ Uniform(0.3, 0.7) random noise around equilibrium T_eq = 0.5
# BC: T = 0.5 on all walls (isothermal enclosure at equilibrium)
# Parameters: D = 0.1 m²/s (thermal diffusivity), dt = 0.01, t_end = 2.0
# ============================================================

# HIT top-level variable: thermal diffusivity D [m²/s].
# Increasing D speeds up all mode relaxations proportionally.
# The fundamental-mode decay time is τ₁ = 1/γ₁ = 1/(D·2π²) ≈ 0.51 s.
D = 0.1

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  # T is temperature, initialised with uniform random noise in [0.3, 0.7].
  # The equilibrium value is T_eq = 0.5 (the mean of the random IC interval).
  # RandomIC seeds a different (reproducible) random pattern for each element.
  [T]
    order  = FIRST
    family = LAGRANGE
    [InitialCondition]
      # RandomIC draws T uniformly from [min, max] at each node.
      # The deviation amplitude is ±0.2 about T_eq = 0.5.
      type = RandomIC
      min  = 0.3
      max  = 0.7
    []
  []
[]

[Kernels]
  # ∂T/∂t — accumulation term (Fourier spectral coefficient rate of change).
  # ADTimeDerivative uses automatic differentiation for exact Jacobian entries.
  [time_deriv]
    type     = ADTimeDerivative
    variable = T
  []

  # D · ∇²T — diffusion term, responsible for mode-selective decay.
  # ADMatDiffusion reads the diffusivity from the material property 'D_coeff'.
  # Each mode k decays at rate D·|k|², so high-frequency noise dies first.
  [diffusion]
    type        = ADMatDiffusion
    variable    = T
    diffusivity = D_coeff   # material property defined in [Materials]
  []
[]

[BCs]
  # T = T_eq = 0.5 on all four walls: isothermal enclosure at equilibrium.
  # This BC pins the boundary and drives the interior to relax toward 0.5.
  # The wall acts as a thermal reservoir at the Nyquist equilibrium temperature.
  [walls]
    type     = ADDirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0.5
  []
[]

[Materials]
  # Constant thermal diffusivity D = 0.1 m²/s.
  # ADGenericConstantMaterial is the AD-compatible form of GenericConstantMaterial.
  # The property name 'D_coeff' is referenced by the diffusion kernel above.
  [thermal]
    type        = ADGenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D}'   # thermal diffusivity from top-level variable
  []
[]

[Postprocessors]
  # Spatial average of T over the whole domain.
  # Starts near 0.5 (mean of uniform IC) and stays there — conservation check.
  # Small drift toward 0.5 is expected as the boundary BC pulls the interior.
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []

  # Maximum T anywhere in the domain.
  # Starts near 0.7 (IC max) and decays monotonically toward 0.5.
  # The gap (max_T − 0.5) measures the surviving high-T fluctuation amplitude.
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []

  # Minimum T anywhere in the domain.
  # Starts near 0.3 (IC min) and rises monotonically toward 0.5.
  # The gap (0.5 − min_T) measures the surviving low-T fluctuation amplitude.
  [min_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = min
  []

  # L2 norm of T over the domain: ‖T‖ = √(∫T² dV).
  # Because T → 0.5 (not 0), this approaches √0.25 = 0.5 at steady state.
  # The decay of fluctuations is better read from (max_T − min_T); see above.
  [L2_norm]
    type     = ElementL2Norm
    variable = T
  []
[]

[Executioner]
  type = Transient

  # NEWTON with exact AD Jacobians.
  # The linear diffusion problem converges in 1–2 Newton iterations per step.
  solve_type          = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Fixed timestep dt = 0.01 s, 200 steps to t = 2.0 s.
  # dt is well below the diffusion stability limit for the finest resolved mode:
  #   dt_crit ≈ h²/(4D) = (1/30)²/(0.4) ≈ 0.0028 s  (explicit reference only;
  #   implicit time integration has no such stability restriction).
  # By t = 2.0 s the fundamental mode has decayed by exp(−1.97×2) ≈ e⁻³·⁹ ≈ 0.02.
  dt       = 0.01
  end_time = 2.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

[Outputs]
  exodus = true   # full T(x,y,t) field at every timestep for ParaView animation
  csv    = true   # avg_T, max_T, min_T, L2_norm time-history for plotting
[]
