# ============================================================
# Case 23: Magnetic Diffusion into a Conducting Slab
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 6 §6.2
#
# A step-applied magnetic field diffuses into a conductor:
#   ∂B/∂t = D_m · ∇²B      (magnetic diffusion equation)
#
# Mathematically identical to the heat equation (Case 03).
# Analytical solution: B(x,t) = erfc(x / 2√(D_m·t))
#
# Domain: [0,1] × [0,0.04] quasi-1D slab, 50×2 mesh
# IC: B = 0 (field-free conductor at t=0)
# BC: B = 1 at x=0 (step-applied field), B = 0 at x=1 (far field)
# Parameters: D_m = 0.01  (magnetic diffusivity = 1/(μ₀σ))
# ============================================================

# HIT top-level variable: magnetic diffusivity D_m = 1/(μ₀·σ)
# For copper: μ₀ = 4π×10⁻⁷ H/m, σ ≈ 6×10⁷ S/m → D_m ≈ 0.013 m²/s
# Here we use D_m = 0.01 (dimensionless / normalised units).
D_m = 0.01

[Mesh]
  type = GeneratedMesh
  dim  = 2
  # 50 elements along x (the diffusion direction) for good spatial resolution.
  # 2 elements along y makes the domain quasi-1D while remaining 2D-valid for MOOSE.
  nx   = 50
  ny   = 2
  xmin = 0
  xmax = 1        # slab thickness = 1 (normalised)
  ymin = 0
  ymax = 0.04     # thin strip; physics is independent of y
[]

[Variables]
  # B is the magnetic flux density (Tesla, or normalised by B_applied = 1).
  # Initial condition B=0 everywhere: the conductor is field-free before t=0.
  [B]
    initial_condition = 0.0
  []
[]

[Kernels]
  # ∂B/∂t term — rate of change of magnetic flux density at each point.
  # ADTimeDerivative uses automatic differentiation for exact Jacobians.
  [B_time]
    type     = ADTimeDerivative
    variable = B
  []

  # D_m · ∇²B term — magnetic diffusion (identical in form to Fourier heat conduction).
  # ADMatDiffusion reads the diffusivity from the material property named 'diffusivity'.
  [B_diff]
    type        = ADMatDiffusion
    variable    = B
    diffusivity = diffusivity   # material property name (defined in [Materials])
  []
[]

[BCs]
  # Left boundary (x=0): step-applied external field B = 1 at t=0+.
  # In Melcher's notation this is the "skin depth" driving boundary condition.
  [B_left]
    type     = DirichletBC
    variable = B
    boundary = left
    value    = 1.0   # B_applied = 1 (normalised)
  []

  # Right boundary (x=1): far-field condition B = 0.
  # The conductor is thick enough that the field has not yet penetrated to x=1
  # for most of the simulation (penetration depth δ ≈ 2√(D_m·t) reaches 1 only
  # at t ≈ 25, well beyond end_time = 20).
  [B_right]
    type     = DirichletBC
    variable = B
    boundary = right
    value    = 0.0
  []
[]

[Materials]
  # Constant magnetic diffusivity D_m = 1/(μ₀·σ).
  # ADGenericConstantMaterial is the AD-compatible version of GenericConstantMaterial.
  # The property name 'diffusivity' is the default name expected by ADMatDiffusion.
  [mag_diff]
    type        = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '${D_m}'   # magnetic diffusivity from top-level variable
  []
[]

[Postprocessors]
  # Spatial average of B over the whole slab — rises monotonically as the field
  # penetrates.  At steady state (t→∞) the profile becomes linear: avg_B → 0.5.
  [avg_B]
    type     = ElementAverageValue
    variable = B
  []

  # Maximum B anywhere in the domain.  Should stay at 1.0 (left-boundary value)
  # once the first time step has been taken; confirms the BC is active.
  [max_B]
    type       = ElementExtremeValue
    variable   = B
    value_type = max
  []

  # Minimum B anywhere in the domain.  Decreases toward 0 at the right boundary.
  # Useful for checking that the far-field BC is being enforced.
  [min_B]
    type       = ElementExtremeValue
    variable   = B
    value_type = min
  []
[]

[Executioner]
  type = Transient

  # NEWTON with exact AD Jacobians — fast convergence even for nonlinear extensions.
  solve_type          = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Fixed timestep: dt = 0.2 gives 100 steps to t = 20.
  # Penetration depth at each output: δ(t) ≈ 2√(D_m·t) = 2√(0.01·t)
  #   t=2  → δ ≈ 0.28    t=10 → δ ≈ 0.63    t=20 → δ ≈ 0.89
  [TimeStepper]
    type = ConstantDT
    dt   = 0.2
  []

  start_time = 0.0
  end_time   = 20.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true   # full B(x,y,t) field at every timestep for ParaView animation
  csv    = true   # avg_B, max_B, min_B time-history for plotting
[]
