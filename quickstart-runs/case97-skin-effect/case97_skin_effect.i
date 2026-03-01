# ============================================================
# Case 97: Skin Effect — Magnetic Diffusion into Conducting Slab
# MIT 6.641, Lec 9 — Magnetic Diffusion
# Prof. Markus Zahn, Spring 2005
#
# An oscillating magnetic field H₀ sin(ωt) is applied to the
# surface (x = 0) of a conducting slab. The field diffuses
# into the conductor according to:
#
#   ∂H/∂t = D_m ∂²H/∂x²    where D_m = 1/(μσ)
#
# The sinusoidal steady-state solution is:
#   H(x,t) = H₀ exp(−x/δ) sin(ωt − x/δ)
#
# where the skin depth is:
#   δ = √(2D_m/ω) = √(2/(μσω))
#
# Applications: electromagnetic shielding, induction heating,
# eddy current braking, transformer lamination design.
#
# Domain: [0, 1] × [0, 0.04] quasi-1D slab, 100×2 mesh
# IC: H = 0 (field-free conductor)
# BC: H(0,t) = sin(ωt) at x = 0, H(L,t) = 0 at x = 1
# Parameters: D_m = 0.01, ω = 2π → δ = √(2·0.01/(2π)) ≈ 0.0564
# ============================================================

D_m   = 0.01       # magnetic diffusivity 1/(μσ)
omega = 6.2831853  # angular frequency ω = 2π (period T = 1 s)
# Skin depth δ = √(2·D_m/ω) = √(0.02/6.283) ≈ 0.0564

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100      # fine resolution along diffusion direction
  ny   = 2        # quasi-1D
  xmin = 0
  xmax = 1.0
  ymin = 0
  ymax = 0.04
[]

[Variables]
  # Magnetic field H [A/m] (normalised, H₀ = 1).
  [H]
    initial_condition = 0.0
  []
[]

[Kernels]
  # ∂H/∂t — rate of change of magnetic field.
  [H_time]
    type     = ADTimeDerivative
    variable = H
  []

  # D_m ∂²H/∂x² — magnetic diffusion.
  [H_diff]
    type        = ADMatDiffusion
    variable    = H
    diffusivity = diffusivity
  []
[]

[BCs]
  # Left boundary (x = 0): oscillating applied field H = sin(ωt).
  # This drives the skin effect — the field penetrates into the slab.
  [H_oscillating]
    type     = FunctionDirichletBC
    variable = H
    boundary = left
    function = H_applied
  []

  # Right boundary (x = 1): far-field H = 0.
  # The slab is thick compared to the skin depth (1/0.056 ≈ 18 skin depths).
  [H_far]
    type     = DirichletBC
    variable = H
    boundary = right
    value    = 0.0
  []
[]

[Functions]
  # Oscillating boundary field H₀ sin(ωt).
  [H_applied]
    type       = ParsedFunction
    expression = 'sin(${omega} * t)'
  []
[]

[Materials]
  # Magnetic diffusivity D_m = 1/(μσ).
  [mag_diff]
    type        = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '${D_m}'
  []
[]

[Postprocessors]
  # Field at x = δ ≈ 0.056: should have amplitude ≈ 1/e ≈ 0.368.
  [H_at_skin_depth]
    type     = PointValue
    variable = H
    point    = '0.056 0.02 0'
  []

  # Field at x = 2δ ≈ 0.113: should have amplitude ≈ e^{-2} ≈ 0.135.
  [H_at_2delta]
    type     = PointValue
    variable = H
    point    = '0.113 0.02 0'
  []

  # Field at x = 3δ ≈ 0.169: amplitude ≈ e^{-3} ≈ 0.050.
  [H_at_3delta]
    type     = PointValue
    variable = H
    point    = '0.169 0.02 0'
  []

  # Surface field at x = 0 (should follow sin(ωt)).
  [H_surface]
    type     = PointValue
    variable = H
    point    = '0.0 0.02 0'
  []

  # Domain average: approaches zero in SSS (symmetric oscillation).
  [avg_H]
    type     = ElementAverageValue
    variable = H
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Run for 5 periods (5 s) with 50 steps per period.
  # This allows the initial transient to die out and reach
  # sinusoidal steady state by about t = 3 s.
  dt       = 0.02
  end_time = 5.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
