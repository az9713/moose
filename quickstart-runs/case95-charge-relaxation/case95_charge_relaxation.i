# ============================================================
# Case 95: Dielectric Relaxation Time
# MIT 6.641, Lec 7 — Ohmic Conduction / Charge Relaxation
# Prof. Markus Zahn, Spring 2005
#
# Free charge in a weakly conducting dielectric decays
# exponentially with the dielectric relaxation time:
#
#   ∂ρ/∂t + (σ/ε)ρ = 0   →   ρ(t) = ρ₀ exp(−t/τ_e)
#
# where τ_e = ε/σ is the dielectric relaxation time.
#
# This is the fundamental timescale for ESD protection,
# charge dissipation in insulating materials, and the
# transition from displacement-current to conduction-current
# dominated behaviour.
#
# Domain: unit square [0,1]², 30×30 mesh
# IC: Gaussian charge blob at centre
# BC: zero-flux (natural Neumann) — charge decays in place
# Parameters: ε = 8.854e-12 F/m, σ = 1e-10 S/m → τ_e = 0.0885 s
#             (normalised: σ/ε = 11.3 → τ_e ≈ 0.0885 s)
# ============================================================

# Dielectric relaxation time ratio σ/ε.
# For a typical insulating polymer: ε_r ≈ 3 → ε = 3 × 8.854e-12,
# σ ≈ 1e-10 S/m → τ_e = ε/σ = 2.66e-11 / 1e-10 = 0.266 s
# For a quick simulation we normalise: σ/ε = 5.0 → τ_e = 0.2 s
sigma_over_eps = 5.0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 30
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Variables]
  # Free charge density ρ [C/m³].
  # Governed by the local ODE ∂ρ/∂t + (σ/ε)ρ = 0 at every point.
  # No spatial diffusion — the relaxation is purely temporal.
  [rho]
    order  = FIRST
    family = LAGRANGE
    [InitialCondition]
      # Gaussian blob centred at (0.5, 0.5) with width 0.1.
      # Peak value of 1.0 provides a non-trivial spatial distribution.
      type     = FunctionIC
      function = 'exp(-((x-0.5)^2 + (y-0.5)^2) / 0.01)'
    []
  []
[]

[Kernels]
  # ∂ρ/∂t: time derivative (accumulation term).
  [rho_time]
    type     = TimeDerivative
    variable = rho
  []

  # (σ/ε)ρ: linear decay (reaction term).
  # CoefReaction adds ∫ coeff·ρ·ψ dV to the residual.
  # The coefficient is σ/ε = 5.0, so the decay timescale is τ = 1/5 = 0.2 s.
  [rho_decay]
    type        = CoefReaction
    variable    = rho
    coefficient = ${sigma_over_eps}
  []
[]

# No BCs needed: zero-flux (natural Neumann) is the default.
# Physically, charge cannot leave the domain — it decays in place.
[BCs]
[]

[Postprocessors]
  # Domain-average charge: should decay as exp(−t/τ_e) = exp(−5t).
  [avg_rho]
    type     = ElementAverageValue
    variable = rho
  []

  # Peak charge density: also decays as exp(−5t) since the PDE
  # is a pointwise ODE (no spatial coupling).
  [max_rho]
    type       = ElementExtremeValue
    variable   = rho
    value_type = max
  []

  # Total charge integral: ∫ ρ dV.
  # Should decay as Q₀ exp(−t/τ_e).
  [total_charge]
    type     = ElementIntegralVariablePostprocessor
    variable = rho
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # 50 steps over 1.0 s = 5τ_e.
  # By t = 5τ_e the charge has decayed to exp(−5) ≈ 0.67%.
  dt       = 0.02
  end_time = 1.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
