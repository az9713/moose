# ============================================================
# Case 31: Driven Resonant Cavity — Frequency Response and Q
# Haus, Electromagnetic Noise and Quantum Optical Measurements
#   (Springer, 2000), Ch. 3
#
# Helmholtz equation for the z-component of electric field
# inside a 2-D rectangular cavity driven by a localized current
# source (TM polarization):
#
#   ∇²E + k²·E = −J_source(x,y)      in Ω = [0,2]×[0,1]
#   E = 0                              on all PEC walls
#
# Eigenvalues of the bare cavity (J=0, E=0 walls) are
#   k²_mn = (mπ/a)² + (nπ/b)²
# For a=2, b=1 the TM₁₁ mode has
#   k²₁₁ = (π/2)² + π² = π²(1/4+1) = 5π²/4 ≈ 12.337
#
# Driving k² near k²₁₁ produces a large resonant field.
# Off resonance (e.g. k² = 15.0) the response is weak.
#
# MOOSE weak-form mapping
# -----------------------
#   ∫ ∇E·∇v dV          →  Diffusion kernel
#   −k²·∫ E·v dV        →  CoefReaction(coefficient = −k²)
#   −∫ J·v dV           →  BodyForce(function = source_fn)
#
# Sign check:
#   Diffusion contributes +∫ ∇E·∇v dV  (strong form: −∇²E)
#   CoefReaction(−k²) contributes −k²·∫ E·v dV  (strong form: −k²·E)
#   Sum → strong form  −∇²E − k²E = J_source  ⟺  ∇²E + k²E = −J_source  ✓
#
# To explore off-resonance response: change k_squared to 15.0
# and re-run.  The peak |E| drops by roughly two orders of magnitude.
# ============================================================

# HIT top-level variable — change this single value to sweep frequency.
# k_squared = 12.337 is the TM₁₁ resonance; 12.3 is slightly below
# to avoid an exact singular system while still showing large response.
k_squared = 12.3

[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    nx    = 40      # 40 elements in x (domain width = 2)
    ny    = 20      # 20 elements in y (domain height = 1)
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1
  []
[]

[Variables]
  # E is the z-component of electric field [V/m] (TM polarization).
  # PEC boundary condition requires E_z = 0 on all metallic walls.
  [E]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # Localized current source modelling a small antenna or probe.
  # A narrow Gaussian centred at (0.7, 0.4) with half-width √0.02 ≈ 0.14 m
  # couples efficiently to the TM₁₁ mode because the mode amplitude is
  # non-zero there (the mode has a single interior maximum).
  [source_fn]
    type       = ParsedFunction
    expression = 'exp(-((x-0.7)^2+(y-0.4)^2)/0.02)'
  []
[]

[Kernels]
  # −∇²E term: Diffusion contributes ∫ ∇E·∇v dV to the residual.
  # In the strong form this is −∇²E (after integration by parts).
  [laplacian]
    type     = Diffusion
    variable = E
  []

  # −k²·E term: CoefReaction with a NEGATIVE coefficient contributes
  # coeff·∫ E·v dV = −k²·∫ E·v dV.  In the strong form this is −k²·E.
  # Combined with the Diffusion kernel the full strong form is
  #   −∇²E − k²E = 0,  i.e.  ∇²E + k²E = 0 (homogeneous).
  # The BodyForce below adds the source to make it inhomogeneous.
  [reaction]
    type        = CoefReaction
    variable    = E
    coefficient = ${fparse -k_squared}   # must be negative: see sign check above
  []

  # Source term J_source: BodyForce adds −∫ f·v dV to the residual,
  # which in the strong form contributes +f.  The full equation becomes
  #   −∇²E − k²E = J_source  →  ∇²E + k²E = −J_source  ✓
  [source]
    type     = BodyForce
    variable = E
    function = source_fn
  []
[]

[BCs]
  # Perfect Electric Conductor (PEC) boundary condition: E_z = 0
  # on all four cavity walls.  This models a metallic enclosure.
  [pec_walls]
    type     = DirichletBC
    variable = E
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Postprocessors]
  # Peak field magnitude — the key diagnostic for resonance.
  # Near k²₁₁ ≈ 12.337 this is large; off resonance it is small.
  [max_E]
    type       = ElementExtremeValue
    variable   = E
    value_type = max
  []

  # Domain-average of E — measures the net driven response.
  # For a mode with mixed sign (TM₁₁ has one lobe), this is smaller
  # than the peak but still indicates the resonant enhancement.
  [avg_E]
    type     = ElementAverageValue
    variable = E
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # LU direct solver: ideal for a single-variable Helmholtz problem.
  # Near resonance the system matrix is nearly singular and iterative
  # solvers may struggle; LU handles it robustly.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true   # spatial field E(x,y) — visualise the cavity mode shape
  csv    = true   # max_E and avg_E — compare across k_squared sweeps
[]
