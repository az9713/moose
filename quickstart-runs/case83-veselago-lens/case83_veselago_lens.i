# ============================================================
# Case 83: Veselago Flat Lens — Point Source Focusing
# Veselago, Soviet Physics Uspekhi 10 (1968) 509–514
# Pendry, Physical Review Letters 85 (2000) 3966–3969
#
# Lectures: 1 (Left-Handed Materials, Perfect Lens)
#           10 (Negative Refraction in Photonic Crystals)
#
# ============================================================
# HISTORICAL AND PHYSICAL BACKGROUND
# ============================================================
#
# In 1968 Victor Veselago predicted that a medium with BOTH
# ε_r < 0 and μ_r < 0 simultaneously would support
# electromagnetic wave propagation with a NEGATIVE refractive
# index:
#
#   n = -sqrt(ε_r μ_r)   (choosing the negative root for LHM)
#
# For an ideal material with ε_r = -1, μ_r = -1 we have n = -1.
#
# The most striking prediction: a flat slab of n = -1 material
# of thickness d acts as a PERFECT LENS. A point source at
# distance s < d from the first face of the slab is focused to
# an IMAGE inside the slab (at depth s from the entry face) and
# again outside the slab (at distance d - s from the exit face).
#
# In 2000, John Pendry showed that the "perfect lens" is not
# merely a geometric optics result — it applies even to
# evanescent fields. For an ideal lossless n = -1 slab, the
# transmission coefficient for any wave vector k is exactly 1
# in amplitude, meaning ALL spatial frequencies (including
# evanescent ones) are reconstructed in the image. This would
# give resolution beyond the classical diffraction limit.
#
# In practice:
#   (1) Loss regularises the divergence — the perfect lens
#       requires Im(ε_r) → 0, and any non-zero loss limits
#       the maximum spatial frequency that is restored.
#   (2) The flat lens requires the source within one wavelength
#       of the slab face (s < λ₀) to capture evanescent fields.
#   (3) Real metamaterials (split-ring resonators) have both
#       non-zero loss and dispersion, limiting performance.
#
# This case demonstrates the FOCUSING EFFECT for propagating
# (non-evanescent) waves, which is already visible in the
# geometric optics limit (Snell's law with n = -1) and is
# captured by the frequency-domain Helmholtz FEM.
#
# ============================================================
# GEOMETRY
# ============================================================
#
#        x = -6        x = -1       x = +1       x = +6
#          |           |<-- LHM -->|              |
#          |           |   n = -1  |              |
#          |           |           |              |
#          |           |           |              |
#  source *             slab d=2    image (ext)
#  at (-2,0)           enters here  at (+2,0)
#                      internal focus at (0,0)
#
# Domain: [-6, 6] × [-4, 4]  (12 × 8 units = 3λ₀ × 2λ₀ with λ₀ = 4)
# LHM slab: -1 < x < +1  (thickness d = 2 = λ₀/2)
# Point source: (-2, 0)   (distance s = 1 from slab face)
#
# Predicted focus locations for ideal n = -1, d = 2, s = 1:
#   Internal focus:  x = -1 + s = 0   (inside slab, at slab centre)
#   External focus:  x = +1 + (d - s) = +2   (outside slab)
#
# ============================================================
# GOVERNING EQUATION — 1/μ_r FORMULATION
# ============================================================
#
# For TE polarisation (E_z out of the 2D plane, E = E_z ẑ):
#
#   ∇·(1/μ_r ∇E_z) + k₀² ε_r E_z = −jωμ₀ J_z
#
# This is NOT the same as the simple Helmholtz ∇²E + k₀² n² E = source,
# because n² = ε_r μ_r is the same for n = -1 and n = +1 (both give
# n² = 1 when ε_r = μ_r = ±1). The ONLY way to distinguish negative-
# index from positive-index in the Helmholtz equation is through the
# 1/μ_r coefficient in the diffusion term:
#
#   Vacuum:    1/μ_r = +1   →  normal Laplacian
#   LHM slab:  1/μ_r = -1   →  REVERSED sign on ∇² term
#
# The interface condition that follows from this weak form is:
#
#   Continuity: E_z continuous,  (1/μ_r) ∂E_z/∂n continuous
#
# At the vacuum→LHM interface:
#   (1/1) × dE/dx|vac = (1/(-1)) × dE/dx|LHM
#   ⟹  dE/dx|LHM = −dE/dx|vac   (SIGN REVERSAL)
#
# This sign reversal is the mathematical signature of negative
# refraction: the phase velocity inside the LHM is in the OPPOSITE
# direction to the energy velocity, bending the wavefront backward.
#
# The same 1/μ_r formulation was used in Case 74 (1D LHM slab).
# Here we extend it to 2D with a point source.
#
# ============================================================
# LOSS REGULARISATION
# ============================================================
#
# An ideal lossless n = -1 lens is a singular problem: for a
# perfect impedance match (ε_r = -1, μ_r = -1 exactly), the
# transmission coefficient diverges for evanescent waves, and
# the FEM system matrix becomes singular or nearly singular.
#
# We regularise by adding a small imaginary part to ε_r:
#
#   ε_r = -1 + j × δ   with δ = 0.3
#
# This gives:
#   ε_r' = Re(ε_r) = -1    (LHM condition maintained)
#   ε_r'' = Im(ε_r) = 0.3   (moderate loss → well-conditioned system)
#
# The loss couples the real and imaginary parts of E_z:
#
# E_r equation:
#   ∇·(1/μ_r ∇E_r) + k₀² ε_r' E_r + k₀² ε_r'' E_i = 0
#
# E_i equation:
#   ∇·(1/μ_r ∇E_i) + k₀² ε_r' E_i − k₀² ε_r'' E_r = source
#
# (The source −jωμ₀J has only an imaginary part in the phasor
#  convention: RHS of E_r = 0, RHS of E_i = source amplitude.)
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# λ₀ = 4 units  →  k₀ = 2π/λ₀ = π/2 = 1.5707963... rad/unit
# k₀² = (π/2)² = π²/4 = 2.4674011... unit⁻²
#
# Vacuum (|x| > 1):
#   1/μ_r       = +1.0
#   k₀² ε_r'    = +2.4674011   (ε_r = +1)
#   k₀² ε_r''   =  0           (no loss)
#
# LHM slab (-1 < x < +1):
#   1/μ_r       = -1.0
#   k₀² ε_r'    = -2.4674011   (ε_r = -1)
#   k₀² ε_r''   = +0.7402203   (δ = 0.3, so k₀²×0.3 = 0.74022)
#
# Point source Gaussian:
#   σ = 0.15 (narrowband, ~ λ₀/27)
#   A = 50   (amplitude chosen to give visible field levels)
#   J(x,y) = A × exp(−((x+2)² + y²) / (2σ²))
#           = 50 × exp(−((x+2)² + y²) / 0.045)
#
# ============================================================
# WEAK FORM AND KERNEL SIGN CONVENTIONS
# ============================================================
#
# MOOSE kernel sign conventions (reproduced here for clarity):
#
#   ADMatDiffusion(diffusivity = D):
#     Residual: +∫ D ∇E · ∇v dV
#     Strong form equivalent: −∇·(D ∇E)
#
#   ADMatReaction(reaction_rate = r):
#     Residual: −r × ∫ E v dV
#     Strong form equivalent: −r × E
#
#   ADMatCoupledForce(v = u, mat_prop_coef = c):
#     Residual: −c × ∫ u × test dV
#     Strong form equivalent: −c × u
#
#   BodyForce(function = f):
#     Residual: −∫ f v dV
#     Strong form equivalent: −f (i.e. adds −f to the LHS)
#
# Assembling E_r equation:
#   ADMatDiffusion(inv_mu_r)   →  −∇·(1/μ_r ∇E_r)
#   ADMatReaction(k0sq_eps_pr) →  −k₀² ε_r' E_r
#   ADMatCoupledForce(E_i, k0sq_eps_pp) →  −k₀² ε_r'' E_i
#   Total strong form: −∇·(1/μ_r ∇E_r) − k₀²ε_r'E_r − k₀²ε_r''E_i = 0
#   ↔  ∇·(1/μ_r ∇E_r) + k₀²ε_r'E_r + k₀²ε_r''E_i = 0  ✓
#
# Assembling E_i equation:
#   ADMatDiffusion(inv_mu_r)   →  −∇·(1/μ_r ∇E_i)
#   ADMatReaction(k0sq_eps_pr) →  −k₀² ε_r' E_i
#   ADMatCoupledForce(E_r, neg_k0sq_eps_pp) → +k₀² ε_r'' E_r
#   BodyForce(source_fn) → −f(x,y) [strong form adds -f to LHS]
#   Total: −∇·(1/μ_r ∇E_i) − k₀²ε_r'E_i + k₀²ε_r''E_r = −f
#   ↔  ∇·(1/μ_r ∇E_i) + k₀²ε_r'E_i − k₀²ε_r''E_r = +f  ✓
#
# The BodyForce adds −∫f·v to the residual, which in strong form
# means the equation LHS = f. For source f > 0, this represents
# the −jωμ₀J source with f acting as the source amplitude.
#
# ============================================================
# BOUNDARY CONDITIONS (Absorbing Robin BC)
# ============================================================
#
# All four boundaries use the EMRobinBC with profile_func_real = 0
# (absorbing only — no incident wave injection). The source is
# an interior point source, not an incident plane wave from a port.
#
# The Robin ABC is:  ∂E/∂n + jk₀E = 0
# This absorbs outgoing cylindrical waves reasonably well for
# domains that are several wavelengths from the source/focus.
#
# The 6-unit clearance from the source to the left boundary
# (left boundary at x = -6, source at x = -2) is 4 units = 1 λ₀
# of clearance, adequate for a first-order ABC.
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
k     = 1.5707963267948966   # k₀ = π/2 rad/unit (λ₀ = 4 units)
E0    = 0                    # no incident wave from port (source is interior)
theta = 0                    # incidence angle placeholder for EMRobinBC

[Mesh]
  # 2D domain [-6, 6] × [-4, 4].
  # 120 × 80 = 9,600 quad elements.
  # Spatial resolution: 120 elements over 12 units = 10 elements/unit.
  # Free-space wavelength λ₀ = 4 units → 40 elements per λ₀ (well resolved).
  # Inside the LHM slab λ_LHM = λ₀/|n| = 4/1 = 4 units → same resolution.
  # The Gaussian source with σ = 0.15 is resolved by ~1.5 elements per σ
  # on this mesh; this is marginal but sufficient to capture the gross
  # focusing effect. A finer mesh (240×160) would give a crisper focus.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 120
    ny   = 80
    xmin = -6
    xmax = 6
    ymin = -4
    ymax = 4
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E_z = E_real + j E_imag.
  # For the point source as written (−jωμ₀J source), the source lies entirely
  # in the E_imag equation, so E_real will be driven by cross-coupling through
  # the slab. In the lossless limit (δ → 0) E_real would be zero and E_imag
  # would carry all the field energy; with loss both components are non-zero.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part. The source (BodyForce) acts directly on E_imag,
  # so the primary field pattern appears in E_imag. The focusing in the slab
  # shows most clearly in E_magnitude_sq = E_real² + E_imag².
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # inv_mu_fn: 1/μ_r(x,y) — diffusion coefficient
  #
  # This is the single most important function for the flat-lens effect.
  # In the vacuum regions (|x| > 1):  1/μ_r = +1.0  (normal diffusion)
  # In the LHM slab  (-1 < x < +1):   1/μ_r = -1.0  (reversed diffusion)
  #
  # The NEGATIVE diffusion coefficient in the slab is what makes Snell's
  # law give a negative refraction angle. Combined with ε_r = -1, the
  # impedance Z = sqrt(μ_r/ε_r) = sqrt((-1)/(-1)) = +1 matches vacuum
  # impedance, so there is (in the lossless limit) no reflection at the
  # vacuum–LHM interface — all power is transmitted into the slab.
  # ------------------------------------------------------------------
  [inv_mu_fn]
    type       = ParsedFunction
    expression = 'if(x > -1.0 & x < 1.0, -1.0, 1.0)'
  []

  # ------------------------------------------------------------------
  # k0sq_eps_pr_fn: k₀² × Re(ε_r(x,y)) — real part reaction coefficient
  #
  # Vacuum (|x| > 1):  k₀² × (+1) = +2.4674011
  # LHM slab:          k₀² × (-1) = -2.4674011
  #
  # The negative value in the slab is consistent with n = -1 physics.
  # ADMatReaction uses residual = −rate × E, so rate = +k₀²ε_r' gives
  # the correct Helmholtz contribution −k₀²ε_r'E in both regions.
  # ------------------------------------------------------------------
  [k0sq_eps_pr_fn]
    type       = ParsedFunction
    expression = 'if(x > -1.0 & x < 1.0, -2.4674011002723395, 2.4674011002723395)'
  []

  # ------------------------------------------------------------------
  # k0sq_eps_pp_fn: +k₀² × Im(ε_r(x,y)) — imaginary cross-coupling
  #
  # Applied to the E_real equation via ADMatCoupledForce(v = E_imag).
  # In vacuum (|x| > 1):  ε_r'' = 0    →  coefficient = 0
  # In LHM slab:           ε_r'' = 0.3  →  coefficient = k₀² × 0.3
  #                                        = 2.4674011 × 0.3 = 0.74022033
  #
  # This term couples E_imag into E_real in the lossy slab region.
  # Loss regularises the near-perfect transmission at the LHM interfaces,
  # preventing the system matrix from becoming singular.
  # ------------------------------------------------------------------
  [k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(x > -1.0 & x < 1.0, 0.74022033, 0.0)'
  []

  # ------------------------------------------------------------------
  # neg_k0sq_eps_pp_fn: −k₀² × Im(ε_r(x,y)) — negative imaginary cross-coupling
  #
  # Applied to the E_imag equation via ADMatCoupledForce(v = E_real).
  # The sign is reversed compared to the E_real equation, reflecting
  # the conjugate relationship between the two field components.
  # This asymmetric cross-coupling (±ε_r'') keeps the system
  # well-posed and arises from splitting E = E_r + jE_i with
  # complex ε_r = ε_r' + jε_r''.
  # ------------------------------------------------------------------
  [neg_k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(x > -1.0 & x < 1.0, -0.74022033, 0.0)'
  []

  # ------------------------------------------------------------------
  # source_fn: Gaussian point source at (-2, 0)
  #
  # Models a localised current source J_z(x,y):
  #   J(x,y) = A × exp(−((x+2)² + y²) / (2σ²))
  # with A = 50, σ = 0.15, so 2σ² = 0.045.
  #
  # The source is located at (-2, 0), a distance s = 1 unit from
  # the left face of the LHM slab at x = -1.
  # σ = 0.15 ≈ λ₀/27 ≈ 1.5 mesh elements — a nearly point-like source.
  #
  # The source enters only the E_imag equation (see sign convention above).
  # BodyForce adds −∫f·v to the residual; the strong-form equation
  # E_imag equation then has source = +f on the RHS:
  #   ∇·(1/μ_r ∇E_i) + k₀²ε_r'E_i − k₀²ε_r''E_r = f(x,y)
  # This corresponds to the imaginary part of −jωμ₀J with f = ωμ₀|J|.
  # ------------------------------------------------------------------
  [source_fn]
    type       = ParsedFunction
    expression = '50.0 * exp(-((x + 2.0) * (x + 2.0) + y * y) / 0.045)'
  []

  # cosTheta = 1.0 for normal incidence (used by EMRobinBC).
  # Since all four boundary EMRobinBCs are pure absorbing (no injection),
  # cosTheta appears in the k_eff = k₀ cos(θ) term of the ABC.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ------------------------------------------------------------------
  # Wrap the four ParsedFunctions as AD material properties.
  # ADGenericFunctionMaterial evaluates each function at every quadrature
  # point and returns an ADReal, which enables automatic differentiation
  # for the Newton Jacobian. This is critical for convergence:
  # the Jacobian of ADMatDiffusion and ADMatReaction depends on the
  # material property values, which vary discontinuously at x = ±1.
  # ------------------------------------------------------------------

  # 1/μ_r(x,y): diffusion coefficient for ADMatDiffusion
  # Negative in the LHM slab — this is the key that makes the flat lens work.
  [inv_mu_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'inv_mu_r'
    prop_values = 'inv_mu_fn'
  []

  # k₀² Re(ε_r(x,y)): reaction coefficient for both field equations
  [k0sq_eps_pr_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pr'
    prop_values = 'k0sq_eps_pr_fn'
  []

  # +k₀² Im(ε_r(x,y)): cross-coupling from E_imag into E_real equation
  [k0sq_eps_pp_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pp'
    prop_values = 'k0sq_eps_pp_fn'
  []

  # −k₀² Im(ε_r(x,y)): cross-coupling from E_real into E_imag equation
  [neg_k0sq_eps_pp_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'neg_k0sq_eps_pp'
    prop_values = 'neg_k0sq_eps_pp_fn'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇·(1/μ_r ∇E_r) + k₀² ε_r' E_r + k₀² ε_r'' E_i = 0
  #
  # Three kernels: ADMatDiffusion + ADMatReaction + ADMatCoupledForce
  # No BodyForce for E_real (source is purely imaginary in phasor form).
  # ==================================================================

  # Diffusion term: −∇·(1/μ_r ∇E_r) in strong form.
  # With inv_mu_r = +1 in vacuum: standard Laplacian.
  # With inv_mu_r = -1 in LHM: REVERSED Laplacian.
  # This sign reversal is the FEM implementation of negative μ_r.
  # The weak form naturally enforces (1/μ_r) ∂E/∂n continuity at the
  # vacuum–LHM interfaces at x = ±1 without any special interface kernels.
  [diff_real]
    type        = ADMatDiffusion
    variable    = E_real
    diffusivity = inv_mu_r
  []

  # Reaction term: −k₀² ε_r'(x,y) E_r in strong form.
  # In vacuum: −k₀² E_r (positive k² for wave propagation).
  # In LHM:   +k₀² E_r (negative k₀²ε_r' = negative × negative = positive).
  # The reversed sign of k₀²ε_r' in the LHM is the wave-equation
  # counterpart of negative index of refraction.
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling: −k₀² ε_r''(x,y) E_i in strong form.
  # This term is non-zero only inside the lossy LHM slab (|x| < 1).
  # In the lossless limit (ε_r'' → 0) this kernel contributes nothing
  # and E_real/E_imag decouple in the bulk (coupling only through BCs).
  # With small loss δ = 0.3, the coupling regularises
  # the near-singular LHM transmission.
  [cross_real]
    type          = ADMatCoupledForce
    variable      = E_real
    v             = E_imag
    mat_prop_coef = k0sq_eps_pp
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇·(1/μ_r ∇E_i) + k₀² ε_r' E_i − k₀² ε_r'' E_r = f(x,y)
  #
  # Three volume kernels + one source BodyForce.
  # The source acts only on E_imag because the current source J is real
  # and the phasor source −jωμ₀J has only an imaginary component.
  # ==================================================================

  # Diffusion term for E_imag: same 1/μ_r coefficient as E_real.
  # The negative diffusion in the LHM slab applies equally to both
  # field components — the slab material properties are the same.
  [diff_imag]
    type        = ADMatDiffusion
    variable    = E_imag
    diffusivity = inv_mu_r
  []

  # Reaction term for E_imag: same coefficient k₀²ε_r'(x,y).
  # Both equations share the same diagonal Helmholtz operator.
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling: +k₀² ε_r''(x,y) E_r in strong form.
  # The negative sign of neg_k0sq_eps_pp = −k₀²ε_r'' means
  # ADMatCoupledForce contributes −(−k₀²ε_r'')×E_r = +k₀²ε_r''×E_r,
  # which is the correct sign for the E_imag equation.
  [cross_imag]
    type          = ADMatCoupledForce
    variable      = E_imag
    v             = E_real
    mat_prop_coef = neg_k0sq_eps_pp
  []

  # Point source: BodyForce on E_imag only.
  # Residual: −∫ source_fn × test dV
  # Strong form: E_imag equation RHS = +source_fn
  # The Gaussian source at (−2, 0) with σ = 0.15 drives cylindrical
  # waves that propagate outward, enter the LHM slab at x = -1,
  # and are refocused by the negative-index material.
  [source_imag]
    type     = BodyForce
    variable = E_imag
    function = source_fn
  []
[]

[BCs]
  # ==================================================================
  # Absorbing boundary conditions on all four domain faces.
  #
  # EMRobinBC (Jin, "FEM in Electromagnetics", 3rd Ed., Eq. 9.60):
  #   ∂E/∂n + j k₀ cos(θ) E = 0   (absorbing only, no injection)
  #
  # profile_func_real = 0: sets the incident wave amplitude to zero
  # on all boundaries. The cylindrical waves from the interior source
  # are absorbed without significant spurious reflections.
  #
  # The domain is 3λ₀ wide (x from -6 to 6) and 2λ₀ tall (y from -4
  # to 4). The source at (-2, 0) and the expected focus at (+2, 0)
  # are each 4 units (1 λ₀) from the nearest boundary. The first-order
  # ABC has reflection coefficient |R| = |cos(θ)-1|/|cos(θ)+1| which
  # is small at near-normal incidence. For cylindrical waves from a
  # point source, the angular spread increases near the source, but
  # at 1λ₀ clearance the outgoing waves are approximately plane-wave
  # like and the first-order ABC is adequate.
  #
  # All 8 BC blocks use sign = negative, consistent with Cases 74–77.
  # ==================================================================

  # --- Left boundary (x = -6) ---
  [abc_real_left]
    type              = EMRobinBC
    variable          = E_real
    boundary          = left
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_left]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = left
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Right boundary (x = +6) ---
  [abc_real_right]
    type              = EMRobinBC
    variable          = E_real
    boundary          = right
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_right]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = right
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Top boundary (y = +4) ---
  [abc_real_top]
    type              = EMRobinBC
    variable          = E_real
    boundary          = top
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_top]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = top
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Bottom boundary (y = -4) ---
  [abc_real_bottom]
    type              = EMRobinBC
    variable          = E_real
    boundary          = bottom
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_bottom]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = bottom
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []
[]

[AuxVariables]
  # |E_z|² = E_real² + E_imag² — time-averaged field intensity.
  # This is proportional to the cycle-averaged Poynting flux magnitude.
  # Visualising E_intensity in ParaView reveals:
  #   - Peak near the source at (-2, 0): the driving current location
  #   - Secondary peak at (0, 0): internal focus inside the LHM slab
  #   - Third peak at (+2, 0): external (image) focus
  # The presence of these three local maxima along y = 0 at
  # x = {-2, 0, +2} is the signature of Veselago flat-lens focusing.
  # Their ratio (focus peak / source peak) measures the focusing quality.
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E_z|² nodewise from the primary field components.
  # ParsedAux evaluates the expression at each mesh node post-solve.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Point values along the optical axis y = 0.
  #
  # The flat-lens prediction: for source at (−2, 0), slab from −1 to +1,
  # the field intensity should peak at (0, 0) [internal focus] and at
  # (+2, 0) [external focus]. The ratio of these peaks to the source
  # amplitude quantifies the focusing efficiency degraded by the loss δ.
  #
  # Additional probes map the full x-axis and verify that:
  #   - Far left (x = -5): field decays as 1/r from the source
  #   - Left vacuum (x = -1.5): standing wave between source and slab
  #   - Slab interior: complex refracted field with sign-reversed phase
  #   - Right vacuum: refocused diverging wave from the image point
  #   - Far right (x = +5): field decays as 1/r from the image
  # ==================================================================

  # Source location (-2, 0) — primary field peak
  [E_real_source]
    type     = PointValue
    variable = E_real
    point    = '-2 0 0'
  []
  [E_imag_source]
    type     = PointValue
    variable = E_imag
    point    = '-2 0 0'
  []
  [E_intensity_source]
    type     = PointValue
    variable = E_intensity
    point    = '-2 0 0'
  []

  # Internal focus (0, 0) — expected intensity peak inside slab
  # For ideal n = -1 (zero loss): intensity here equals source intensity.
  # With δ = 0.3 loss, the peak is reduced but should be clearly visible.
  [E_real_int_focus]
    type     = PointValue
    variable = E_real
    point    = '0 0 0'
  []
  [E_imag_int_focus]
    type     = PointValue
    variable = E_imag
    point    = '0 0 0'
  []
  [E_intensity_int_focus]
    type     = PointValue
    variable = E_intensity
    point    = '0 0 0'
  []

  # External (image) focus (+2, 0) — the flat-lens image point
  # Veselago (1968): a point source at distance s = 1 from the slab
  # face is imaged at distance (d − s) = 1 on the far side.
  # The image quality (intensity ratio to source) quantifies the
  # regularised Pendry perfect-lens performance for propagating waves.
  [E_real_ext_focus]
    type     = PointValue
    variable = E_real
    point    = '2 0 0'
  []
  [E_imag_ext_focus]
    type     = PointValue
    variable = E_imag
    point    = '2 0 0'
  []
  [E_intensity_ext_focus]
    type     = PointValue
    variable = E_intensity
    point    = '2 0 0'
  []

  # Far field check (+4, 0) — should be weaker than the focus
  # If focusing is working, the field at (+4, 0) diverges from (+2, 0)
  # and intensity should drop below the focus peak, confirming a genuine
  # local maximum at the image point rather than a monotonic decay.
  [E_intensity_far]
    type     = PointValue
    variable = E_intensity
    point    = '4 0 0'
  []

  # Off-axis probe at (0, 1) — inside slab, off the optical axis.
  # Should be weaker than the on-axis internal focus, confirming that
  # the focusing is directional (converging toward y = 0) not uniform.
  [E_intensity_off_axis]
    type     = PointValue
    variable = E_intensity
    point    = '0 1 0'
  []

  # Left vacuum mid-point (-1.5, 0) — between source and slab face.
  # Expect superposition of outgoing wave from source and reflection
  # from the slab. With perfect impedance match the reflection should
  # be near zero (pure transmission), so field here is primarily
  # the outgoing cylindrical wave from the source.
  [E_intensity_left_vac]
    type     = PointValue
    variable = E_intensity
    point    = '-1.5 0 0'
  []

  # Right vacuum mid-point (+1.5, 0) — between slab face and focus.
  # Expect converging cylindrical wave toward (+2, 0); intensity
  # increases as we approach the focus from the left.
  [E_intensity_right_vac]
    type     = PointValue
    variable = E_intensity
    point    = '1.5 0 0'
  []

  # Total field energy in the domain: ∫ |E_z|² dA.
  # Provides a global measure of the field level; useful for
  # monitoring convergence and comparing lossless/lossy cases.
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner — required for this coupled system.
  #
  # Reasons full = true is necessary:
  #   (1) EMRobinBC at all 4 boundaries: introduces off-diagonal Jacobian
  #       blocks between E_real and E_imag through the j k₀ coupling term.
  #   (2) ADMatCoupledForce kernels in the LHM slab (|x| < 1): add
  #       additional off-diagonal blocks throughout the slab volume.
  #   (3) Negative diffusivity (inv_mu_r = -1 in slab): the system
  #       matrix is indefinite (not positive-definite), so iterative
  #       solvers without careful preconditioning would diverge.
  #
  # Without full=true, the SMP only builds diagonal blocks (one for
  # E_real, one for E_imag) and misses the strong cross-coupling in the
  # LHM slab and at the boundaries. Newton iterations diverge.
  #
  # The direct LU solver (MUMPS/SuperLU via PETSc) handles the
  # indefinite system exactly and is fast at this problem size:
  # 120 × 80 mesh × 2 variables = 19,282 DOFs (FIRST-order nodes).
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 2D coupled Helmholtz system with ~19k DOFs
  # is solved exactly by LU in a few seconds. Direct methods are preferable
  # here because:
  #   (1) The Helmholtz operator is indefinite (k₀² term changes sign of
  #       the otherwise positive-definite stiffness matrix).
  #   (2) In the LHM slab the diffusion coefficient is negative, making the
  #       local stiffness sub-matrix negative-definite — global spectrum
  #       has mixed signs that confuse GMRES or conjugate-gradient methods.
  #   (3) The off-diagonal blocks from the loss cross-coupling and Robin BCs
  #       make the system non-symmetric, ruling out standard CG.
  #
  # LU gives the exact solution in a single iteration of Newton (the system
  # is linear — E_z appears at most linearly in all kernels and BCs).
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output — E_real, E_imag, E_intensity on the full mesh.
  # Visualise in ParaView:
  #   - E_intensity: reveals the focusing pattern as a series of bright spots
  #     along the optical axis at x = {-2, 0, +2}
  #   - E_real, E_imag: show the phase structure of the refracted field;
  #     inside the LHM slab (−1 < x < 1) the phase progression is reversed
  #     compared to the vacuum regions, demonstrating backward phase velocity
  #   - Use Threshold on E_intensity to isolate the focus peaks
  #   - The LHM slab boundaries at x = ±1 are visible as sharp changes
  #     in the field gradient, reflecting the 1/μ_r interface condition
  exodus = true

  # csv: postprocessor values — intensity at source, internal focus, and
  # external focus, plus far-field and off-axis probes and total energy.
  # The key result: E_intensity_int_focus / E_intensity_source and
  # E_intensity_ext_focus / E_intensity_source should both be O(1) for
  # good focusing (they degrade toward zero as loss δ increases).
  csv    = true
[]
