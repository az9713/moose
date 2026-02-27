# ============================================================
# Case 79: Oblique Incidence — Snell's Law and Total Internal Reflection
# Griffiths, "Introduction to Electrodynamics", 4th Ed. (2017), Ch. 9
# Hecht, "Optics", 5th Ed. (2016), Ch. 4
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed. (2014), Ch. 9
#
# Lectures 3-4: Green's functions for layered media, TE/TM decomposition,
#               Fresnel coefficients
#
# PHYSICAL SCENARIO
# -----------------
# A TE-polarised plane wave propagates in medium 1 (εᵣ = 4, n₁ = 2) and
# strikes a planar dielectric interface at x = 0. Medium 2 on the right
# has εᵣ = 1 (n₂ = 1, free space). The angle of incidence is θ = 20°,
# which is BELOW the critical angle θ_c = arcsin(n₂/n₁) = arcsin(0.5) = 30°.
# The transmitted wave refracts according to Snell's law:
#
#   n₁ sin(θ₁) = n₂ sin(θ₂)
#   2 × sin(20°) = 1 × sin(θ₂)
#   θ₂ = arcsin(2 sin 20°) = arcsin(0.6840) = 43.16°
#
# When θ₁ > θ_c = 30° (e.g. θ₁ = 45°), total internal reflection (TIR)
# occurs: no propagating wave exists in medium 2. Instead, an evanescent
# wave penetrates into medium 2 with exponential amplitude decay:
#
#   E_evan ∝ exp(−κ x) exp(j k_y y)   where κ = k₀ √(n₁² sin²θ₁ − n₂²)
#
# For θ₁ = 45° (TIR case): κ = k₀ √(4×0.5 − 1) = k₀ = π/2 ≈ 1.5708 m⁻¹
#   penetration depth = 1/κ = 2/π ≈ 0.637 m ≈ λ₀/2π (sub-wavelength)
#
# GOVERNING EQUATION
# ------------------
# For TE polarisation (E field has only a z-component, E = E_z ẑ), the
# 2D scalar Helmholtz equation holds in each medium:
#
#   ∇²E_z + k₀² εᵣ(x) E_z = 0
#
# where εᵣ(x) is piecewise constant:
#   εᵣ = 4  for x < 0   (medium 1, n₁ = 2)
#   εᵣ = 1  for x ≥ 0   (medium 2, n₂ = 1)
#
# WAVELENGTH AND WAVENUMBER CHOICE
# ---------------------------------
# λ₀ = 4 m   (same as cases 77–78 for comparability)
# k₀ = 2π/λ₀ = π/2 = 1.5707963...  rad/m
# k₀² = π²/4 = 2.4674011...  m⁻²
#
# Wave vectors (θ₁ = 20°):
#   k_x = k₀ n₁ cos(θ₁) = π cos(20°) = 2.9521314...  rad/m
#   k_y = k₀ n₁ sin(θ₁) = π sin(20°) = 1.0744879...  rad/m
#
# SCATTERED-FIELD FORMULATION (Background Field Decomposition)
# ------------------------------------------------------------
# The oblique incident plane wave
#   E_inc(x,y) = E₀ exp(j(k_x·x + k_y·y))
# has a phase variation along y that cannot be straightforwardly injected
# via the standard normal-incidence EMRobinBC port.
#
# Instead, we decompose the total field as:
#   E_z = E_inc + E_scat
#
# and solve for the SCATTERED field E_scat. Substituting into the
# Helmholtz equation:
#
#   ∇²(E_inc + E_scat) + k₀² εᵣ (E_inc + E_scat) = 0
#
# In medium 1, E_inc satisfies its own Helmholtz equation exactly
# (∇²E_inc + k₀²×4×E_inc = 0), so:
#
#   ∇²E_scat + k₀² εᵣ E_scat = −(∇²E_inc + k₀² εᵣ E_inc)
#                              = −k₀²(εᵣ − 4) E_inc
#
# This source term is NONZERO only in medium 2 (x > 0, εᵣ = 1 ≠ 4):
#   source = −k₀²(1 − 4) E_inc = +3k₀² E_inc   (in x > 0)
#   source = 0                                    (in x < 0)
#
# MOOSE BodyForce kernel sign convention:
#   MOOSE residual: ∫∇E_scat·∇v − ∫k₀²εᵣ E_scat v − ∫f v = 0
#   Strong form solved: ∇²E_scat + k₀²εᵣ E_scat = −f
#   We want RHS = +3k₀² E_inc, so f = −3k₀² E_inc
#   BodyForce function = −3k₀² × cos/sin(k_x x + k_y y)  [NEGATIVE SIGN]
#
# After solving for E_scat, the total field is:
#   E_total = E_inc + E_scat
# reconstructed as an AuxVariable for visualization.
#
# BOUNDARY CONDITIONS ON E_SCAT
# -------------------------------
# The scattered field E_scat is a purely outgoing/evanescent wave:
#   - Reflected wave propagating in medium 1 (x < 0)
#   - Transmitted wave propagating in medium 2 (x > 0)  [or evanescent in TIR]
#
# The four domain boundaries use EMRobinBC as first-order absorbing BCs
# (Sommerfeld radiation condition) with profile_func_real = 0 (no injection).
# This absorbs the outgoing scattered field without spurious reflections.
#
# FRESNEL COEFFICIENTS (analytical verification)
# -----------------------------------------------
# For TE polarisation at θ₁ = 20°, the Fresnel amplitude coefficients are:
#
#   r_s = (n₁ cos θ₁ − n₂ cos θ₂) / (n₁ cos θ₁ + n₂ cos θ₂)
#       = (2×cos20° − 1×cos43.16°) / (2×cos20° + 1×cos43.16°)
#       = (1.8794 − 0.7294) / (1.8794 + 0.7294) = 0.4408
#
#   t_s = 2 n₁ cos θ₁ / (n₁ cos θ₁ + n₂ cos θ₂) = 1.4408
#
# Power reflectance R = |r_s|² = 0.1943
# Power transmittance T = (n₂ cos θ₂)/(n₁ cos θ₁) × |t_s|² = 0.8057
# R + T = 1.0000 (energy conservation ✓)
#
# The FEM solution is verified by checking that the total field at
# representative points in each medium matches the Fresnel prediction.
#
# DOMAIN AND MESH
# ----------------
# The 2D domain is [−5, 5] × [−5, 5] m (10×10 m = 2.5λ₀ × 2.5λ₀).
# The dielectric interface is at x = 0 (vertical line).
# Medium 1 (εᵣ = 4) occupies x < 0; medium 2 (εᵣ = 1) occupies x > 0.
#
# Mesh: 100×100 QUAD4 elements:
#   λ₁ = λ₀/n₁ = 4/2 = 2 m → 100 elements / 10 m = 10 elements/m
#   → 10 × λ₁ = 20 elements per wavelength in medium 1  ✓
#   λ₂ = λ₀/n₂ = 4/1 = 4 m → 10 elements/m = 40 elements/λ₂  ✓
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (all referenced with ${} below)
# -----------------------------------------------------------
# k0 = 2π/λ₀ = 2π/4 = π/2 = 1.5707963... rad/m  (free-space wavenumber)
# kx = k₀ n₁ cos(θ) = π cos(20°)                 (x-component of k in medium 1)
# ky = k₀ n₁ sin(θ) = π sin(20°)                 (y-component of k, same in both media)
# source_coeff = −3 k₀²                           (BodyForce amplitude for scattered-field formulation)
# E0           = incident field amplitude (= 1, used in source scaling)
# theta        = dummy angle for EMRobinBC cosTheta = cos(0) = 1 (normal-incidence ABC approx)
k0           = 1.5707963267948966   # free-space wavenumber k₀ = π/2 [rad/m]
kx           = 2.9521314340786907   # k₀ n₁ cos(20°) = π cos(20°) [rad/m]
ky           = 1.0744879697552057   # k₀ n₁ sin(20°) = π sin(20°) [rad/m]
source_coeff = -7.4022033008170185  # −3 k₀² [m⁻²] — BodyForce amplitude
E0           = 1                    # incident field amplitude [V/m]
theta        = 0                    # incidence angle [degrees] — used in cosTheta_fn

[Mesh]
  # 2D square domain [−5, 5] × [−5, 5] m.
  # The dielectric interface at x = 0 divides the domain into two halves.
  # No explicit sub-domain meshing is needed — the ParsedFunction for
  # k₀² εᵣ(x) handles the material discontinuity via if(x<0, ..., ...).
  # C⁰ Lagrange elements enforce continuity of E_z at the interface,
  # correctly implementing the physical boundary condition that tangential
  # E_z is continuous across the dielectric boundary.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 100
    ny   = 100
    xmin = -5
    xmax = 5
    ymin = -5
    ymax = 5
  []

  # Rename the four boundary labels to compass-direction names for clarity.
  # GeneratedMeshGenerator produces: left, right, bottom, top.
  [rename]
    type         = RenameBoundaryGenerator
    input        = domain
    old_boundary = 'left right bottom top'
    new_boundary = 'west east south north'
  []
[]

[Variables]
  # E_scat_real — real part of the complex scattered field phasor.
  # The total field is E_z = E_scat + E_inc, where E_inc is known analytically.
  # Solving for E_scat (rather than E_total) allows us to use purely
  # absorbing BCs on all domain boundaries: the scattered field is outgoing.
  [E_scat_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_scat_imag — imaginary part of the complex scattered field phasor.
  # Decoupled from E_scat_real in the bulk (lossless, real εᵣ). The two
  # components couple only through the EMRobinBC absorbing boundary conditions.
  [E_scat_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # coeff_fn: spatially varying Helmholtz coefficient k₀² εᵣ(x)
  #
  # The Helmholtz equation for E_scat is:
  #   ∇²E_scat + k₀² εᵣ(x) E_scat = source
  #
  # Discretised in MOOSE as Diffusion + ADMatReaction with reaction_rate
  # = +k₀² εᵣ (positive rate, since ADMatReaction residual = −rate × E × test,
  # so combined strong form = −∇²E − rate×E = 0 → ∇²E + rate×E = 0).
  #
  # Medium 1 (x < 0, εᵣ = 4): k₀² × 4 = (π/2)² × 4 = π² = 9.8696044 m⁻²
  # Medium 2 (x ≥ 0, εᵣ = 1): k₀² × 1 = (π/2)² × 1 = π²/4 = 2.4674011 m⁻²
  # ------------------------------------------------------------------
  [coeff_fn]
    type       = ParsedFunction
    # k₀²×4 = 9.8696044 in medium 1 (εᵣ=4); k₀²×1 = 2.4674011 in medium 2 (εᵣ=1)
    expression = 'if(x<0, 9.8696044010893580, 2.4674011002723395)'
  []

  # ------------------------------------------------------------------
  # source_real_fn: real part of the BodyForce source term
  #
  # From the scattered-field formulation:
  #   source = −3 k₀² E_inc = −3 k₀² exp(j(k_x x + k_y y))
  # Real part:
  #   source_real = −3 k₀² cos(k_x x + k_y y)   [only for x > 0]
  #
  # Note the NEGATIVE sign: the MOOSE BodyForce residual is −∫f v,
  # so the strong form equation becomes ∇²E + k₀²εᵣ E = −f.
  # We want the physical RHS to be +3k₀² E_inc, so f = −3k₀² E_inc.
  # The function value is therefore −3k₀² × cos(k_x x + k_y y).
  # ------------------------------------------------------------------
  [source_real_fn]
    type       = ParsedFunction
    # f = source_coeff × cos(kx x + ky y) × H(x)
    # = −3k₀² × cos(...) in medium 2; zero in medium 1
    expression = 'if(x>0, ${source_coeff} * cos(${kx}*x + ${ky}*y), 0)'
  []

  # ------------------------------------------------------------------
  # source_imag_fn: imaginary part of the BodyForce source term
  #
  # source_imag = −3 k₀² sin(k_x x + k_y y)   [only for x > 0]
  # Same sign reasoning as for source_real_fn.
  # ------------------------------------------------------------------
  [source_imag_fn]
    type       = ParsedFunction
    expression = 'if(x>0, ${source_coeff} * sin(${kx}*x + ${ky}*y), 0)'
  []

  # ------------------------------------------------------------------
  # cosTheta_fn: cos(θ) factor for EMRobinBC
  #
  # The EMRobinBC implements ∂E/∂n + jk₀cos(θ)E = 0 for purely absorbing
  # operation (profile_func_real = 0). For the scattered field formulation,
  # we use a normal-incidence approximation (cosTheta = 1.0) at all four
  # boundaries. This is the standard first-order Sommerfeld ABC, which is
  # exact for waves hitting the boundary at normal incidence.
  # The scattered wave hits the left/right boundaries nearly normally for
  # a domain large enough (5 m = 1.25λ₀ from the interface to each boundary).
  # theta = 0 → cos(0) = 1.0.
  # ------------------------------------------------------------------
  [cosTheta_fn]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []

  # ------------------------------------------------------------------
  # Einc_real_fn, Einc_imag_fn: analytical incident field
  #
  # E_inc(x,y) = E₀ exp(j(k_x x + k_y y))
  # Real part:  E₀ cos(k_x x + k_y y)
  # Imag part:  E₀ sin(k_x x + k_y y)
  #
  # Used by the AuxKernels to reconstruct the total field:
  #   E_total = E_scat + E_inc
  # and to generate the visualisation of the incident beam separately.
  # ------------------------------------------------------------------
  [Einc_real_fn]
    type       = ParsedFunction
    expression = '${E0} * cos(${kx}*x + ${ky}*y)'
  []

  [Einc_imag_fn]
    type       = ParsedFunction
    expression = '${E0} * sin(${kx}*x + ${ky}*y)'
  []
[]

[Materials]
  # Wrap coeff_fn as an AD material property so ADMatReaction can consume it.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each quadrature
  # point and returns an ADReal, enabling automatic differentiation through
  # the spatially varying permittivity. Essential for Newton to converge.
  [coeff_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_material'
    prop_values = 'coeff_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Scattered-field Helmholtz equation for E_scat_real:
  #   ∇²E_scat_real + k₀²εᵣ(x) E_scat_real = −3k₀² cos(k_x x + k_y y) H(x)
  #
  # MOOSE residual (weak form):
  #   ∫∇E_scat_real·∇v − k₀²εᵣ∫E_scat_real v − ∫f_real v = 0
  #
  # where Diffusion gives the first term, ADMatReaction gives the second,
  # and BodyForce gives the third (f_real = source_real_fn as defined above).
  # ==================================================================

  # Laplacian term for E_scat_real: contributes −∇²E_scat_real to the
  # strong-form residual, i.e. +∇²E_scat_real to the Helmholtz operator.
  [diffusion_real]
    type     = Diffusion
    variable = E_scat_real
  []

  # Reaction term: ADMatReaction adds −k₀²εᵣ × E_scat_real × test to the
  # residual. With reaction_rate = +k₀²εᵣ the strong-form contribution is
  # −k₀²εᵣ E_scat_real, completing the Helmholtz operator on the left side.
  # lossless (εᵣ'' = 0) → no coupling between real and imaginary parts here.
  [field_real]
    type          = ADMatReaction
    variable      = E_scat_real
    reaction_rate = coeff_material
  []

  # Source term: BodyForce adds −∫f v to the residual (residual = −∫f·test).
  # f = source_real_fn = −3k₀² cos(k_x x + k_y y) H(x)
  # → MOOSE strong form: ∇²E + k₀²εᵣ E = −f = +3k₀² cos(...) H(x)  ✓
  [source_real]
    type     = BodyForce
    variable = E_scat_real
    function = source_real_fn
  []

  # ==================================================================
  # Scattered-field Helmholtz equation for E_scat_imag:
  #   ∇²E_scat_imag + k₀²εᵣ(x) E_scat_imag = −3k₀² sin(k_x x + k_y y) H(x)
  #
  # Same operator as E_scat_real; source is the imaginary part of E_inc.
  # The two components are decoupled in the bulk (lossless medium),
  # coupling only through the absorbing boundary conditions.
  # ==================================================================

  [diffusion_imag]
    type     = Diffusion
    variable = E_scat_imag
  []

  [field_imag]
    type          = ADMatReaction
    variable      = E_scat_imag
    reaction_rate = coeff_material
  []

  [source_imag]
    type     = BodyForce
    variable = E_scat_imag
    function = source_imag_fn
  []
[]

[BCs]
  # ==================================================================
  # ALL FOUR BOUNDARIES: First-order absorbing boundary condition (ABC)
  #
  # For the scattered field formulation, we need purely absorbing BCs:
  # the scattered field (reflected in medium 1, transmitted/evanescent
  # in medium 2) radiates outward and should exit without reflection.
  #
  # EMRobinBC with profile_func_real = 0 implements:
  #   ∂E_scat/∂n + j k₀ E_scat = 0   (outgoing wave condition)
  #
  # This is the first-order Sommerfeld ABC, exact for normally incident
  # plane waves. The scattered field hits the left boundary (reflected
  # beam, nearly normal to x = −5 face) and the right boundary (transmitted
  # beam, hitting x = +5 at the refracted angle θ₂ ≈ 43° from normal).
  # For the top/bottom boundaries, the evanescent/scattered waves arrive
  # at larger angles, but the domain is wide enough (5m = 1.25λ₀ from
  # the beam centre to each lateral boundary) to keep residual reflections
  # small compared to the scattered field near the interface.
  #
  # sign = negative: consistent with cases 32, 74, 75, 77.
  # ==================================================================

  # ------ WEST boundary (x = −5): absorb reflected scattered wave ------
  [west_real]
    type              = EMRobinBC
    variable          = E_scat_real
    boundary          = west
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  [west_imag]
    type              = EMRobinBC
    variable          = E_scat_imag
    boundary          = west
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  # ------ EAST boundary (x = +5): absorb transmitted scattered wave ------
  # The refracted wave in medium 2 propagates at θ₂ ≈ 43° from normal.
  # The first-order ABC is not perfectly transparent at this angle,
  # but residual reflections are small for the domain size used.
  [east_real]
    type              = EMRobinBC
    variable          = E_scat_real
    boundary          = east
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  [east_imag]
    type              = EMRobinBC
    variable          = E_scat_imag
    boundary          = east
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  # ------ NORTH boundary (y = +5): absorb laterally radiating scattered wave ------
  [north_real]
    type              = EMRobinBC
    variable          = E_scat_real
    boundary          = north
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  [north_imag]
    type              = EMRobinBC
    variable          = E_scat_imag
    boundary          = north
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  # ------ SOUTH boundary (y = −5): absorb laterally radiating scattered wave ------
  [south_real]
    type              = EMRobinBC
    variable          = E_scat_real
    boundary          = south
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []

  [south_imag]
    type              = EMRobinBC
    variable          = E_scat_imag
    boundary          = south
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta_fn
    profile_func_real = 0
    field_real        = E_scat_real
    field_imaginary   = E_scat_imag
    sign              = negative
  []
[]

[AuxVariables]
  # E_total_real — real part of the total electric field phasor.
  # E_total = E_scat + E_inc
  # Computed from the scattered-field solution plus the known incident field.
  [E_total_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_total_imag — imaginary part of the total field phasor.
  [E_total_imag]
    order  = FIRST
    family = LAGRANGE
  []

  # E_total_sq — squared magnitude |E_total|² = E_r² + E_i².
  # Proportional to the time-averaged power density (intensity).
  # Visualising |E_total|² in ParaView shows:
  #   - Incident beam pattern with oblique phase fronts in medium 1 (x < 0)
  #   - Refracted (transmitted) beam in medium 2 (x > 0) at a steeper angle
  #   - Interface standing-wave pattern from superposition of incident and reflected waves
  #   - Phase fronts in medium 2 rotated to θ₂ ≈ 43° consistent with Snell's law
  [E_total_sq]
    order  = FIRST
    family = LAGRANGE
  []

  # E_scat_sq — scattered field intensity |E_scat|².
  # Reveals the reflected and transmitted (scattered) components separately
  # from the total field, useful for extracting Fresnel coefficients numerically.
  [E_scat_sq]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Total field real part: E_total_real = E_scat_real + E₀ cos(k_x x + k_y y)
  # The incident field E_inc is added back analytically at each node.
  [total_real_aux]
    type              = ParsedAux
    variable          = E_total_real
    expression        = 'E_scat_real + ${E0} * cos(${kx}*x + ${ky}*y)'
    coupled_variables = 'E_scat_real'
    use_xyzt          = true
  []

  # Total field imaginary part: E_total_imag = E_scat_imag + E₀ sin(k_x x + k_y y)
  [total_imag_aux]
    type              = ParsedAux
    variable          = E_total_imag
    expression        = 'E_scat_imag + ${E0} * sin(${kx}*x + ${ky}*y)'
    coupled_variables = 'E_scat_imag'
    use_xyzt          = true
  []

  # Total field intensity: |E_total|² = E_r² + E_i²
  [total_sq_aux]
    type              = ParsedAux
    variable          = E_total_sq
    expression        = 'E_total_real^2 + E_total_imag^2'
    coupled_variables = 'E_total_real E_total_imag'
  []

  # Scattered field intensity: |E_scat|²
  [scat_sq_aux]
    type              = ParsedAux
    variable          = E_scat_sq
    expression        = 'E_scat_real^2 + E_scat_imag^2'
    coupled_variables = 'E_scat_real E_scat_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Field values at diagnostic points for Snell's law verification
  #
  # The refracted angle θ₂ ≈ 43.16° means the transmitted beam's phase
  # fronts in medium 2 have slopes:
  #   Δy / Δx_for_one_phase_cycle = ky / kx2 = ky / (k₀ cos θ₂)
  #
  # Comparing field values at points along the interface (x=0) and
  # within each medium validates the scattered-field formulation.
  # ==================================================================

  # --- On the interface (x = 0) at y = 0 ---
  # The total field at the interface should match Fresnel coefficients.
  # The transmitted amplitude t_s = 1.4408, so |E_total|² at x=0⁺ (medium 2)
  # should be approximately t_s² ≈ 2.08 when probed in the beam axis.
  [E_total_real_interface]
    type     = PointValue
    variable = E_total_real
    point    = '0 0 0'
  []

  [E_total_imag_interface]
    type     = PointValue
    variable = E_total_imag
    point    = '0 0 0'
  []

  # --- Deep in medium 1 (x = −2, y = 0): incident + reflected superposition ---
  # At x = −2 m in medium 1, the total field is the superposition of the
  # incident wave (amplitude 1) and reflected wave (amplitude r_s ≈ 0.4408).
  # The standing-wave envelope amplitude ranges from 1−|r_s| to 1+|r_s|,
  # i.e. from 0.559 to 1.441. The actual phase depends on exact position.
  [E_total_real_med1]
    type     = PointValue
    variable = E_total_real
    point    = '-2 0 0'
  []

  [E_total_imag_med1]
    type     = PointValue
    variable = E_total_imag
    point    = '-2 0 0'
  []

  # --- In medium 2 (x = +2, y = 0): transmitted refracted wave ---
  # In medium 2, the total field is approximately the transmitted wave
  # (amplitude t_s ≈ 1.4408) propagating at angle θ₂ ≈ 43.16°.
  # At x = +2, y = 0, the phase is k_x2 × 2 + k_y2 × 0 ≈ 1.1458 × 2 = 2.29 rad.
  [E_total_real_med2]
    type     = PointValue
    variable = E_total_real
    point    = '2 0 0'
  []

  [E_total_imag_med2]
    type     = PointValue
    variable = E_total_imag
    point    = '2 0 0'
  []

  # --- Scattered field amplitude at x = −1 (reflected wave in medium 1) ---
  # The scattered field in medium 1 is the reflected plane wave:
  #   E_scat = r_s × exp(j(−k_x x + k_y y))   (wave propagating in −x + k_y y direction)
  # At x = −1, y = 0: |E_scat| = |r_s| = 0.4408.
  [E_scat_real_reflected]
    type     = PointValue
    variable = E_scat_real
    point    = '-1 0 0'
  []

  [E_scat_imag_reflected]
    type     = PointValue
    variable = E_scat_imag
    point    = '-1 0 0'
  []

  # --- Scattered field amplitude at x = +1 (transmitted wave in medium 2) ---
  # In medium 2, E_scat ≈ (t_s − 1) × E_inc_med2 + E_inc_med2
  # Actually: E_total = t_s E_trans, E_scat = E_total − E_inc
  # At x = +1, y = 0: E_inc = cos(k_x × 1) + j sin(k_x × 1)
  # |E_scat| is not simply t_s − 1 because E_inc and E_trans have different wavevectors.
  [E_scat_real_transmitted]
    type     = PointValue
    variable = E_scat_real
    point    = '1 0 0'
  []

  [E_scat_imag_transmitted]
    type     = PointValue
    variable = E_scat_imag
    point    = '1 0 0'
  []

  # --- Total field intensity at three x-positions along y = 0 ---
  # Used to construct an intensity profile across the interface.
  [intensity_m3_y0]
    type     = PointValue
    variable = E_total_sq
    point    = '-3 0 0'
  []

  [intensity_m1_y0]
    type     = PointValue
    variable = E_total_sq
    point    = '-1 0 0'
  []

  [intensity_p1_y0]
    type     = PointValue
    variable = E_total_sq
    point    = '1 0 0'
  []

  [intensity_p3_y0]
    type     = PointValue
    variable = E_total_sq
    point    = '3 0 0'
  []

  # --- Average intensity in each medium ---
  # ElementAverageValue of E_total_sq over the full domain.
  # For an oblique incident beam traversing the interface, the average
  # intensity should be of order unity (incident amplitude E₀ = 1).
  [avg_intensity]
    type     = ElementAverageValue
    variable = E_total_sq
  []

  # --- Phase-front check in medium 2: Snell's law verification ---
  # The refracted wavevector in medium 2 is (k_x2, k_y) with:
  #   k_y2 = k_y = 1.0744879...  (phase-matching at interface)
  #   k_x2 = k₀ n₂ cos θ₂ = 1.1458...
  # The refraction angle θ₂ = atan(k_y/k_x2) = atan(1.0745/1.1458) ≈ 43.16°
  # This is consistent with Snell's law n₁ sin θ₁ = n₂ sin θ₂.
  #
  # Probing at two y-offsets along x = 2.0 verifies the transverse phase:
  #   Δφ = k_y × Δy = 1.0745 × Δy
  # For Δy = 1 m: expected phase difference ≈ 1.0745 rad.
  [E_total_real_p2_y0]
    type     = PointValue
    variable = E_total_real
    point    = '2 0 0'
  []

  [E_total_real_p2_y1]
    type     = PointValue
    variable = E_total_real
    point    = '2 1 0'
  []

  [E_total_real_p2_ym1]
    type     = PointValue
    variable = E_total_real
    point    = '2 -1 0'
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner.
  # The EMRobinBC at the four domain boundaries introduces off-diagonal
  # Jacobian coupling between E_scat_real and E_scat_imag (through the
  # imaginary factor j in the ABC condition ∂E/∂n + jk₀E = 0).
  # Without full=true, the preconditioner ignores these off-diagonal blocks
  # and Newton convergence degrades significantly at the boundary DOFs.
  # Using full=true assembles the complete Jacobian for the coupled
  # (E_scat_real, E_scat_imag) system.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 2D system has 100×100×2 = 20,000 DOFs.
  # The Helmholtz operator is indefinite (k₀²εᵣ term makes it non-coercive),
  # which makes iterative Krylov solvers unreliable without specialised
  # preconditioners. Direct LU (PETSc PCLU) is robust for this problem size
  # and completes in seconds on a modern workstation.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # Tight tolerances: the scattered-field Helmholtz system is linear in
  # E_scat, so Newton converges in 1–2 iterations. Setting tight tolerances
  # ensures the LU factorisation is fully utilised.
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output (E_scat_real, E_scat_imag, E_total_real,
  #   E_total_imag, E_total_sq, E_scat_sq) on the full mesh.
  # Visualise in ParaView:
  #   - E_total_sq: intensity pattern showing incident beam (left) and
  #     refracted beam (right) with oblique phase fronts
  #   - E_total_real: real part of the total field, showing the standing
  #     wave in medium 1 (x < 0) and the single traveling wave in medium 2
  #   - The interface at x = 0 is visible as a change in field scale
  #     because n₁ = 2 → n₂ = 1 changes the wavelength abruptly
  #   - Plot along y = 0 (horizontal cut) to see the field amplitude jump
  #     at x = 0 from 1 (incident) to ~1.44 (transmitted), consistent with
  #     Fresnel transmission t_s = 1.4408 for TE at 20°
  exodus = true

  # csv: postprocessor values for Fresnel verification.
  # The intensity values at x = ±1, ±2, ±3 along y = 0 allow construction
  # of the 1D amplitude profile across the interface for comparison with
  # the Fresnel prediction:
  #   Medium 1: |E|² oscillates between (1−|r|)² and (1+|r|)² = [0.31, 2.07]
  #   Medium 2: |E|² ≈ t² |cos θ_correction|  (single wave, no standing wave)
  csv    = true
[]
