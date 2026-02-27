# ============================================================
# Case 75: Lossy Drude Slab — Skin Depth and Attenuation
# Griffiths, "Introduction to Electrodynamics", 4th Ed., Ch. 9
# Cheng, "Field and Wave Electromagnetics", 2nd Ed., Ch. 8
#
# A plane EM wave at frequency f = 20 MHz propagates in the +x
# direction and enters a thin metallic slab described by the
# Drude free-electron model. The complex permittivity causes
# exponential attenuation (skin depth effect) and cross-coupling
# between the real and imaginary parts of the electric field.
#
# Physical layout:
#   x = 0        x = 30        x = 60        x = 100
#   |             |<-- slab -->|              |
#   PEC          interface    interface      port (Robin BC)
#   (E=0)        εᵣ=1→Drude   Drude→εᵣ=1   E_inc injected
#
# ============================================================
# Drude Model Physics
# ============================================================
#
# The Drude free-electron model gives the complex permittivity:
#
#   εᵣ(ω) = 1 - ωₚ² / (ω² + j γ ω)
#
# where:
#   ω  = 2π × 20 MHz = 1.25664×10⁸ rad/s  (operating frequency)
#   ωₚ = 2ω = 2.51327×10⁸ rad/s            (plasma frequency)
#   γ  = 0.3ω = 3.76991×10⁷ s⁻¹           (electron collision frequency)
#
# Substituting and simplifying:
#   εᵣ = 1 - ωₚ² / (ω² + j 0.3 ω²)
#      = 1 - (2ω)² / (ω²(1 + j 0.3))
#      = 1 - 4 / (1 + j 0.3)
#
# Rationalising the denominator (multiply by (1 - j 0.3)/(1 - j 0.3)):
#   4 / (1 + j 0.3) = 4(1 - j 0.3) / (1 + 0.09) = 4(1 - j 0.3) / 1.09
#
# Therefore:
#   εᵣ = 1 - 4(1 - j 0.3)/1.09
#      = 1 - 4/1.09 + j × 4 × 0.3/1.09
#      = 1 - 3.66972 + j × 1.10092
#      = -2.66972 + j × 1.10092
#
# Physical interpretation:
#   Re(εᵣ) = -2.66972  →  negative (below plasma frequency, evanescent-like)
#   Im(εᵣ) = +1.10092  →  positive (loss, exp(+jωt) convention)
#
# The skin depth (characteristic attenuation length) is:
#   δ = 1 / Im(k_z)
# where k_z = k₀ sqrt(εᵣ) = k₀ sqrt(-2.66972 + j 1.10092).
# Numerically: k_z = 0.1383 + j 0.6983 m⁻¹
#   → δ = 1/0.6983 = 1.432 m
# The 30 m slab is ~21 skin depths thick: the field amplitude entering
# at x = 30 is attenuated by a factor of exp(-21) ≈ 7 × 10⁻¹⁰ by x = 60.
# Practically the field is completely extinguished inside the slab; the
# transmitted region (x > 60) sees only the evanescent tail.
#
# ============================================================
# Governing Equations (Frequency Domain, 1D, Normal Incidence)
# ============================================================
#
# Starting from the Helmholtz equation with complex permittivity:
#
#   d²E/dx² + k₀² εᵣ(x) E = 0
#
# with k₀ = 2πf/c = 0.41888 rad/m and εᵣ(x) piecewise:
#   εᵣ = 1              for 0 ≤ x ≤ 30  (vacuum)
#   εᵣ = εᵣ' + j εᵣ''  for 30 < x ≤ 60 (Drude slab)
#   εᵣ = 1              for 60 < x ≤ 100 (vacuum)
#
# Writing E = E_real + j E_imag and separating real/imaginary parts:
#
#   d²E_real/dx² + k₀² εᵣ'(x) E_real + k₀² εᵣ''(x) E_imag = 0   ...(R)
#   d²E_imag/dx² + k₀² εᵣ'(x) E_imag - k₀² εᵣ''(x) E_real = 0   ...(I)
#
# The cross-coupling terms (±k₀² εᵣ''·E_imag / ±k₀² εᵣ''·E_real)
# are the hallmark of the lossy problem: they vanish when εᵣ'' = 0
# (lossless dielectric, Case 32) but are essential here.
#
# ============================================================
# MOOSE Weak Form Implementation
# ============================================================
#
# Each equation uses three kernels:
#
#   Diffusion:
#     Residual: +∫ (dE/dx)(dv/dx) dx  →  strong form: −d²E/dx²
#
#   ADMatReaction(rate = k₀² εᵣ'(x)):
#     Residual: −rate × ∫ E v dx      →  strong form: −k₀² εᵣ' E
#     Combined: −d²E/dx² − k₀² εᵣ' E = 0  ✓
#
#   ADMatCoupledForce(v = coupled_var, mat_prop_coef = c):
#     Residual: −c × ∫ coupled_var × test dx  →  strong form: −c × coupled_var
#
#   So for equation (R), setting c = k₀² εᵣ'' gives:
#     strong form: −d²E_real/dx² − k₀² εᵣ' E_real − k₀² εᵣ'' E_imag = 0
#     ↔  d²E_real/dx² + k₀² εᵣ' E_real + k₀² εᵣ'' E_imag = 0  ✓
#
#   For equation (I), setting c = −k₀² εᵣ'' gives:
#     strong form: −d²E_imag/dx² − k₀² εᵣ' E_imag + k₀² εᵣ'' E_real = 0
#     ↔  d²E_imag/dx² + k₀² εᵣ' E_imag − k₀² εᵣ'' E_real = 0  ✓
#
# ============================================================
# Numerical Values (precomputed)
# ============================================================
#
# k₀ = 2π × 20×10⁶ / 3×10⁸ = 2π/15 = 0.41887902... rad/m
# k₀² = (2π/15)² = 4π²/225 = 0.175460 m⁻²
#
# In vacuum (x ≤ 30 and x > 60):
#   k₀² × εᵣ'  = 0.175460 × 1       = +0.175460
#   k₀² × εᵣ'' = 0.175460 × 0       =  0  (no cross-coupling)
#
# In Drude slab (30 < x ≤ 60):
#   k₀² × εᵣ'  = 0.175460 × (-2.66972) = -0.46843
#   k₀² × εᵣ'' = 0.175460 × (+1.10092) = +0.19317
#   k₀² × (-εᵣ'') = -0.19317
#
# The sign flip in the E_real coefficient (positive in vacuum, negative in slab)
# reflects the Drude plasma's sub-plasma-frequency behaviour: rather than
# propagating (positive k²), the field is evanescent (negative k²) in the slab.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# k = 2π × 20 MHz / c = 2π × 20e6 / 3e8 = 0.41888 rad/m
k     = 0.41887902047863906   # free-space wavenumber [rad/m]
L     = 100                   # domain length [m]
E0    = 1                     # incident field amplitude [V/m]
theta = 0                     # incidence angle [degrees] — normal incidence

[Mesh]
  # 1D domain from x = 0 (PEC wall) to x = L (vacuum port).
  # 500 elements give ~5 elements per free-space wavelength (λ₀ = 15 m)
  # and good resolution of the skin-depth variation inside the slab.
  # The Drude slab occupies x ∈ [30, 60], centred in the domain.
  [slab]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 500
    xmin = 0
    xmax = ${L}
  []

  # Rename boundary labels to physically meaningful names.
  # 'metal'  → x = 0: Perfect Electric Conductor (PEC), E = 0
  # 'port'   → x = L: open port where the incident wave is injected
  [rename]
    type         = RenameBoundaryGenerator
    input        = slab
    old_boundary = 'left right'
    new_boundary = 'metal port'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E(x) = E_real + j E_imag.
  # Both components use first-order Lagrange elements. Lagrange elements enforce
  # C⁰ continuity across the vacuum–slab interfaces at x = 30 and x = 60,
  # consistent with tangential E-field continuity at material boundaries.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Cross-coupled to E_real throughout the slab via the ±k₀² εᵣ'' terms.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # coeff_real_fn: the reaction coefficient k₀² Re(εᵣ(x))
  #
  # Used as the rate in ADMatReaction for BOTH E_real and E_imag.
  # Value: +0.175460 in vacuum, −0.46843 in the Drude slab.
  # The negative value in the slab (below plasma frequency) causes the
  # Helmholtz operator to become an elliptic damping operator — the
  # field solution is evanescent (exponentially decaying) rather than
  # oscillatory in that region.
  # ------------------------------------------------------------------
  [coeff_real_fn]
    type       = ParsedFunction
    expression = 'if(x>30 & x<=60, -0.46843, 0.17546)'
  []

  # ------------------------------------------------------------------
  # coeff_imag_fn: the cross-coupling coefficient +k₀² Im(εᵣ(x))
  #
  # Applied in the E_real equation via ADMatCoupledForce(v = E_imag).
  # Residual contribution: −c × E_imag × test
  # Strong form contribution: −c × E_imag
  # For equation (R): d²E_real/dx² + k₀²εᵣ'E_real + k₀²εᵣ''E_imag = 0
  #   → c = +k₀²εᵣ'' = +0.19317 in slab, 0 in vacuum  ✓
  # ------------------------------------------------------------------
  [coeff_imag_fn]
    type       = ParsedFunction
    expression = 'if(x>30 & x<=60, 0.19317, 0)'
  []

  # ------------------------------------------------------------------
  # coeff_neg_imag_fn: the cross-coupling coefficient −k₀² Im(εᵣ(x))
  #
  # Applied in the E_imag equation via ADMatCoupledForce(v = E_real).
  # Residual contribution: −c × E_real × test
  # Strong form contribution: −c × E_real
  # For equation (I): d²E_imag/dx² + k₀²εᵣ'E_imag − k₀²εᵣ''E_real = 0
  #   → c = −k₀²εᵣ'' = −0.19317 in slab, 0 in vacuum  ✓
  # ------------------------------------------------------------------
  [coeff_neg_imag_fn]
    type       = ParsedFunction
    expression = 'if(x>30 & x<=60, -0.19317, 0)'
  []

  # cosTheta = cos(0°) = 1 for normal incidence.
  # EMRobinBC uses k_eff = k₀ cos(θ); at normal incidence k_eff = k₀.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap the three ParsedFunctions as AD material properties so the
  # AD kernels (ADMatReaction, ADMatCoupledForce) can consume them.
  # ADGenericFunctionMaterial evaluates each function at every quadrature
  # point and returns an ADReal, enabling automatic differentiation through
  # the spatially varying Drude permittivity.

  # k₀² Re(εᵣ(x)) — reaction coefficient for both equations
  [coeff_real_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_real_material'
    prop_values = 'coeff_real_fn'
  []

  # +k₀² Im(εᵣ(x)) — cross-coupling from E_imag into E_real equation
  [coeff_imag_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_imag_material'
    prop_values = 'coeff_imag_fn'
  []

  # −k₀² Im(εᵣ(x)) — cross-coupling from E_real into E_imag equation
  [coeff_neg_imag_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_neg_imag_material'
    prop_values = 'coeff_neg_imag_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Equation (R) for E_real:
  #   d²E_real/dx² + k₀² εᵣ'(x) E_real + k₀² εᵣ''(x) E_imag = 0
  #
  # Three kernels: Diffusion + ADMatReaction(real part) +
  #                ADMatCoupledForce(imaginary cross-coupling)
  # ==================================================================

  # Laplacian term: contributes −d²E_real/dx² to the strong form.
  [diffusion_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction term: adds −k₀² εᵣ'(x) E_real to the strong form.
  # Combined with diffusion: −d²E_real/dx² − k₀² εᵣ' E_real = 0
  #   ↔  d²E_real/dx² + k₀² εᵣ' E_real = 0   (bulk of equation R)
  [reaction_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = coeff_real_material
  []

  # Cross-coupling term from E_imag into the E_real equation.
  # ADMatCoupledForce residual = −mat_prop_coef × v × test
  #
  # With coeff_imag_material = +k₀² εᵣ'' = +0.19317 in the slab:
  #   residual contribution = −(+0.19317) × E_imag × test
  #
  # The full residual for equation (R) is then:
  #   R = Diffusion + ADMatReaction + ADMatCoupledForce
  #     = [+(dE_real/dx)(dv/dx)] + [−k₀²εᵣ' × E_real × v] + [−k₀²εᵣ'' × E_imag × v]
  # Setting R = 0 and integrating by parts recovers the strong form:
  #   d²E_real/dx² + k₀²εᵣ' E_real + k₀²εᵣ'' E_imag = 0  ✓
  [cross_real_from_imag]
    type         = ADMatCoupledForce
    variable     = E_real
    v            = E_imag
    mat_prop_coef = coeff_imag_material
  []

  # ==================================================================
  # Equation (I) for E_imag:
  #   d²E_imag/dx² + k₀² εᵣ'(x) E_imag − k₀² εᵣ''(x) E_real = 0
  #
  # Three kernels: Diffusion + ADMatReaction(real part) +
  #                ADMatCoupledForce(negative imaginary cross-coupling)
  # ==================================================================

  # Laplacian term: contributes −d²E_imag/dx² to the strong form.
  [diffusion_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction term: adds −k₀² εᵣ'(x) E_imag to the strong form.
  [reaction_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = coeff_real_material
  []

  # Cross-coupling term: adds +k₀² εᵣ''(x) E_real to the strong form.
  # With coeff_neg_imag_material = −k₀² εᵣ'' = −0.19317 in slab:
  #   ADMatCoupledForce residual = −(−k₀²εᵣ'') × E_real × test
  # Strong form contribution: −(coeff_neg_imag_material) × E_real
  #                         = −(−k₀²εᵣ'') × E_real = +k₀²εᵣ'' E_real
  # Full equation (I): d²E_imag/dx² + k₀²εᵣ'E_imag − k₀²εᵣ''E_real = 0  ✓
  [cross_imag_from_real]
    type         = ADMatCoupledForce
    variable     = E_imag
    v            = E_real
    mat_prop_coef = coeff_neg_imag_material
  []
[]

[BCs]
  # ------------------------------------------------------------------
  # PEC wall at x = 0 (boundary 'metal')
  # ------------------------------------------------------------------
  # A Perfect Electric Conductor forces the tangential E to zero.
  # In phasor form: E(0) = 0 ↔ E_real(0) = 0 AND E_imag(0) = 0.
  # Physical reasoning: the PEC supplies whatever current is needed to
  # cancel the incident field at its surface. In contrast to the Case 32
  # geometry (slab at x ∈ [0,25] against PEC), here the PEC backing at
  # x = 0 creates a standing wave in the vacuum gap (0 to 30) and forces
  # the field to zero at the wall. The slab (30 to 60) then attenuates
  # the transmitted field before it reaches the second vacuum region.
  [metal_real]
    type     = DirichletBC
    variable = E_real
    boundary = metal
    value    = 0
  []

  [metal_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = metal
    value    = 0
  []

  # ------------------------------------------------------------------
  # Robin port condition at x = L = 100 (boundary 'port')
  # ------------------------------------------------------------------
  # The EMRobinBC implements the first-order Sommerfeld-style port
  # condition (Jin, "FEM in Electromagnetics", 3rd Ed., Eq. 9.60):
  #
  #   ∂E/∂x + j k₀ cos(θ) E = 2 j k₀ cos(θ) E₀ exp(j k₀ cos(θ) x)
  #
  # This simultaneously:
  #   (1) Injects the incident plane wave E_inc = E₀ exp(+j k₀ x)
  #   (2) Absorbs outgoing (reflected) waves without spurious reflections
  #
  # The port couples E_real ↔ E_imag through the imaginary unit j k₀:
  #   j k₀ E = j k₀ (E_real + j E_imag) = -k₀ E_imag + j k₀ E_real
  # This off-diagonal coupling in the port Jacobian is the reason full
  # SMP preconditioning is required.
  #
  # Parameters:
  #   coeff_real = k₀ = 0.41888 rad/m
  #   func_real  = cosTheta = 1.0 (normal incidence)
  #   profile_func_real = E₀ = 1 V/m (unit amplitude)
  #   sign = negative  (consistent with the INL benchmark slab_reflection.i)
  [port_real]
    type              = EMRobinBC
    variable          = E_real
    boundary          = port
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [port_imag]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = port
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []
[]

[Postprocessors]
  # ------------------------------------------------------------------
  # ReflectionCoefficient — power reflectance |R|² at the port
  # ------------------------------------------------------------------
  # Reconstructs the complex reflection coefficient from the total
  # field at x = L by subtracting the known incident wave:
  #   R_complex = (E_total(L) − E_inc(L)) / E_ref(L)
  #   |R|² = |R_complex|²
  # where E_inc(L) = E₀ exp(+j k₀ L) and E_ref(L) = E₀ exp(−j k₀ L).
  #
  # Because the Drude slab is ~21 skin depths thick, it is effectively
  # opaque: the transmitted field (x > 60) is negligible (< 10⁻⁹ of E₀).
  # The skin depth δ = 1/Im(k₀√εᵣ) = 1/0.6983 = 1.432 m, and the wave
  # is essentially fully reflected at the vacuum–slab interface at x = 60.
  # Expected |R|² ≈ 1.0 (near-total reflection from the thick Drude slab).
  # Any deviation from unity is due to ohmic absorption in the slab skin.
  [reflection_coefficient]
    type                     = ReflectionCoefficient
    k                        = ${k}
    theta                    = ${theta}
    length                   = ${L}
    incoming_field_magnitude = ${E0}
    field_real               = E_real
    field_imag               = E_imag
    boundary                 = port
  []

  # ------------------------------------------------------------------
  # Field values at key spatial locations (skin depth diagnostics)
  # ------------------------------------------------------------------

  # Field at x = 15: inside the vacuum gap between PEC and slab.
  # Expect a standing wave pattern from interference between the
  # incident field and the reflection at the vacuum–slab interface (x=30).
  [E_real_vacuum_gap]
    type     = PointValue
    variable = E_real
    point    = '15 0 0'
  []

  [E_imag_vacuum_gap]
    type     = PointValue
    variable = E_imag
    point    = '15 0 0'
  []

  # Field at x = 35: just inside the Drude slab (5 m from entrance).
  # Should show moderate attenuation relative to the slab entrance at x=30.
  [E_real_slab_near]
    type     = PointValue
    variable = E_real
    point    = '35 0 0'
  []

  [E_imag_slab_near]
    type     = PointValue
    variable = E_imag
    point    = '35 0 0'
  []

  # Field at x = 45: middle of the Drude slab.
  # Attenuation is more pronounced here; |E|² ≪ |E at x=35|² verifies
  # exponential decay consistent with the skin depth calculation.
  [E_real_slab_mid]
    type     = PointValue
    variable = E_real
    point    = '45 0 0'
  []

  [E_imag_slab_mid]
    type     = PointValue
    variable = E_imag
    point    = '45 0 0'
  []

  # Field at x = 55: deep inside the slab, 5 m from exit.
  # Should be nearly extinguished if slab thickness ~ several skin depths.
  [E_real_slab_deep]
    type     = PointValue
    variable = E_real
    point    = '55 0 0'
  []

  [E_imag_slab_deep]
    type     = PointValue
    variable = E_imag
    point    = '55 0 0'
  []

  # Field at x = 80: in the transmitted vacuum region beyond the slab.
  # The Drude slab is ~21 skin depths thick (δ = 1.43 m, slab = 30 m),
  # so the amplitude is attenuated by exp(-30/1.43) ≈ exp(-21) ≈ 10⁻⁹.
  # Numerically this will be near machine epsilon — effectively zero.
  # The incident wave is almost entirely reflected at x = 60.
  [E_real_transmitted]
    type     = PointValue
    variable = E_real
    point    = '80 0 0'
  []

  [E_imag_transmitted]
    type     = PointValue
    variable = E_imag
    point    = '80 0 0'
  []

  # Field at the port x = 100 for standing-wave diagnostics.
  [E_real_at_port]
    type     = PointValue
    variable = E_real
    point    = '${L} 0 0'
  []

  [E_imag_at_port]
    type     = PointValue
    variable = E_imag
    point    = '${L} 0 0'
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner: essential for the coupled
  # (E_real, E_imag) system because:
  #   (1) The Robin BC at the port introduces off-diagonal Jacobian
  #       blocks between E_real and E_imag through j k₀.
  #   (2) The ADMatCoupledForce kernels introduce additional off-diagonal
  #       blocks throughout the slab (the ±k₀² εᵣ'' cross-coupling).
  # Without full=true, the preconditioner ignores all these couplings
  # and Newton convergence degrades severely in the lossy slab region.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorization: the 1D coupled Helmholtz system (1000 DOFs
  # for 500 elements × 2 variables) is small enough that a direct solve
  # is fast and numerically exact. Iterative solvers can struggle with
  # near-singular Helmholtz operators, particularly below the plasma
  # frequency where the real part of k² is negative.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true   # Full 1D field (E_real, E_imag) along x — visualise in ParaView
  csv    = true   # Reflection coefficient and spatial field samples
[]
