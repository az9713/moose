# ============================================================
# Case 32: EM Wave Reflection from a Lossless Dielectric Slab
# Haus, Electromagnetic Noise and Quantum Optical Measurements
#   (Springer, 2000), Ch. 1
#
# A plane wave at frequency f = 20 MHz propagates in the +x
# direction and strikes a lossless dielectric slab (εᵣ = 4)
# backed by a Perfect Electric Conductor (PEC) wall at x = 0.
#
# Physical layout:
#   x = 0          x = 25          x = 75
#   |<---- slab --->|<---- vacuum --->|
#   PEC           interface        port (Robin BC)
#   (E=0)         εᵣ = 4→1        E_inc injected
#
# Governing equation (frequency domain, 1D, normal incidence):
#
#   d²E/dx² + k₀² εᵣ(x) E = 0
#
# where k₀ = 2πf/c = 0.41888 rad/m is the free-space wavenumber
# and εᵣ(x) is piecewise constant:
#   εᵣ = 4  for 0 ≤ x ≤ 25  (slab,   n = 2)
#   εᵣ = 1  for 25 < x ≤ 75 (vacuum, n = 1)
#
# Real/imaginary splitting (lossless case):
# -----------------------------------------
# Because εᵣ is purely real, the governing equation splits cleanly
# into two decoupled real equations:
#
#   d²E_real/dx² + k₀² εᵣ(x) E_real = 0
#   d²E_imag/dx² + k₀² εᵣ(x) E_imag = 0
#
# For a LOSSY slab (complex εᵣ = εᵣ' - jεᵣ''), there would be
# cross-coupling terms:
#   d²E_real/dx² + k₀² εᵣ'(x) E_real + k₀² εᵣ''(x) E_imag = 0
#   d²E_imag/dx² + k₀² εᵣ'(x) E_imag - k₀² εᵣ''(x) E_real = 0
# The ADMatCoupledForce kernels in the reference benchmark (slab_reflection.i)
# implement those cross-coupling terms via JinSlabCoeffFunc. Here we drop them
# entirely because εᵣ'' = 0 (lossless), which simplifies the system greatly.
#
# Robin port boundary condition at x = L:
# ----------------------------------------
# The EMRobinBC implements the first-order Sommerfeld-style port condition
# (Jin, "FEM in Electromagnetics", 3rd Ed., Eq. 9.60):
#
#   ∂E/∂x + j k₀ E = 2 j k₀ E_inc e^{j k₀ x}   at x = L
#
# Split into real and imaginary parts, this couples E_real and E_imag
# at the port boundary only. The incident wave amplitude is E0 = 1 V/m.
# The port simultaneously:
#   (1) injects the incident wave E_inc = E0 exp(j k₀ x)
#   (2) absorbs outgoing (reflected) waves without non-physical reflections
#
# Analytical result (Fresnel, lossless slab, normal incidence):
# --------------------------------------------------------------
# For a semi-infinite slab (no PEC backing), the Fresnel amplitude
# reflection coefficient is:
#   R_Fresnel = (n₁ - n₂)/(n₁ + n₂) = (1 - 2)/(1 + 2) = -1/3
#   |R_Fresnel| = 1/3 ≈ 0.333,   |R|² = 1/9 ≈ 0.111
#
# The PEC backing adds a secondary reflection, and the finite slab thickness
# (25 m ≈ 1.67λ_slab with λ_slab = λ₀/n = 15/2 = 7.5 m) introduces
# Fabry-Pérot interference between the two slab interfaces. The MOOSE
# ReflectionCoefficient postprocessor extracts |R|² from the computed
# complex field at the port boundary.
#
# Expected behaviour: |R|² output from ReflectionCoefficient should be
# close to the Fabry-Pérot value for this slab-on-PEC geometry. The
# pure Fresnel value (1/9 ≈ 0.111) applies only in the limit of an
# impedance-matched or non-reflecting backing.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# k = 2π × 20 MHz / c = 2π × 20e6 / 3e8 = 0.41888 rad/m
# Changing 'k' rescales the operating frequency proportionally.
# Changing 'L' varies the total domain length (keep ≥ 2–3 wavelengths
# from the slab interface at x=25 to the port to allow the field to
# resolve its spatial variation cleanly).
k     = 0.41887902047863906   # free-space wavenumber [rad/m]
L     = 75                    # domain length [m] = 5 × λ₀ = 5 × 15 m
E0    = 1                     # incident field amplitude [V/m]
theta = 0                     # incidence angle [degrees] — normal incidence

[Mesh]
  # 1D domain from x = 0 (PEC wall) to x = L (vacuum port).
  # 300 elements gives ~3 elements per free-space wavelength (λ₀=15m)
  # and ~6 per slab wavelength (λ_slab = 7.5m), sufficient for FIRST-order
  # Lagrange to resolve both the incident and standing-wave fields.
  [slab]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 300
    xmin = 0
    xmax = ${L}
  []

  # Rename the auto-generated boundary labels to physically meaningful names.
  # 'metal'  → x = 0: Perfect Electric Conductor (PEC), E = 0
  # 'vacuum' → x = L: open port where the incident wave is injected
  [rename]
    type         = RenameBoundaryGenerator
    input        = slab
    old_boundary = 'left right'
    new_boundary = 'metal vacuum'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E(x) = E_real + j E_imag.
  # Both components use first-order Lagrange elements on the 1D mesh.
  # Lagrange elements enforce C⁰ continuity across the εᵣ interface at x=25,
  # which is appropriate since the tangential E field must be continuous there.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Decoupled from E_real in the bulk (lossless) but coupled through the Robin BC.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # coeff_real_fn: the reaction coefficient  k₀² εᵣ(x)
  #
  # The Helmholtz equation  d²E/dx² + k₀² εᵣ E = 0  is written in
  # weak form using Diffusion + ADMatReaction.  Diffusion contributes
  # +∫ (dE/dx)(dv/dx) dx  →  strong form  −d²E/dx²
  # ADMatReaction(rate) contributes  −rate × ∫ E v dx  →  strong form  −rate×E
  #
  # For the two terms to sum to the Helmholtz residual:
  #   −d²E/dx² − rate×E = 0  ⟺  d²E/dx² + rate×E = 0
  # So rate must be POSITIVE: rate = +k₀² εᵣ.
  #
  # Numerical values:
  #   Slab   (0 ≤ x ≤ 25): +k₀² × 4 = +0.70184
  #   Vacuum (x > 25):      +k₀² × 1 = +0.17546
  # ------------------------------------------------------------------
  [coeff_real_fn]
    type       = ParsedFunction
    expression = 'if(x<=25, 0.70184, 0.17546)'
  []

  # cosTheta = cos(0°) = 1.0 for normal incidence.
  # EMRobinBC multiplies the wavenumber by this factor to handle oblique incidence:
  #   k_eff = k₀ cos(θ)   (component of k in the direction of domain traversal)
  # At normal incidence θ = 0 so cos(θ) = 1 and k_eff = k₀.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap coeff_real_fn as an AD material property so ADMatReaction can consume it.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each quadrature point
  # and returns an ADReal, enabling automatic differentiation through the coefficient.
  [coeff_real_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_real_material'
    prop_values = 'coeff_real_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Equation for E_real:
  #   d²E_real/dx² + k₀² εᵣ(x) E_real = 0
  # Weak form:
  #   ∫ (dE_real/dx)(dv/dx) dx + ∫ (−k₀² εᵣ) E_real v dx = 0
  # ==================================================================

  # Diffusion kernel: adds +∫ ∇E_real · ∇v dV to the residual.
  # In 1D and strong form this is −d²E_real/dx².
  [diffusion_real]
    type     = Diffusion
    variable = E_real
  []

  # ADMatReaction: adds coeff × ∫ E_real v dx to the residual.
  # coeff = −k₀² εᵣ(x) so the strong form contribution is −k₀² εᵣ E_real.
  # Combined with Diffusion: strong form = −d²E_real/dx² − k₀² εᵣ E_real = 0  ✓
  [field_real]
    type          = ADMatReaction
    reaction_rate = coeff_real_material
    variable      = E_real
  []

  # NOTE: No ADMatCoupledForce kernel for E_real here.
  # In the LOSSY benchmark (slab_reflection.i) there would be a term
  #   +k₀² εᵣ''(x) E_imag
  # from the imaginary part of the permittivity. Since εᵣ'' = 0 (lossless),
  # that cross-coupling vanishes and E_real equation is self-contained.

  # ==================================================================
  # Equation for E_imag:
  #   d²E_imag/dx² + k₀² εᵣ(x) E_imag = 0
  # Identical structure to E_real — both use the same real coefficient.
  # ==================================================================

  [diffusion_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Same real coefficient as for E_real: lossless → no cross-term.
  [field_imag]
    type          = ADMatReaction
    reaction_rate = coeff_real_material
    variable      = E_imag
  []

  # NOTE: No ADMatCoupledForce kernel for E_imag either.
  # The lossy benchmark would have a term −k₀² εᵣ''(x) E_real here.
  # Dropped for the lossless slab.
[]

[BCs]
  # ------------------------------------------------------------------
  # PEC wall at x = 0 (boundary 'metal')
  # ------------------------------------------------------------------
  # A Perfect Electric Conductor forces the tangential electric field to zero.
  # In phasor form: E(0) = 0  ⟹  E_real(0) = 0  AND  E_imag(0) = 0.
  # Homogeneous Dirichlet BCs on both field components implement this.
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
  # Robin port condition at x = L (boundary 'vacuum')
  # ------------------------------------------------------------------
  # The EMRobinBC implements the first-order absorbing/port condition
  # (Sommerfeld-style, scalar form for plane waves, Jin Eq. 9.60):
  #
  #   ∂E/∂x + j k₀ cos(θ) E = 2 j k₀ cos(θ) E0 exp(j k₀ cos(θ) x)
  #
  # In the source code (EMRobinBC.C) the complex kernel is:
  #   common = j × coeff_real × func_real   (with coeff_imag=0, func_imag=0)
  #   lhs = common × E     (field-dependent, couples real ↔ imag)
  #   rhs = 2 × common × profile_func × exp(common × x)   (incident wave)
  #   residual = sign × test × (rhs − lhs)
  #
  # For our parameters:
  #   coeff_real = k₀ = 0.41888
  #   func_real  = cosTheta = cos(0°) = 1.0
  #   profile_func_real = E0 = 1 (uniform amplitude, plane wave profile)
  #   sign = negative  (consistent with the upstream reference benchmark)
  #
  # The port condition couples E_real and E_imag through the imaginary unit j:
  #   common = j k₀ is purely imaginary, so:
  #     lhs_real = Re(j k₀ E) = −k₀ E_imag
  #     lhs_imag = Im(j k₀ E) = +k₀ E_real
  # This coupling (off-diagonal in the BC Jacobian) is the only place the
  # two field components interact in the lossless slab problem.
  [vacuum_real]
    type              = EMRobinBC
    variable          = E_real
    boundary          = vacuum
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [vacuum_imag]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = vacuum
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
  # ReflectionCoefficient computes the power reflection coefficient |R|²
  # from the total field at the port boundary.
  #
  # The postprocessor reconstructs the reflected amplitude by:
  #   E_total(L) = E_inc(L) + R × E_ref(L)
  #   R_complex = (E_total(L) − E_inc(L)) / E_ref(L)
  #   |R|² = |R_complex|²
  # where E_inc(L) = E0 exp(+j k₀ L) and E_ref(L) = E0 exp(−j k₀ L).
  #
  # For a lossless PEC-backed slab with εᵣ = 4 and thickness d = 25 m:
  #   n = sqrt(4) = 2,  λ_slab = λ₀/n = 15/2 = 7.5 m
  #   slab thickness in wavelengths: d/λ_slab = 25/7.5 ≈ 3.33 λ_slab
  # The Fabry-Pérot resonance conditions modify |R|² from the simple
  # single-interface Fresnel value of 1/9 ≈ 0.111.
  # The MOOSE-computed value is the exact result for this specific geometry.
  [reflection_coefficient]
    type                     = ReflectionCoefficient
    k                        = ${k}
    theta                    = ${theta}
    length                   = ${L}
    incoming_field_magnitude = ${E0}
    field_real               = E_real
    field_imag               = E_imag
    boundary                 = vacuum
  []

  # Field amplitude at the port for diagnostics.
  # The magnitude sqrt(E_real² + E_imag²) at x = L reveals whether the
  # standing-wave pattern inside the domain is correctly resolved.
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
  # Full Single Matrix Preconditioner: essential for the coupled (E_real, E_imag)
  # system because the Robin BC introduces off-diagonal Jacobian blocks between
  # the two field components. Without full=true the preconditioner ignores these
  # couplings and convergence degrades significantly near the port.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorization: the 1D Helmholtz system (600 DOFs for 300 elements
  # with two variables) is small enough that a direct solve is fast and exact.
  # Iterative solvers can struggle with Helmholtz near resonance; LU avoids that.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true   # Full 1D field (E_real, E_imag) along x — visualise in ParaView
  csv    = true   # reflection_coefficient and port field values
[]
