# ============================================================
# Case 80: Multilayer Dielectric Stack — Bragg Mirror
# Saleh & Teich, Fundamentals of Photonics (Wiley, 2019), Ch. 7
# Yariv & Yeh, Photonics: Optical Electronics in Modern
#   Communications (Oxford, 2007), Ch. 6
# Born & Wolf, Principles of Optics, 7th Ed., Ch. 1
#
# A 1D distributed Bragg reflector (DBR) consisting of 5 periods
# of alternating high/low permittivity layers, each cut to quarter-
# wave optical thickness at the design frequency f₀ = 20 MHz.
# A PEC wall sits at x = 0. The periodic stack extends from x = 0
# to x = L_stack = 24.684 m. A vacuum half-space from x = L_stack
# to x = 75 m separates the stack from the Robin port.
#
# Physical layout:
#
#  x=0   1.875 4.937 6.812 9.874 11.749 14.811 16.686 19.747 21.622 24.684   x=75
#  |H1---|L1---|H2---|L2---|H3---|L3----|H4----|L4----|H5----|L5----|  vacuum  |
#  PEC                                                             stack/vac  port
#  (E=0)                   eps_r = 4 / 1.5 / 4 / ...             eps_r = 1
#
# Layer design — quarter-wave optical thickness:
# --------------------------------------------------
# At the design wavelength λ₀ = c/f₀ = 15 m, each layer is chosen so
# that its optical thickness equals one quarter of λ₀:
#
#   n_H × d_H = λ₀/4     →   d_H = 15/(4×2)       = 1.875  m
#   n_L × d_L = λ₀/4     →   d_L = 15/(4×√1.5)    ≈ 3.062  m
#
# where n_H = √(ε_H) = √4 = 2.0 and n_L = √(ε_L) = √1.5 ≈ 1.2247.
#
# Why quarter-wave thickness creates a stopband:
# -----------------------------------------------
# Each H→L (or L→H) interface produces a partial Fresnel reflection.
# For a quarter-wave layer, the round-trip phase through the layer is
# exactly π, so reflections from successive interfaces arrive back at
# the input face with a phase shift that makes them interfere
# CONSTRUCTIVELY.  This cumulative constructive interference builds up
# a very high reflectivity in the "stopband" — a bandwidth centred on
# f₀ where the mirror is highly reflective.
#
# Transfer-matrix analogy:
# -------------------------
# Each layer is described by its 2×2 transfer matrix M_j (Abeles method):
#
#   M_j = [ cos(φ_j)        −j sin(φ_j)/n_j ]
#         [ −j n_j sin(φ_j)  cos(φ_j)        ]
#
# where φ_j = (2π/λ₀) n_j d_j is the single-pass phase.  At design:
# φ_H = φ_L = π/2, so cos(π/2) = 0 and sin(π/2) = 1, giving:
#
#   M_H = [ 0        −j/n_H ]    M_L = [ 0        −j/n_L ]
#         [ −j n_H   0      ]          [ −j n_L   0      ]
#
# The combined transfer matrix for one period is:
#   M_period = M_L × M_H = [−n_L/n_H  0       ]
#                           [0         −n_H/n_L]
#
# After N periods and a PEC termination (where E = 0), the amplitude
# reflection coefficient is |R| = 1 — perfect reflection at the design
# wavelength for any number of periods N ≥ 1.  This means the Bragg
# mirror on a PEC wall has |R|² = 1 exactly at f₀.
#
# MOOSE verification:
# --------------------
# The finite-element solution of
#   d²E/dx² + k₀² ε_r(x) E = 0
# with PEC at x = 0 and Robin port at x = 75 should produce
# |R|² ≈ 1.000 (up to discretisation error) at the design frequency.
# This provides a quantitative check of the EM solver against the
# analytic transfer-matrix result.
#
# Governing equation (frequency domain, 1D, normal incidence):
#
#   d²E/dx² + k₀² ε_r(x) E = 0
#
# where k₀ = 2πf₀/c = 0.41888 rad/m and ε_r(x) is piecewise constant:
#
#   ε_r = 4   in H layers (n_H = 2)
#   ε_r = 1.5 in L layers (n_L = √1.5 ≈ 1.2247)
#   ε_r = 1   in the vacuum gap (x > 24.684 m) and at the port
#
# Real/imaginary splitting (lossless case):
# ------------------------------------------
# Because ε_r is purely real, the complex Helmholtz equation splits into
# two decoupled real equations — one for the real part and one for the
# imaginary part of the phasor:
#
#   d²E_real/dx² + k₀² ε_r(x) E_real = 0
#   d²E_imag/dx² + k₀² ε_r(x) E_imag = 0
#
# Both equations share the same spatially-varying coefficient k₀² ε_r(x).
# The two components couple only through the Robin port BC (see [BCs] block).
#
# Robin port boundary condition at x = L:
# -----------------------------------------
# The EMRobinBC injects the incident wave and absorbs the outgoing
# (reflected) wave without non-physical re-reflections:
#
#   ∂E/∂x + j k₀ E = 2 j k₀ E0 exp(j k₀ x)   at x = L
#
# The ReflectionCoefficient postprocessor then extracts |R|² from the
# total complex field at the port.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# k = 2π × 20 MHz / c = 2π × 20e6 / 3e8 = 0.41888 rad/m
# Same k as case 32 (dielectric slab) and case 33 (coupled resonators).
# Changing 'k' shifts the operating frequency; the layer thicknesses
# remain fixed so the quarter-wave condition is met only at k₀.
# Changing 'L' moves the Robin port; keep ≥ 2 free-space wavelengths
# of vacuum between the stack (at x=24.684 m) and the port.
k     = 0.41887902047863912   # free-space wavenumber [rad/m], f₀ = 20 MHz
L     = 75                    # domain length [m] = 5 × λ₀ = 5 × 15 m
E0    = 1                     # incident field amplitude [V/m]
theta = 0                     # incidence angle [degrees] — normal incidence

[Mesh]
  # 1D domain from x = 0 (PEC wall) to x = L (vacuum port).
  #
  # Resolution requirements:
  #   - Free-space wavelength:   λ₀     = 15.0  m  (10 el/λ → dx ≈ 1.5 m)
  #   - H-layer wavelength:      λ_H    = 7.5   m  (d_H/dx → ~18 el per H layer)
  #   - L-layer wavelength:      λ_L    ≈ 12.25 m  (d_L/dx → ~30 el per L layer)
  #   - H-layer thickness:       d_H    = 1.875 m
  #
  # nx = 750 gives dx = 0.1 m, which provides:
  #   • ~18.75 elements per H layer (thinnest layer, 1.875 m)
  #   • ~30.6  elements per L layer (3.062 m)
  #   • ~75 elements per free-space wavelength
  # This is more than adequate for first-order Lagrange elements on
  # this smooth piecewise-constant coefficient problem.
  [bragg]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 750
    xmin = 0
    xmax = ${L}
  []

  # Rename the auto-generated boundary labels to physically meaningful names.
  # 'metal'  → x = 0: Perfect Electric Conductor (PEC), E = 0
  # 'port'   → x = L: open port where the incident wave is injected
  [rename]
    type         = RenameBoundaryGenerator
    input        = bragg
    old_boundary = 'left right'
    new_boundary = 'metal port'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E(x) = E_real + j E_imag.
  # Both components use first-order Lagrange elements on the 1D mesh.
  # Lagrange elements enforce C⁰ continuity across every ε_r interface,
  # which is correct: the tangential E field must be continuous at each
  # H/L and stack/vacuum interface.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Decoupled from E_real in the bulk (lossless) but coupled through the
  # Robin port BC at x = L via the complex impedance term.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # coeff_real_fn: the spatially-varying reaction coefficient k₀² ε_r(x)
  #
  # The Helmholtz equation  d²E/dx² + k₀² ε_r(x) E = 0  is cast in
  # weak form using two MOOSE kernels:
  #
  #   Diffusion:    +∫ (dE/dx)(dv/dx) dx   →  strong form  −d²E/dx²
  #   ADMatReaction: −rate × ∫ E v dx       →  strong form  −rate × E
  #
  # For the residual to equal the Helmholtz equation (= 0):
  #   −d²E/dx² − rate × E = 0   ⟺   d²E/dx² + rate × E = 0
  # the rate must be POSITIVE: rate = +k₀² ε_r(x).
  #
  # Numerical values of k₀² ε_r for each region:
  #   H layers (ε_r = 4):   k₀² × 4   = 0.70184 rad²/m²
  #   L layers (ε_r = 1.5): k₀² × 1.5 = 0.26319 rad²/m²
  #   Vacuum   (ε_r = 1):   k₀² × 1   = 0.17546 rad²/m²
  #
  # Layer boundaries (precise to 6 decimal places):
  #   Period 1:  H [0.000000, 1.875000)  L [1.875000,  4.936862)
  #   Period 2:  H [4.936862, 6.811862)  L [6.811862,  9.873724)
  #   Period 3:  H [9.873724, 11.748724) L [11.748724, 14.810587)
  #   Period 4:  H [14.810587, 16.685587) L [16.685587, 19.747449)
  #   Period 5:  H [19.747449, 21.622449) L [21.622449, 24.684311)
  #   Vacuum:    [24.684311, 75.0]
  #
  # The nested if() expression encodes the 10-layer piecewise-constant
  # profile in a single ParsedFunction — no additional materials or
  # subdomain splitting is needed.
  # ------------------------------------------------------------------
  [coeff_real_fn]
    type       = ParsedFunction
    expression = 'if(x<1.875000, 0.701839, if(x<4.936862, 0.263189, if(x<6.811862, 0.701839, if(x<9.873724, 0.263189, if(x<11.748724, 0.701839, if(x<14.810587, 0.263189, if(x<16.685587, 0.701839, if(x<19.747449, 0.263189, if(x<21.622449, 0.701839, if(x<24.684311, 0.263189, 0.175460))))))))))'
  []

  # cosTheta = cos(0°) = 1.0 for normal incidence.
  # The EMRobinBC multiplies k₀ by this factor to accommodate oblique incidence:
  #   k_eff = k₀ cos(θ)   (x-component of the wavevector)
  # At normal incidence (θ = 0) we have cos(θ) = 1 and k_eff = k₀.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap coeff_real_fn as an AD material property so ADMatReaction can consume it.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each quadrature point
  # (using the quadrature point x-coordinate) and returns an ADReal, enabling
  # automatic differentiation through the spatially varying coefficient.
  # The same material property is shared by both the E_real and E_imag kernels,
  # which is correct since both equations have the same ε_r profile.
  [coeff_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_real_material'
    prop_values = 'coeff_real_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Equation for E_real:
  #   d²E_real/dx² + k₀² ε_r(x) E_real = 0
  #
  # Weak form (multiply by test function v, integrate by parts):
  #   ∫ (dE_real/dx)(dv/dx) dx − ∫ k₀² ε_r E_real v dx = 0
  # ==================================================================

  # Diffusion kernel: adds +∫ ∇E_real · ∇v dV to the residual.
  # In 1D this corresponds to the strong-form operator −d²E_real/dx².
  [diffusion_real]
    type     = Diffusion
    variable = E_real
  []

  # ADMatReaction: adds −rate × ∫ E_real v dx to the residual.
  # Combined with Diffusion:
  #   strong form = −d²E_real/dx² − k₀² ε_r(x) E_real = 0  ✓
  # Note the NEGATIVE sign convention in ADMatReaction: the kernel
  # adds (−rate × E) to the residual, not (+rate × E).
  [reaction_real]
    type          = ADMatReaction
    reaction_rate = coeff_real_material
    variable      = E_real
  []

  # NOTE: No cross-coupling term for E_real.
  # In a lossy medium with complex ε_r = ε_r' − j ε_r'', there would be
  # a term +k₀² ε_r''(x) E_imag coupling the two equations (implemented
  # via ADMatCoupledForce in the reference lossy benchmark). Since ε_r is
  # purely real here (ε_r'' = 0), that term vanishes and E_real decouples
  # from E_imag in the bulk.

  # ==================================================================
  # Equation for E_imag:
  #   d²E_imag/dx² + k₀² ε_r(x) E_imag = 0
  #
  # Identical structure to the E_real equation — the same ε_r profile
  # appears in both equations (lossless: no imaginary part of ε_r).
  # ==================================================================

  [diffusion_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Same reaction coefficient as for E_real: the Bragg stack is lossless
  # so no cross-term between E_real and E_imag appears in the bulk.
  [reaction_imag]
    type          = ADMatReaction
    reaction_rate = coeff_real_material
    variable      = E_imag
  []

  # NOTE: No cross-coupling term for E_imag either.
  # The lossy benchmark would include −k₀² ε_r''(x) E_real here.
  # Dropped for this lossless Bragg mirror.
[]

[BCs]
  # ------------------------------------------------------------------
  # PEC wall at x = 0 (boundary 'metal')
  # ------------------------------------------------------------------
  # A Perfect Electric Conductor forces the tangential E field to zero.
  # In phasor notation: E(0) = 0 for all time, which means:
  #   E_real(0) = 0   AND   E_imag(0) = 0
  # Both homogeneous Dirichlet BCs must be applied simultaneously.
  # This BC is the physical backing for the Bragg mirror; the PEC
  # reflection produces the phase reversal that, combined with the
  # quarter-wave layers, leads to constructive interference.
  [pec_real]
    type     = DirichletBC
    variable = E_real
    boundary = metal
    value    = 0
  []

  [pec_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = metal
    value    = 0
  []

  # ------------------------------------------------------------------
  # Robin port condition at x = L (boundary 'port')
  # ------------------------------------------------------------------
  # The EMRobinBC implements the first-order Sommerfeld-style port BC
  # (Jin, "FEM in Electromagnetics", 3rd Ed., Eq. 9.60):
  #
  #   ∂E/∂x + j k₀ cos(θ) E = 2 j k₀ cos(θ) E0 exp(j k₀ x)
  #
  # This boundary condition simultaneously:
  #   (1) Injects the incident plane wave E_inc = E0 exp(+j k₀ x)
  #   (2) Absorbs the outgoing reflected wave without spurious reflections
  #
  # At normal incidence the incident wave at x = L has:
  #   E_inc(L) = E0 exp(j k₀ L) = exp(j × 0.41888 × 75)
  #
  # The coupling between E_real and E_imag through the Robin BC arises
  # from the factor j = √(−1):
  #   The term j k₀ E = j k₀ (E_real + j E_imag)
  #                   = −k₀ E_imag + j k₀ E_real
  # so the real part of the BC involves E_imag and vice versa. This
  # off-diagonal coupling is the ONLY interaction between the two field
  # components in the lossless Bragg mirror problem.
  #
  # The 'sign = negative' parameter matches the convention in the
  # upstream reference benchmark (slab_reflection.i) and ensures the
  # residual sign is consistent with the Diffusion kernel convention.
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
  # ReflectionCoefficient: power reflection coefficient |R|²
  # ------------------------------------------------------------------
  # This postprocessor reconstructs the complex reflection coefficient
  # from the total field at the port boundary x = L.
  #
  # The field at the port is the superposition of the incident and
  # reflected waves (in vacuum, k_vac = k₀):
  #   E_total(L) = E_inc(L) + R × E_refl(L)
  #              = E0 exp(+j k₀ L) + R × E0 exp(−j k₀ L)
  #
  # The postprocessor solves for R (complex) and returns |R|².
  #
  # Expected result — analytic transfer-matrix prediction:
  # -------------------------------------------------------
  # At the design frequency f₀ = 20 MHz, every layer has exactly
  # quarter-wave optical thickness. The combined period matrix for
  # one H-L pair is diagonal: M_period = diag(−n_L/n_H, −n_H/n_L).
  # After N = 5 periods backed by PEC (E = 0 BC), the transfer-matrix
  # calculation gives |R| = 1 EXACTLY — perfect reflection.
  #
  # The MOOSE FEM solution should return |R|² ≈ 1.000 (deviations
  # below ~10⁻⁴ are expected from mesh discretisation only).
  #
  # Contrast with case 32 (single slab, |R|² < 1 due to partial
  # transmission and Fabry-Pérot oscillations): the periodic quarter-
  # wave stack eliminates all transmission by design.
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
  # Field diagnostics at selected points
  # ------------------------------------------------------------------
  # These PointValue postprocessors sample the electric field at key
  # locations to verify the spatial structure of the solution:
  #
  #   x = 0    (PEC wall):     both E_real and E_imag should be ~ 0
  #   x = 12.4 (mid-stack):   standing wave amplitude inside the mirror
  #   x = 50.0 (in vacuum):   reduced field amplitude outside the mirror
  #   x = 75.0 (port):        total field at the port
  #
  # In the vacuum region (x > 24.684 m), with |R|² ≈ 1, the field is
  # a nearly perfect standing wave:
  #   E_total = E_inc + R × E_refl ≈ 2j sin(k₀ x) (for R = −1)
  # so E_real ≈ 0 and E_imag ≈ 2 sin(k₀ x) at the design frequency.
  # Inside the stack, the field is greatly attenuated (evanescent-like
  # in the stopband sense) — the energy is reflected, not transmitted.

  # At the PEC wall (Dirichlet BC enforces E = 0 here)
  [E_real_at_pec]
    type     = PointValue
    variable = E_real
    point    = '0 0 0'
  []

  [E_imag_at_pec]
    type     = PointValue
    variable = E_imag
    point    = '0 0 0'
  []

  # Near the middle of the stack (Period 3 H layer, x = 10.81 m)
  # This probes the standing-wave pattern inside the mirror structure.
  [E_real_mid_stack]
    type     = PointValue
    variable = E_real
    point    = '10.81 0 0'
  []

  [E_imag_mid_stack]
    type     = PointValue
    variable = E_imag
    point    = '10.81 0 0'
  []

  # In the vacuum gap (x = 50 m), ~1.7 free-space wavelengths beyond
  # the stack. Field amplitude here reflects how much energy escapes.
  # For |R|² = 1 the standing-wave amplitude should be ~2 E0 = 2 V/m.
  [E_real_in_vacuum]
    type     = PointValue
    variable = E_real
    point    = '50 0 0'
  []

  [E_imag_in_vacuum]
    type     = PointValue
    variable = E_imag
    point    = '50 0 0'
  []

  # At the port boundary (x = L = 75 m)
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
  # Full Single Matrix Preconditioner: required for the coupled
  # (E_real, E_imag) system. Although the two Helmholtz equations are
  # decoupled in the bulk, the Robin BC at x = L introduces off-diagonal
  # Jacobian blocks coupling E_real and E_imag. Without full=true, the
  # preconditioner ignores these off-diagonal blocks and convergence
  # degrades near the port boundary.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 1D Helmholtz system (1500 DOFs for 750
  # elements with two variables) is small enough for an exact direct solve.
  # LU is preferred over iterative solvers here because the Helmholtz
  # operator can be indefinite near resonance, which causes GMRES/CG to
  # stagnate. LU solves in one step with no iteration.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true   # Full 1D field (E_real, E_imag) along x — visualise in ParaView
  csv    = true   # reflection_coefficient and field samples at diagnostic points
[]
