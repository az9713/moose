# ============================================================
# Case 57: Xenon-135 Poisoning Transient — Reactor Flux Decay
# Simulates the coupled evolution of neutron flux and Xe-135
# concentration following a reactor power reduction.
#
# Simplified 0-D point kinetics in a 1D slab geometry:
#   d(phi)/dt = (nu_Sf - Sigma_a - sigma_Xe * Xe) * v * phi - D * laplacian(phi)
#   d(Xe)/dt  = gamma_Xe * Sigma_f * phi + lambda_I * I - lambda_Xe * Xe - sigma_Xe * Xe * phi
#   d(I)/dt   = gamma_I * Sigma_f * phi - lambda_I * I
#
# Simplified model (no spatial diffusion for Xe/I, they're volumetric):
#   Neutron flux phi: diffusion + reaction + Xe poisoning
#   Iodine-135 I: production from fission, decay to Xe
#   Xenon-135 Xe: production from I decay + direct fission yield,
#                  destruction by decay + neutron absorption
#
# Parameters (scaled for educational clarity):
#   D = 1.0 cm (diffusion coefficient)
#   nu_Sf = 0.10 /cm, Sigma_a = 0.08 /cm (k_inf ~ 1.25 without Xe)
#   sigma_Xe = 2.0e6 barns = 2.0e-18 cm^2 (huge cross section)
#   gamma_I = 0.061 (I-135 fission yield)
#   gamma_Xe = 0.003 (direct Xe-135 fission yield)
#   lambda_I = 2.87e-5 /s (I-135 half-life ~ 6.7 hr)
#   lambda_Xe = 2.09e-5 /s (Xe-135 half-life ~ 9.2 hr)
#   Sigma_f = 0.05 /cm (fission cross section)
#
# To keep things tractable in a short educational run:
#   - Scale time in hours (3600 s = 1 hr)
#   - Use scaled rate constants
#   - Domain: [0, 20] cm slab, 40 elements
#   - Run for 24 hours to see Xe buildup and flux depression
#
# The classic "Xe pit" effect: after power reduction, Xe-135
# builds up from I-135 decay, causing further flux depression.
# ============================================================

# Time in hours for readable output
# lambda_I = 2.87e-5 /s * 3600 s/hr = 0.1033 /hr
# lambda_Xe = 2.09e-5 /s * 3600 s/hr = 0.0752 /hr

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 2
    xmin = 0
    xmax = 20.0
    ymin = 0
    ymax = 0.5
  []
[]

[Variables]
  [phi]
    # Neutron flux (arbitrary units, normalized to ~1 at center)
  []
  [Xe]
    # Xenon-135 concentration (atoms/cm^3, scaled)
  []
  [I135]
    # Iodine-135 concentration (atoms/cm^3, scaled)
  []
[]

[ICs]
  # Initial steady-state flux profile (cosine shape)
  [phi_ic]
    type = FunctionIC
    variable = phi
    function = 'cos(3.14159265*(x - 10)/20)'
  []
  # Initial Xe at steady-state equilibrium:
  # Xe_eq = (gamma_I + gamma_Xe) * Sigma_f * phi / (lambda_Xe + sigma_Xe_eff * phi)
  # With scaled values: ~ 0.5 (set by function)
  [Xe_ic]
    type = FunctionIC
    variable = Xe
    function = '0.5 * cos(3.14159265*(x - 10)/20)'
  []
  # Initial I135 at steady-state: I_eq = gamma_I * Sigma_f * phi / lambda_I
  [I_ic]
    type = FunctionIC
    variable = I135
    function = '0.3 * cos(3.14159265*(x - 10)/20)'
  []
[]

[Kernels]
  # === Neutron flux equation ===
  # d(phi)/dt
  [phi_time]
    type = TimeDerivative
    variable = phi
  []
  # -D * laplacian(phi)
  [phi_diff]
    type = MatDiffusion
    variable = phi
    diffusivity = D_neutron
  []
  # +Sigma_a * phi (absorption, loss term)
  [phi_abs]
    type = CoefReaction
    variable = phi
    coefficient = 0.08   # Sigma_a
  []
  # -nu_Sf * phi (fission source, gain term)
  [phi_fission]
    type = CoefReaction
    variable = phi
    coefficient = -0.10   # -nu_Sigma_f (negative = source)
  []
  # Xe poisoning: +sigma_Xe_eff * Xe * phi (approx as linear: Xe acts on phi)
  # Model as: flux reduction proportional to Xe concentration
  # Use CoupledForce: residual = -coef * v * test => coef = -sigma_Xe_eff
  [phi_xe_poison]
    type = CoupledForce
    variable = phi
    v = Xe
    coef = -0.05   # -sigma_Xe_eff (negative of CoupledForce convention = +loss)
  []

  # === Iodine-135 equation ===
  # d(I)/dt
  [I_time]
    type = TimeDerivative
    variable = I135
  []
  # -lambda_I * I (decay, loss)
  [I_decay]
    type = CoefReaction
    variable = I135
    coefficient = 0.1033    # lambda_I in /hr
  []
  # +gamma_I * Sigma_f * phi (production from fission)
  [I_production]
    type = CoupledForce
    variable = I135
    v = phi
    coef = 0.00305    # gamma_I * Sigma_f = 0.061 * 0.05
  []

  # === Xenon-135 equation ===
  # d(Xe)/dt
  [Xe_time]
    type = TimeDerivative
    variable = Xe
  []
  # -lambda_Xe * Xe (radioactive decay, loss)
  [Xe_decay]
    type = CoefReaction
    variable = Xe
    coefficient = 0.0752   # lambda_Xe in /hr
  []
  # +lambda_I * I (production from I-135 decay)
  [Xe_from_I]
    type = CoupledForce
    variable = Xe
    v = I135
    coef = 0.1033    # lambda_I (source for Xe)
  []
  # +gamma_Xe * Sigma_f * phi (direct fission yield to Xe)
  [Xe_direct]
    type = CoupledForce
    variable = Xe
    v = phi
    coef = 0.00015   # gamma_Xe * Sigma_f = 0.003 * 0.05
  []
  # Xe burnup by neutron absorption: -sigma_Xe_eff * Xe * phi
  # Approximate as linear loss: -sigma_eff * phi_avg * Xe
  # This is handled by the phi equation coupling above
  # For Xe equation, use a small linear loss representing burnup at avg flux
  [Xe_burnup]
    type = CoefReaction
    variable = Xe
    coefficient = 0.04   # approximate burnup rate at nominal flux
  []
[]

[Materials]
  [neutron_diff]
    type = GenericConstantMaterial
    prop_names = 'D_neutron'
    prop_values = '1.0'
  []
[]

[BCs]
  # Vacuum BCs for neutron flux
  [left_phi]
    type = DirichletBC
    variable = phi
    boundary = left
    value = 0
  []
  [right_phi]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0
  []
  # Xe and I have no spatial BCs (volumetric only, no diffusion)
  # But we need dummy BCs or natural (zero flux) BCs — natural BC is fine
[]

[Postprocessors]
  # Peak neutron flux
  [phi_max]
    type = ElementExtremeValue
    variable = phi
    value_type = max
  []
  # Average Xe concentration
  [Xe_avg]
    type = ElementAverageValue
    variable = Xe
  []
  # Average I-135 concentration
  [I_avg]
    type = ElementAverageValue
    variable = I135
  []
  # Flux at center
  [phi_center]
    type = PointValue
    variable = phi
    point = '10 0.25 0'
  []
  # Xe at center
  [Xe_center]
    type = PointValue
    variable = Xe
    point = '10 0.25 0'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15
  l_max_its = 50

  # 24 hours of reactor operation
  dt = 0.5        # 0.5 hour time steps
  end_time = 24.0  # 24 hours total
[]

[Outputs]
  exodus = true
  csv    = true
[]
