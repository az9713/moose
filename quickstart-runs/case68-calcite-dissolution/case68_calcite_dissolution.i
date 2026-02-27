# ============================================================
# Case 68: Calcite Dissolution — Reactive Geochemistry
# Models the batch dissolution of calcite (CaCO3) in water.
# Calcium and bicarbonate are released into solution as the
# mineral dissolves. The pH rises as the system equilibrates.
#
# Reaction network:
#   Aqueous equilibrium (instantaneous):
#     H+ + HCO3- ⇌ CO2(aq)          log10(Keq) = 6.3447
#     HCO3- ⇌ H+ + CO3²⁻            log10(Keq) = -10.3288
#     Ca²+ + HCO3- ⇌ H+ + CaCO3(aq)  log10(Keq) = -7.0017
#     Ca²+ + HCO3- ⇌ CaHCO3+         log10(Keq) = -1.0467
#     Ca²+ ⇌ H+ + CaOH+              log10(Keq) = -12.85
#     H2O ⇌ H+ + OH⁻                 log10(Keq) = -13.9951
#
#   Kinetic (mineral dissolution):
#     CaCO3(s) → Ca²+ + HCO3- - H+
#     Keq = 10^1.8487, k = 6.46e-7 mol/m²/s
#     A = 0.1 m²/L (reactive surface area)
#
# Primary species: Ca²+, H+, HCO3-
# Secondary (equilibrium): CO2(aq), CO3²⁻, CaCO3(aq),
#                           CaHCO3+, CaOH+, OH⁻
# Kinetic mineral: CaCO3(s)
#
# Initial conditions:
#   Ca²+ = 1e-5 M, H+ = 1e-6 M (pH~6), HCO3- = 1e-5 M
#   CaCO3(s) = 0.05 mol/L (excess solid present)
#
# Expected: Ca²+ increases, pH rises as CaCO3 dissolves
#
# Domain: 1x1 mesh (batch reactor — no spatial variation)
# Run for 100 s to observe dissolution kinetics
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 1
  ny   = 1
[]

[Variables]
  [ca++]
    initial_condition = 1.0e-5
  []
  [h+]
    initial_condition = 1.0e-6
  []
  [hco3-]
    initial_condition = 1.0e-5
  []
[]

[AuxVariables]
  [caco3_s]
    initial_condition = 0.05    # excess solid
  []
  [ph]
  []
[]

[AuxKernels]
  [ph_aux]
    type = PHAux
    h_conc = h+
    variable = ph
  []
[]

[ReactionNetwork]
  [AqueousEquilibriumReactions]
    primary_species = 'ca++ hco3- h+'
    secondary_species = 'co2_aq co3-- caco3_aq cahco3+ caoh+ oh-'
    reactions = 'h+ + hco3- = co2_aq 6.3447,
                 hco3- - h+ = co3-- -10.3288,
                 ca++ + hco3- - h+ = caco3_aq -7.0017,
                 ca++ + hco3- = cahco3+ -1.0467,
                 ca++ - h+ = caoh+ -12.85,
                 - h+ = oh- -13.9951'
  []
  [SolidKineticReactions]
    primary_species = 'ca++ hco3- h+'
    kin_reactions = 'ca++ + hco3- - h+ = caco3_s'
    secondary_species = caco3_s
    log10_keq = 1.8487
    reference_temperature = 298.15
    system_temperature = 298.15
    gas_constant = 8.314
    specific_reactive_surface_area = 0.1
    kinetic_rate_constant = 6.456542e-7
    activation_energy = 1.5e4
  []
[]

[Kernels]
  [ca_time]
    type = PrimaryTimeDerivative
    variable = ca++
  []
  [h_time]
    type = PrimaryTimeDerivative
    variable = h+
  []
  [hco3_time]
    type = PrimaryTimeDerivative
    variable = hco3-
  []
[]

[Materials]
  [porous]
    type = GenericConstantMaterial
    prop_names  = 'porosity diffusivity conductivity'
    prop_values = '0.25     1e-9        1.0'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [h_plus]
    type = ElementAverageValue
    variable = h+
    execute_on = 'initial timestep_end'
  []
  [ca_plus]
    type = ElementAverageValue
    variable = ca++
    execute_on = 'initial timestep_end'
  []
  [hco3_minus]
    type = ElementAverageValue
    variable = hco3-
    execute_on = 'initial timestep_end'
  []
  [co2_aq_val]
    type = ElementAverageValue
    variable = co2_aq
    execute_on = 'initial timestep_end'
  []
  [oh_val]
    type = ElementAverageValue
    variable = oh-
    execute_on = 'initial timestep_end'
  []
  [co3_val]
    type = ElementAverageValue
    variable = co3--
    execute_on = 'initial timestep_end'
  []
  [caco3_s_val]
    type = ElementAverageValue
    variable = caco3_s
    execute_on = 'initial timestep_end'
  []
  [ph_val]
    type = ElementAverageValue
    variable = ph
    execute_on = 'initial timestep_end'
  []
  [calcite_vf]
    type = TotalMineralVolumeFraction
    variable = caco3_s
    molar_volume = 36.934e-6
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  end_time = 100
  dt = 10
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8
[]

[Outputs]
  csv = true
[]
