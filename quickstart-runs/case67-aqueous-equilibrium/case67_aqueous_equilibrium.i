# ============================================================
# Case 67: Aqueous Equilibrium — CO2-H2O System with pH
# Models the batch equilibration of dissolved CO2 in water.
# Three instantaneous aqueous equilibrium reactions determine
# the speciation of the carbonate system:
#
#   H+ + HCO3-  ⇌  CO2(aq)    log10(Keq) = 6.3447
#   HCO3-       ⇌  H+ + CO3²⁻  log10(Keq) = -10.3288
#   H2O         ⇌  H+ + OH⁻    log10(Keq) = -13.9951
#
# Primary species: H+ and HCO3-
# Secondary species: CO2(aq), CO3²⁻, OH⁻
#
# The ReactionNetwork/AqueousEquilibriumReactions action auto-
# creates AqueousEquilibriumRxnAux for each secondary species
# and coupling kernels for primary species equations.
#
# Initial conditions: [H+] = 1e-5 M, [HCO3-] = 1e-5 M
# Expected equilibrium: pH ≈ 6.6 (slightly acidic from CO2)
#
# Tracked quantities:
#   - pH via PHAux = -log10([H+])
#   - Total concentrations via TotalConcentrationAux
#   - All species concentrations
#
# Domain: 1x1 unit mesh (batch — no spatial variation)
# This is a zero-dimensional chemistry problem solved as a
# relaxation to equilibrium.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 1
  ny   = 1
[]

[Variables]
  [h+]
    initial_condition = 1e-5
  []
  [hco3-]
    initial_condition = 1e-5
  []
[]

[AuxVariables]
  [ph]
  []
  [total_h+]
  []
  [total_hco3-]
  []
[]

[AuxKernels]
  [ph_aux]
    type = PHAux
    variable = ph
    h_conc = h+
  []
  [total_h_aux]
    type = TotalConcentrationAux
    variable = total_h+
    primary_species = h+
    v = 'oh- co3-- co2_aq'
    sto_v = '-1 1 1'
  []
  [total_hco3_aux]
    type = TotalConcentrationAux
    variable = total_hco3-
    primary_species = hco3-
    v = 'co2_aq co3--'
    sto_v = '1 1'
  []
[]

[ReactionNetwork]
  [AqueousEquilibriumReactions]
    primary_species = 'hco3- h+'
    secondary_species = 'co2_aq co3-- oh-'
    reactions = 'hco3- + h+ = co2_aq 6.3447,
                 hco3- - h+ = co3-- -10.3288,
                 - h+ = oh- -13.9951'
  []
[]

[Kernels]
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
    prop_names  = 'diffusivity porosity conductivity'
    prop_values = '1e-7        0.25     1.0'
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
  [co3_val]
    type = ElementAverageValue
    variable = co3--
    execute_on = 'initial timestep_end'
  []
  [oh_val]
    type = ElementAverageValue
    variable = oh-
    execute_on = 'initial timestep_end'
  []
  [ph_val]
    type = ElementAverageValue
    variable = ph
    execute_on = 'initial timestep_end'
  []
  [total_h_val]
    type = ElementAverageValue
    variable = total_h+
    execute_on = 'initial timestep_end'
  []
  [total_hco3_val]
    type = ElementAverageValue
    variable = total_hco3-
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8
  end_time = 1
  dt = 0.1
[]

[Outputs]
  csv = true
[]
