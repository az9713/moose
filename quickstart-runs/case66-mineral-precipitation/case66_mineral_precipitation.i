# ============================================================
# Case 66: Mineral Precipitation — Two-Species Kinetic Reaction
# Two dissolved species A and B diffuse from opposite ends of
# a porous medium. Where they meet, they react kinetically to
# form a solid mineral precipitate: A + B → mineral(s).
#
# This uses the chemical_reactions module's SolidKineticReactions
# action, which auto-creates:
#   - AuxVariable for mineral concentration
#   - KineticDisPreConcAux (updates mineral concentration)
#   - CoupledBEKinetic kernels (coupling to primary species)
#
# The kinetic rate follows an Arrhenius law:
#   rate = k * exp(-Ea/(R*T)) * A_s * (1 - Q/Keq)
# where Q is the ion activity product and Keq is the
# equilibrium constant.
#
# Parameters:
#   D = 5e-4 m²/s, porosity = 0.4
#   log10(Keq) = -6, k = 1e-8 mol/m²/s
#   Ea = 15 kJ/mol, T = 298.15 K
#   A left = 0.01 mol/L, B right = 0.01 mol/L
#
# Domain: [0, 1] x [0, 0.2], 40x4 elements
# Run for 50 s to see mineral band form where A meets B
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 40
  ny   = 4
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.2
[]

[Variables]
  [a]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
  [b]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[ReactionNetwork]
  [SolidKineticReactions]
    primary_species = 'a b'
    secondary_species = mineral
    kin_reactions = 'a + b = mineral'
    log10_keq = '-6'
    specific_reactive_surface_area = '1.0'
    kinetic_rate_constant = '1.0e-8'
    activation_energy = '1.5e4'
    gas_constant = 8.314
    reference_temperature = '298.15'
    system_temperature = '298.15'
  []
[]

[Kernels]
  # Species A: time derivative + diffusion
  [a_time]
    type = PrimaryTimeDerivative
    variable = a
  []
  [a_diff]
    type = PrimaryDiffusion
    variable = a
  []
  # Species B: time derivative + diffusion
  [b_time]
    type = PrimaryTimeDerivative
    variable = b
  []
  [b_diff]
    type = PrimaryDiffusion
    variable = b
  []
[]

[BCs]
  # A from left, zero on right
  [a_left]
    type = DirichletBC
    variable = a
    preset = false
    boundary = left
    value = 1.0e-2
  []
  [a_right]
    type = DirichletBC
    variable = a
    preset = false
    boundary = right
    value = 0
  []
  # B from right, zero on left
  [b_left]
    type = DirichletBC
    variable = b
    preset = false
    boundary = left
    value = 0
  []
  [b_right]
    type = DirichletBC
    variable = b
    preset = false
    boundary = right
    value = 1.0e-2
  []
[]

[Materials]
  [porous]
    type = GenericConstantMaterial
    prop_names  = 'diffusivity conductivity porosity'
    prop_values = '5e-4        4e-3         0.4'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  # Species A at midpoint
  [a_mid]
    type = PointValue
    variable = a
    point = '0.5 0.1 0'
  []
  # Species B at midpoint
  [b_mid]
    type = PointValue
    variable = b
    point = '0.5 0.1 0'
  []
  # Mineral concentration at midpoint (where reaction occurs)
  [mineral_mid]
    type = PointValue
    variable = mineral
    point = '0.5 0.1 0'
  []
  # Average mineral
  [mineral_avg]
    type = ElementAverageValue
    variable = mineral
  []
  # Average species concentrations
  [a_avg]
    type = ElementAverageValue
    variable = a
  []
  [b_avg]
    type = ElementAverageValue
    variable = b
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  end_time = 50
  dt       = 5

  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv    = true
[]
