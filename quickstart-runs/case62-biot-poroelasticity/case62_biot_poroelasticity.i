# ============================================================
# Case 62: Biot Poroelasticity — Coupled Consolidation Column
# A saturated porous column is loaded on top. The load is
# initially carried entirely by the pore fluid. As the fluid
# drains through the top surface, the pore pressure dissipates
# and the effective stress in the solid skeleton increases.
#
# This is the classic Terzaghi consolidation problem solved
# with full Biot poroelasticity (coupled displacement + pressure).
#
# Governing equations:
#   Equilibrium: div(sigma) = 0  (sigma = sigma' - alpha*p*I)
#   Flow: S*dp/dt + alpha*d(eps_v)/dt + div(q) = 0
#   where q = -k/mu * grad(p)
#
# Parameters (simplified for educational demo):
#   E = 100 MPa, nu = 0.25, alpha_Biot = 0.6
#   k = 1e-14 m^2, mu = 0.001 Pa.s
#   Applied load = 10 kPa on top
#   Column: 1 m x 10 m
#
# Expected: initial pore pressure = alpha*sigma_applied,
#   decays exponentially with consolidation coefficient
#   c_v = k/(mu*S_confined)
#
# Domain: [0, 1] x [0, 10] m, 3x30 elements
# ============================================================

[GlobalParams]
  PorousFlowDictator = dictator
  displacements = 'disp_x disp_y'
  biot_coefficient = 0.6
[]

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 3
  ny   = 30
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 10
[]

[Variables]
  [porepressure]
    initial_condition = 6000    # alpha * sigma_applied = 0.6 * 10000
  []
  [disp_x]
  []
  [disp_y]
  []
[]

[Kernels]
  # Solid mechanics — stress divergence
  [stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
  []
  # Pore pressure coupling in equilibrium equation
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
  []
  # Fluid mass balance
  [mass_time]
    type = PorousFlowMassTimeDerivative
    variable = porepressure
    fluid_component = 0
  []
  [mass_flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    fluid_component = 0
    gravity = '0 0 0'
  []
  # Volumetric strain coupling in mass balance
  [vol_strain]
    type = PorousFlowMassVolumetricExpansion
    variable = porepressure
    fluid_component = 0
  []
[]

[AuxVariables]
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_yy_aux]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[FluidProperties]
  [water]
    type = SimpleFluidProperties
    density0            = 1000
    viscosity           = 0.001
    thermal_expansion   = 0
    cp                  = 4186
    cv                  = 4186
    bulk_modulus        = 2e9
  []
[]

[Materials]
  # Elastic properties
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e8      # 100 MPa
    poissons_ratio = 0.25
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []

  # PorousFlow materials
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [fluid_props]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.3
    mechanical    = true
    thermal       = false
    fluid         = true
    biot_coefficient = 0.6
    solid_bulk = 1.667e8   # K_drained = E/(3*(1-2*nu)) = 1e8/(3*0.5)
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-14 0 0  0 1e-14 0  0 0 1e-14'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [eff_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [near_zero_thermal_exp]
    type = PorousFlowConstantThermalExpansionCoefficient
    biot_coefficient = 0.6
    drained_coefficient = 0
    fluid_coefficient = 0
  []
  [biot_mod]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.6
    solid_bulk_compliance = 6e-9   # 1/K_s
    fluid_bulk_modulus = 2e9
  []
[]

[BCs]
  # Bottom: fixed
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  # Sides: roller
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_right_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  # Top: applied load via Neumann BC (10 kPa = 10000 Pa downward)
  [top_load]
    type = NeumannBC
    variable = disp_y
    boundary = top
    value = -1e4     # 10 kPa compressive (negative = downward)
  []
  # Top: drained (pore pressure = 0)
  [top_drained]
    type = DirichletBC
    variable = porepressure
    boundary = top
    value = 0
  []
[]

[Postprocessors]
  # Pore pressure at base
  [p_base]
    type = PointValue
    variable = porepressure
    point = '0.5 0 0'
  []
  # Pore pressure at mid-height
  [p_mid]
    type = PointValue
    variable = porepressure
    point = '0.5 5 0'
  []
  # Vertical displacement at top (settlement)
  [disp_y_top]
    type = PointValue
    variable = disp_y
    point = '0.5 10 0'
  []
  # Average pore pressure
  [p_avg]
    type = ElementAverageValue
    variable = porepressure
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  dt = 5e5         # 500000 s time steps (slow consolidation with low k)
  end_time = 2.5e7 # 25 million seconds (~290 days)

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
[]

[Outputs]
  exodus = true
  csv    = true
[]
