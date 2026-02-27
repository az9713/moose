# ============================================================
# Case 59: Terzaghi 1D Consolidation
# Classic soil mechanics problem: a saturated soil column
# with an applied surface load. The excess pore pressure
# dissipates over time as fluid drains from the top.
#
# Governing equation (1D, vertical):
#   d(p)/dt = c_v * d^2(p)/dz^2
# where c_v = k/(mu * S) is the consolidation coefficient,
#   k = permeability, mu = fluid viscosity, S = storage coeff.
#
# Analytical solution (Terzaghi, 1925):
#   p(z,t) = sum_{n=0}^{inf} (4*p0/(2n+1)/pi) *
#            sin((2n+1)*pi*z/(2*H)) * exp(-(2n+1)^2*pi^2*c_v*t/(4*H^2))
#
# Parameters:
#   H = 10 m column height
#   k = 1e-11 m^2 (low-permeability clay)
#   mu = 0.001 Pa.s (water viscosity)
#   Initial excess pore pressure p0 = 100 kPa (from surface load)
#   Fluid bulk modulus = 2e9 Pa
#   Porosity = 0.4
#   Storage S ~ phi/K_f + (1-phi)/K_s ≈ 2e-10 /Pa
#   c_v = k/(mu*S) ≈ 1e-11/(0.001*2e-10) = 0.05 m^2/s
#   Characteristic time: t_c = H^2/c_v = 100/0.05 = 2000 s
#
# BCs: top = drained (p=0), bottom = undrained (zero flux)
# Domain: [0, 10] m, 50 elements, quasi-1D
# Run for 5000 s to see near-complete dissipation
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 2
  ny   = 50
  xmin = 0
  xmax = 0.5
  ymin = 0
  ymax = 10      # 10 m column height
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 1e5   # 100 kPa initial excess pore pressure
  []
[]

[PorousFlowBasicTHM]
  porepressure    = porepressure
  coupling_type   = Hydro       # pressure only (no temperature)
  gravity         = '0 0 0'     # no gravity (pure pressure diffusion)
  fp              = simple_fluid
  multiply_by_density = true
[]

[FluidProperties]
  [simple_fluid]
    type                 = SimpleFluidProperties
    density0             = 1000      # kg/m^3
    viscosity            = 0.001     # Pa.s
    thermal_expansion    = 0
    cp                   = 4186
    cv                   = 4186
    thermal_conductivity = 0.6
  []
[]

[Materials]
  [porosity]
    type          = PorousFlowPorosity
    porosity_zero = 0.4
    mechanical    = false
    thermal       = false
    fluid         = false
  []
  [biot_modulus]
    type                  = PorousFlowConstantBiotModulus
    biot_coefficient      = 1.0
    solid_bulk_compliance = 1e-10
    fluid_bulk_modulus    = 2e9
  []
  [thermal_expansion]
    type                 = PorousFlowConstantThermalExpansionCoefficient
    biot_coefficient     = 1.0
    drained_coefficient  = 0.0
    fluid_coefficient    = 0.0
  []
  [permeability]
    type         = PorousFlowPermeabilityConst
    permeability = '1e-15 0 0  0 1e-15 0  0 0 1e-15'
  []
  # Note: PorousFlowMatrixInternalEnergy and ThermalConductivityIdeal
  # are NOT needed for Hydro-only coupling (they require temperature)
[]

[BCs]
  # Top boundary: drained (excess pore pressure = 0)
  [top_drained]
    type     = DirichletBC
    variable = porepressure
    boundary = top
    value    = 0      # excess pore pressure dissipates to zero at top
  []
  # Bottom: undrained (natural zero-flux BC, no block needed)
  # Left/right: symmetry (natural zero-flux BC)
[]

[Postprocessors]
  # Pore pressure at base (should start at 100 kPa, decay to 0)
  [p_base]
    type     = PointValue
    variable = porepressure
    point    = '0.25 0 0'
  []
  # Pore pressure at mid-height
  [p_mid]
    type     = PointValue
    variable = porepressure
    point    = '0.25 5 0'
  []
  # Average pore pressure (tracks overall consolidation)
  [p_avg]
    type     = ElementAverageValue
    variable = porepressure
  []
  # Degree of consolidation: U = 1 - p_avg/p0
  # (computed in post-processing from p_avg/1e5)
[]

[Executioner]
  type       = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  dt       = 100       # 100 s time steps
  end_time = 5000      # 5000 s total (2.5 * t_c)

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
