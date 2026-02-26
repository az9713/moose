# ============================================================
# Case 19: Darcy Flow with Heat Transport in Porous Media
#
# PorousFlowBasicTHM action sets up a fully-coupled
# thermo-hydro (TH) system in a single-phase saturated
# 2-D porous medium.
#
# Variables:
#   porepressure  [Pa]  -- fluid pressure
#   temperature   [K]   -- mixture (fluid + solid) temperature
#
# Domain: 3 m wide x 2 m tall rectangular slab
# Flow:   pressure gradient left-to-right (1.1 MPa -> 1.0 MPa)
# Heat:   hot fluid injected from left (350 K), domain initially 300 K
# Goal:   watch the thermal plume advect rightward over 5000 s
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 20
  xmin = 0
  xmax = 3      # 3 m wide
  ymin = 0
  ymax = 2      # 2 m tall
[]

# PorousFlowDictator is a UserObject that keeps track of
# PorousFlow variable numbering. All PorousFlow objects
# must reference it through this GlobalParam.
[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 1e6   # 1 MPa initial pressure everywhere
  []
  [temperature]
    initial_condition = 300   # 300 K initial temperature everywhere
    # Scaling improves conditioning when pressure and temperature
    # appear in the same residual vector (different orders of magnitude).
    scaling = 1e-6
  []
[]

# PorousFlowBasicTHM is a MOOSE Action that automatically adds:
#   - PorousFlowDictator UserObject
#   - Mass-balance kernel for fluid (Darcy flow)
#   - Energy-balance kernel for heat (advection + conduction)
#   - Time-derivative kernels for both variables
#   - The PorousFlow material hierarchy linking fluid properties
# This dramatically reduces the number of objects the user must
# specify compared to building the system kernel-by-kernel.
[PorousFlowBasicTHM]
  porepressure    = porepressure
  temperature     = temperature
  coupling_type   = ThermoHydro   # coupled pressure + temperature
  gravity         = '0 0 0'       # no gravity (horizontal flow)
  fp              = simple_fluid  # reference to FluidProperties block
  multiply_by_density = true      # conserve fluid mass, not fluid volume
[]

# SimpleFluidProperties provides an analytical equation of state
# for a compressible liquid.  For water-like properties:
#   density0 = 1000 kg/m^3 at reference conditions
#   viscosity = 0.001 Pa.s (1 centipoise, liquid water at ~20 C)
#   cp = cv = 4186 J/(kg K)  (specific heat of water)
#   thermal_conductivity = 0.6 W/(m K)
[FluidProperties]
  [simple_fluid]
    type                = SimpleFluidProperties
    density0            = 1000      # kg/m^3
    viscosity           = 0.001     # Pa.s
    thermal_expansion   = 0         # no thermal expansion for simplicity
    cp                  = 4186      # J/(kg K)
    cv                  = 4186      # J/(kg K)
    thermal_conductivity = 0.6      # W/(m K)
  []
[]

[Materials]
  # Porosity: 30% void space filled with fluid.
  # PorousFlowPorosity allows coupling to pressure/temperature changes
  # (here mechanical = false, thermal = false → effectively constant).
  [porosity]
    type          = PorousFlowPorosity
    porosity_zero = 0.3
    mechanical    = false
    thermal       = false
    fluid         = false
  []

  # Biot modulus for the mass-balance time derivative.
  # With biot_coefficient = 1 and very small solid_bulk_compliance,
  # the storage is dominated by fluid compressibility.
  [biot_modulus]
    type                 = PorousFlowConstantBiotModulus
    biot_coefficient     = 1.0
    solid_bulk_compliance = 1e-10   # nearly rigid solid skeleton
    fluid_bulk_modulus   = 2e9      # water bulk modulus
  []

  # Thermal expansion coefficient for the mass-balance equation.
  # Set both coefficients to zero → no pore-volume thermal expansion,
  # keeping the problem purely hydro-thermal without poromechanics.
  [thermal_expansion]
    type                = PorousFlowConstantThermalExpansionCoefficient
    biot_coefficient    = 1.0
    drained_coefficient = 0.0     # no skeleton thermal expansion
    fluid_coefficient   = 0.0     # no fluid thermal expansion
  []

  # Isotropic permeability: k = 1e-12 m^2 (~1 millidarcy)
  # This is typical for a sandstone reservoir.
  # The 3x3 tensor is given row-by-row as 9 values.
  [permeability]
    type        = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0  0 1e-12 0  0 0 1e-12'
  []

  # Solid (rock) thermal energy storage.
  # PorousFlowMatrixInternalEnergy contributes rho_s * c_s * T to the
  # energy balance, representing heat stored in the rock matrix.
  [rock_heat]
    type               = PorousFlowMatrixInternalEnergy
    density            = 2600    # kg/m^3 (typical sandstone)
    specific_heat_capacity = 800 # J/(kg K)
  []

  # Effective (bulk) thermal conductivity of the saturated medium,
  # combining fluid and solid contributions via a mixing rule.
  # Given as a 3x3 tensor; here isotropic: lambda = 2.0 W/(m K).
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.0 0 0  0 2.0 0  0 0 2.0'
  []
[]

[BCs]
  # Pressure boundary conditions drive the flow.
  # The pressure difference dP = 0.1 MPa over 3 m gives a Darcy
  # velocity: q = -k/mu * dP/dx = 1e-12/0.001 * 0.1e6/3 ~ 3.3e-8 m/s
  # At this speed, a fluid parcel takes ~3/(3.3e-8) ~ 91000 s to
  # cross the domain; the thermal retardation (heat stored in rock)
  # makes the thermal front even slower.
  [p_inlet]
    type     = DirichletBC
    variable = porepressure
    boundary = left
    value    = 1.1e6   # 1.1 MPa -- high-pressure inlet
  []
  [p_outlet]
    type     = DirichletBC
    variable = porepressure
    boundary = right
    value    = 1.0e6   # 1.0 MPa -- low-pressure outlet
  []

  # Temperature boundary conditions model hot-fluid injection.
  # Left (inlet): hot fluid at 350 K continuously injected.
  # Right (outlet): ambient temperature maintained at 300 K.
  # Top and bottom are adiabatic (no-flux Neumann, automatic).
  [T_inlet]
    type     = DirichletBC
    variable = temperature
    boundary = left
    value    = 350   # K  -- hot injection temperature
  []
  [T_outlet]
    type     = DirichletBC
    variable = temperature
    boundary = right
    value    = 300   # K  -- ambient (cool) outlet temperature
  []
[]

[Postprocessors]
  # Domain-averaged temperature: rises as the thermal plume fills the domain
  [avg_T]
    type     = ElementAverageValue
    variable = temperature
  []

  # Peak temperature in the domain (tracks where the front is hottest)
  [max_T]
    type = ElementExtremeValue
    variable = temperature
  []

  # Average inlet pressure (should stay near 1.1 MPa)
  [p_inlet_avg]
    type     = SideAverageValue
    variable = porepressure
    boundary = left
  []

  # Average outlet pressure (should stay near 1.0 MPa)
  [p_outlet_avg]
    type     = SideAverageValue
    variable = porepressure
    boundary = right
  []
[]

[Executioner]
  type       = Transient
  solve_type = NEWTON   # full Newton with exact Jacobian from PorousFlow

  # LU factorization (via MUMPS if available, otherwise PETSc's built-in LU).
  # For small 2-D problems, direct solvers are robust and fast.
  # If MUMPS is unavailable, replace with:
  #   petsc_options_iname = '-pc_type -pc_hypre_type'
  #   petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  dt       = 100       # 100 s per step
  end_time = 5000      # 5000 s total (~1.4 hours)

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true   # spatial fields (pressure, temperature) for ParaView
  csv    = true   # postprocessor time history
[]
