# MOOSE Modules Reference

**Framework**: MOOSE (Multiphysics Object-Oriented Simulation Environment)
**Language**: C++17, finite-element / finite-volume
**Source**: `modules/` in the MOOSE repository
**Last updated**: 2026-02-24

---

## Table of Contents

1. [How Modules Work](#1-how-modules-work)
2. [Module Dependency Graph](#2-module-dependency-graph)
3. [Which Module Do I Need? — Decision Guide](#3-which-module-do-i-need--decision-guide)
4. [Module Reference](#4-module-reference)
   - [chemical_reactions](#41-chemical_reactions)
   - [combined](#42-combined)
   - [contact](#43-contact)
   - [electromagnetics](#44-electromagnetics)
   - [external_petsc_solver](#45-external_petsc_solver)
   - [fluid_properties](#46-fluid_properties)
   - [fsi](#47-fsi-fluid-structure-interaction)
   - [functional_expansion_tools](#48-functional_expansion_tools)
   - [geochemistry](#49-geochemistry)
   - [heat_transfer](#410-heat_transfer)
   - [level_set](#411-level_set)
   - [misc](#412-misc)
   - [navier_stokes](#413-navier_stokes)
   - [optimization](#414-optimization)
   - [peridynamics](#415-peridynamics)
   - [phase_field](#416-phase_field)
   - [porous_flow](#417-porous_flow)
   - [ray_tracing](#418-ray_tracing)
   - [rdg](#419-rdg-reconstructed-discontinuous-galerkin)
   - [reactor](#420-reactor)
   - [richards](#421-richards)
   - [scalar_transport](#422-scalar_transport)
   - [solid_mechanics](#423-solid_mechanics-formerly-tensor_mechanics)
   - [solid_properties](#424-solid_properties)
   - [stochastic_tools](#425-stochastic_tools)
   - [subchannel](#426-subchannel)
   - [tensor_mechanics](#427-tensor_mechanics-alias)
   - [thermal_hydraulics](#428-thermal_hydraulics)
   - [xfem](#429-xfem)
5. [Cross-Module Coupling Patterns](#5-cross-module-coupling-patterns)

---

## 1. How Modules Work

Each MOOSE module is a self-contained C++ library that registers its own kernels, materials, boundary conditions, user objects, and other MOOSE objects with the framework. Modules are enabled at compile time by setting make variables, for example:

```makefile
HEAT_TRANSFER := yes
SOLID_MECHANICS := yes
```

Dependency resolution is automatic: the file `modules/modules.mk` computes the transitive closure of all required modules. For example, enabling `POROUS_FLOW` automatically enables `SOLID_MECHANICS`, `FLUID_PROPERTIES`, and `CHEMICAL_REACTIONS`.

Each module has a Makefile suffix (used in the object suffix chain) listed in `modules.mk`. The suffixes are:

| Module | Suffix |
|---|---|
| chemical_reactions | cr |
| contact | con |
| electromagnetics | em |
| external_petsc_solver | eps |
| fluid_properties | fp |
| fsi | fsi |
| functional_expansion_tools | fet |
| geochemistry | gc |
| heat_transfer | ht |
| level_set | ls |
| misc | misc |
| navier_stokes | ns |
| optimization | opt |
| peridynamics | pd |
| phase_field | pf |
| porous_flow | pflow |
| ray_tracing | ray |
| rdg | rdg |
| reactor | rct |
| richards | rich |
| scalar_transport | st |
| solid_mechanics | sm |
| solid_properties | sp |
| stochastic_tools | st |
| subchannel | sc |
| thermal_hydraulics | th |
| xfem | xfem |

---

## 2. Module Dependency Graph

The arrows point from a module to the modules it requires. Modules with no incoming arrows are leaf dependencies.

```
misc ─────────────────────────────────────────────────────────────────────────┐
  │                                                                            │
  └──► fluid_properties ──────────────────────────────────────────────────────┤
            │                                                                  │
            ├──► heat_transfer ────────────────────────────────────────────────┤
            │         │                                                        │
            │         ├──► solid_properties ───────────────────────────────────┤
            │         │                                                        │
            │         └──► ray_tracing (leaf) ──────────────────────────────── ┤
            │                                                                  │
            └──► rdg (leaf) ─────────────────────────────────────────────────── ┤

solid_mechanics (leaf) ────────────────────────────────────────────────────────┤
     │                                                                         │
     ├──► contact ──────────────────────────────────────────────────────────── ┤
     ├──► phase_field ──────────────────────────────────────────────────────── ┤
     ├──► peridynamics ─────────────────────────────────────────────────────── ┤
     └──► xfem ────────────────────────────────────────────────────────────── ┤

chemical_reactions (leaf) ────────────────────────────────────────────────────┤

electromagnetics (leaf) ──────────────────────────────────────────────────────┤
functional_expansion_tools (leaf) ────────────────────────────────────────────┤
geochemistry (leaf) ──────────────────────────────────────────────────────────┤
level_set (leaf) ─────────────────────────────────────────────────────────────┤
optimization (leaf) ──────────────────────────────────────────────────────────┤
reactor (leaf) ───────────────────────────────────────────────────────────────┤
richards (leaf) ──────────────────────────────────────────────────────────────┤
stochastic_tools (leaf) ──────────────────────────────────────────────────────┤
external_petsc_solver (leaf) ─────────────────────────────────────────────────┘

Composite modules (depend on multiple leaf/mid-tier modules):

navier_stokes ──► fluid_properties + heat_transfer + rdg
     │
     └──► fsi ──► navier_stokes + solid_mechanics

porous_flow ──► solid_mechanics + fluid_properties + chemical_reactions

thermal_hydraulics ──► navier_stokes + fluid_properties + heat_transfer
                       + rdg + ray_tracing + solid_properties + misc

scalar_transport ──► chemical_reactions + navier_stokes + thermal_hydraulics
                     + fluid_properties + heat_transfer + rdg + ray_tracing
                     + solid_properties + misc

subchannel ──► fluid_properties + heat_transfer + reactor

combined ──► ALL modules
```

**Compact notation** (direct declared dependencies from `modules.mk`):

```
fluid_properties      : misc
heat_transfer         : ray_tracing
solid_properties      : heat_transfer
solid_mechanics       : (none)
contact               : solid_mechanics
peridynamics          : solid_mechanics
phase_field           : solid_mechanics
xfem                  : solid_mechanics
navier_stokes         : fluid_properties, rdg, heat_transfer
fsi                   : navier_stokes, solid_mechanics
porous_flow           : solid_mechanics, fluid_properties, chemical_reactions
thermal_hydraulics    : navier_stokes, fluid_properties, heat_transfer,
                        rdg, ray_tracing, solid_properties, misc
scalar_transport      : chemical_reactions, navier_stokes, thermal_hydraulics,
                        fluid_properties, heat_transfer, rdg, ray_tracing,
                        solid_properties, misc
subchannel            : fluid_properties, heat_transfer, reactor
```

---

## 3. Which Module Do I Need? — Decision Guide

| Physics you want to simulate | Primary module(s) | Also enable |
|---|---|---|
| Heat conduction in a solid | heat_transfer | — |
| Convective heat transfer to fluid | heat_transfer | — |
| Radiation heat transfer (view factors) | heat_transfer | ray_tracing |
| Joule heating (resistive) | heat_transfer, electromagnetics | — |
| Material thermal properties (Cp, k) from library | solid_properties | heat_transfer |
| Linear elastic statics or dynamics | solid_mechanics | — |
| Plasticity, creep, damage | solid_mechanics | — |
| Fracture mechanics (J-integral, SIF) | solid_mechanics | — |
| Cohesive zone fracture | solid_mechanics | — |
| Crack propagation with mesh cutting | xfem | solid_mechanics |
| Frictionless or frictional contact | contact | solid_mechanics |
| Fluid-structure interaction | fsi | navier_stokes, solid_mechanics |
| Incompressible Navier-Stokes (FE) | navier_stokes | fluid_properties |
| Incompressible Navier-Stokes (FV) | navier_stokes | fluid_properties |
| Compressible Euler / RANS flow | navier_stokes | fluid_properties, rdg |
| Equation of state for a real fluid | fluid_properties | misc |
| 1D system-level thermal hydraulics | thermal_hydraulics | (auto) |
| Nuclear subchannel thermal hydraulics | subchannel | fluid_properties, heat_transfer |
| Multiphase porous media (THM) | porous_flow | (auto) |
| Unsaturated Richards flow (legacy) | richards | — |
| Allen-Cahn / Cahn-Hilliard phase field | phase_field | solid_mechanics |
| Grain growth, microstructure evolution | phase_field | solid_mechanics |
| Aqueous chemistry with equilibrium | chemical_reactions | — |
| Reactive geochemistry (databases) | geochemistry | — |
| Reactive transport in porous media | porous_flow + geochemistry | (auto) |
| Scalar advection-diffusion-reaction | scalar_transport | (auto) |
| Meshfree (bond-based or NOS) mechanics | peridynamics | solid_mechanics |
| Design optimization, inverse problems | optimization | — |
| Topology optimization (SIMP) | optimization | solid_mechanics, heat_transfer |
| Uncertainty quantification, surrogates | stochastic_tools | — |
| Interface tracking (level set) | level_set | — |
| Electrostatics / magnetostatics | electromagnetics | — |
| Full-wave EM (Maxwell) | electromagnetics | — |
| Nuclear reactor geometry meshing | reactor | — |
| Functional series expansions (coupling) | functional_expansion_tools | — |
| Line integrals along rays | ray_tracing | — |
| DG hyperbolic conservation laws | rdg | — |
| Solid material property library | solid_properties | heat_transfer |
| Coupling to external PETSc solver | external_petsc_solver | — |
| Everything at once (prototyping) | combined | — |

---

## 4. Module Reference

---

### 4.1 chemical_reactions

**Makefile variable**: `CHEMICAL_REACTIONS := yes`
**Suffix**: `cr`
**Direct dependencies**: none

#### Physics

The `chemical_reactions` module handles aqueous geochemical reactions in a
porous medium. It supports:

- **Primary species transport**: Each dissolved primary species `c_i` satisfies
  an advection-diffusion-reaction equation. The governing equation is

  ```
  phi * d(c_i)/dt + phi * div(-D grad c_i + v c_i)
      + sum_k nu_{ik} * R_k(c)  =  0
  ```

  where `phi` is porosity, `D` is diffusivity, `v` is the pore velocity,
  `nu_{ik}` is the stoichiometric coefficient of species `i` in reaction `k`,
  and `R_k` is the reaction rate.

- **Equilibrium secondary species**: Secondary species concentrations are
  algebraically determined from equilibrium constants. For an equilibrium
  reaction `sum_i nu_i * c_i = c_j`, the equilibrium expression is

  ```
  c_j  =  K_eq * prod_i ( c_i^nu_i )
  ```

  where `K_eq` is the equilibrium constant (stored in log10 form in the input
  as `log10(Keq)`).

- **Kinetic mineral dissolution/precipitation** (solid-phase kinetics): The
  kinetic rate follows a transition-state theory form,

  ```
  R_mineral  =  A_s * k_r * ( 1 - Q/K_eq )
  ```

  where `A_s` is reactive surface area, `k_r` is the rate constant, and `Q` is
  the ion activity product.

- **Langmuir sorption**: Sorption of a species to a solid matrix following the
  Langmuir isotherm is handled via `LangmuirMaterial` and
  `DesorptionFromMatrix`/`DesorptionToPorespace` kernels.

- **Thermochimia coupling**: Optional coupling to the Thermochimia
  thermochemical library for more complex reaction networks via
  `ThermochimicaData` user object and `ThermochimicaUtils`.

#### Key Classes

**Actions**
- `AddPrimarySpeciesAction` — registers primary species variables
- `AddSecondarySpeciesAction` — registers auxiliary variables for secondary species
- `AddCoupledEqSpeciesAction` — adds equilibrium reaction kernels from the
  `ReactionNetwork/AqueousEquilibriumReactions` block
- `AddCoupledSolidKinSpeciesAction` — adds kinetic reaction kernels from the
  `ReactionNetwork/SolidKineticReactions` block
- `ChemicalCompositionAction` — high-level action for multi-species composition

**Kernels**
- `PrimaryTimeDerivative` — `phi * dc/dt`
- `PrimaryDiffusion` — `- div(phi * D * grad c)`
- `PrimaryConvection` — `div(phi * v * c)` driven by Darcy pressure
- `CoupledBEEquilibriumSub` — implicit time derivative contribution from
  equilibrium secondary species
- `CoupledDiffusionReactionSub` — diffusive flux contribution of secondary
  species to primary species equation
- `CoupledConvectionReactionSub` — convective flux contribution
- `CoupledBEKinetic` — backward Euler time derivative for kinetic mineral
- `DesorptionFromMatrix` / `DesorptionToPorespace` — Langmuir desorption terms

**AuxKernels**
- `AqueousEquilibriumRxnAux` — computes equilibrium secondary species
  concentration from primary species
- `KineticDisPreConcAux` / `KineticDisPreRateAux` — kinetic mineral
  concentration and rate output
- `PHAux` — computes pH from H+ concentration
- `TotalConcentrationAux` — sum of primary + secondary contributions
- `EquilibriumConstantAux` — evaluates temperature-dependent `log10(Keq)`

**Boundary Conditions**
- `ChemicalOutFlowBC` — Neumann outflow condition for reactive species

**Materials**
- `LangmuirMaterial` / `MollifiedLangmuirMaterial` — Langmuir isotherm
  sorption coefficients

**UserObjects**
- `ThermochimicaData` — interface to Thermochimia library

**Utilities**
- `EquilibriumConstantFit` — polynomial fit for `log10(Keq)(T)`
- `ReactionNetworkUtils` — parsing of reaction string notation

#### Module Dependencies

None. `chemical_reactions` is a leaf module.

#### Minimal Example Input

```
# Single primary species a transported by diffusion and convection
# with equilibrium reaction: 2a = pa2, log10(Keq) = 1

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
[]

[Variables]
  [a]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxVariables]
  [pressure]
    family = LAGRANGE
    order = FIRST
  []
[]

[ReactionNetwork]
  [AqueousEquilibriumReactions]
    primary_species = a
    reactions = '2a = pa2 1'
    secondary_species = pa2
    pressure = pressure
  []
[]

[Kernels]
  [a_dt]
    type = PrimaryTimeDerivative
    variable = a
  []
  [a_diff]
    type = PrimaryDiffusion
    variable = a
  []
  [a_conv]
    type = PrimaryConvection
    variable = a
    p = pressure
  []
[]

[BCs]
  [outflow]
    type = ChemicalOutFlowBC
    variable = a
    boundary = right
  []
[]

[Materials]
  [porous]
    type = GenericConstantMaterial
    prop_names  = 'diffusivity conductivity porosity'
    prop_values = '1e-4        1e-4        0.2'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  end_time = 100
  dt = 10
[]
```

#### When to Use

Use `chemical_reactions` for reactive transport in a porous medium where
species concentrations are field variables (i.e., spatially varying) and you
need the transport operators (advection, diffusion) coupled to the reaction
network. For zero-dimensional geochemical equilibrium problems without spatial
transport, prefer `geochemistry`. For fully coupled THM reactive transport in
geological systems, use `porous_flow` (which depends on `chemical_reactions`).

---

### 4.2 combined

**Makefile variable**: Built automatically when `ALL_MODULES=yes`
**Suffix**: `comb`

The `combined` module is not a physics module — it simply links all other
modules together into a single executable. It exists so that application
developers can build one binary that includes every registered object in the
MOOSE module ecosystem, without writing their own application that explicitly
enables every individual module.

#### When to Use

- Rapid prototyping across multiple physics
- Testing cross-module interactions
- Users who want a single `moose-opt` or `combined-opt` binary

For production applications, enable only the modules you need: smaller
executables link faster, compile faster, and have a smaller memory footprint.
In a normal `modules.mk` build with `ALL_MODULES=yes`, the `combined`
application is automatically assembled at the end.

---

### 4.3 contact

**Makefile variable**: `CONTACT := yes`
**Suffix**: `con`
**Direct dependencies**: `solid_mechanics`

#### Physics

Mechanical contact between deformable bodies. The kinematic constraint is
the non-penetration condition: the gap function `g_n` between a secondary node
and the primary surface must satisfy

```
g_n >= 0
lambda_n >= 0
lambda_n * g_n = 0
```

where `lambda_n` is the normal contact pressure (Lagrange multiplier or penalty
force). For frictional contact, the Coulomb friction law introduces a tangential
constraint:

```
|lambda_t| <= mu * lambda_n
```

Three major formulations are available:

1. **Mortar contact** (recommended): Uses a mortar Lagrange multiplier field
   defined on a lower-dimensional subdomain. Provides accurate gap integration
   and better convergence than node-to-segment methods.
2. **Penalty contact**: Enforces non-penetration via a penalty force
   `f = k_n * max(0, -g_n) * n`. Simple to implement but perturbs the
   exact solution.
3. **Augmented Lagrangian contact**: Iterates between a Lagrange multiplier
   update and a penalty solve. Handled by `AugmentedLagrangianContactProblem`.

Contact also supports **explicit dynamics** via `ExplicitDynamicsContactConstraint`,
used when time integration uses central differences (no global system solve).

#### Key Classes

**Actions**
- `ContactAction` — high-level action that reads the `[Contact]` block and
  creates all required mortar or penalty objects
- `ExplicitDynamicsContactAction` — sets up explicit dynamics contact

**Constraints**
- `ComputeWeightedGapLMMechanicalContact` — computes the weighted gap for
  mortar normal contact
- `ComputeWeightedGapCartesianLMMechanicalContact` — Cartesian variant
- `ComputeFrictionalForceLMMechanicalContact` — mortar frictional contact
  (frictionless Coulomb)
- `ComputeFrictionalForceCartesianLMMechanicalContact` — Cartesian frictional
- `ComputeDynamicWeightedGapLMMechanicalContact` — mortar dynamic contact
- `ComputeDynamicFrictionalForceLMMechanicalContact` — mortar dynamic
  frictional
- `NormalMortarMechanicalContact` — normal mortar constraint
- `TangentialMortarMechanicalContact` — tangential mortar constraint
- `CartesianMortarMechanicalContact` — full Cartesian mortar
- `MechanicalContactConstraint` — legacy node-to-segment penalty constraint
- `RANFSNormalMechanicalContact` — residual-as-nodal-force-scaling normal
  contact
- `MortarGenericTraction` — generic traction via mortar
- `ExplicitDynamicsContactConstraint` — explicit dynamic velocity correction

**UserObjects**
- `WeightedGapUserObject` / `LMWeightedGapUserObject` — integrate gap on
  mortar interface
- `WeightedVelocitiesUserObject` / `LMWeightedVelocitiesUserObject` — integrate
  velocities for frictional slip computation
- `PenaltyWeightedGapUserObject` — penalty-based gap accumulation
- `PenaltyFrictionUserObject` — penalty-based friction
- `BilinearMixedModeCohesiveZoneModel` — bilinear mixed-mode CZM
- `NodalArea` / `NodalDensity` / `NodalWaveSpeed` — utility UOs for explicit dynamics

**AuxKernels**
- `ContactPressureAux` — output contact pressure from LM
- `MortarFrictionalPressureVectorAux` — frictional traction vector
- `MortarFrictionalStateAux` — slip / stick / open status
- `MortarArchardsLawAux` — Archards wear law output

**Problems and Convergence**
- `AugmentedLagrangianContactProblem` — outer AL loop
- `AugmentedLagrangianContactConvergence` — convergence criterion for AL

**Solvers**
- `ContactLineSearchBase` / `PetscContactLineSearch` — contact-aware line search
- `ContactSlipDamper` — damps large slip increments
- `ContactSplit` — field-split preconditioner separating contact DOFs

#### Module Dependencies

`contact` requires `solid_mechanics`. The solid mechanics displacement
variables and stress divergence kernels must already be set up before adding
contact objects.

#### Minimal Example Input

```
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    strain = FINITE
    block = 'top_block bottom_block'
  []
[]

[Contact]
  [mortar]
    primary   = bottom_top
    secondary = top_bottom
    formulation = mortar
    model = coulomb
    friction_coefficient = 0.4
    c_normal     = 1e4
    c_tangential = 1e4
  []
[]

[BCs]
  [fix_bottom]
    type = DirichletBC
    variable = disp_z
    boundary = bottom_bottom
    value = 0
  []
  [push_top]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_top
    function = '-0.25 * sin(2 * pi / 40 * t)'
  []
[]
```

#### When to Use

Use `contact` whenever two solid bodies come into physical contact and the
contact forces must be computed (not just applied as loads). For thermal
contact resistance across a gap use the `GapConductanceConstraint` from
`heat_transfer` instead (or in addition). For fluid pressure acting on a
deformable interface, use `fsi`. Cohesive zone fracture without contact can
also be handled entirely within `solid_mechanics`.

#### Quickstart Example

**Case 70** (`quickstart-runs/case70-mortar-contact`) demonstrates mortar
frictionless contact between a 2D elastic block and a rigid foundation. It
uses `MeshCollectionGenerator` to assemble the two-body mesh, the `[Contact]`
block to activate `ContactAction`, and `ContactPressureAux` to visualize the
Lagrange multiplier pressure field. The peak contact pressure is compared to
the Hertz analytic solution.

---

### 4.4 electromagnetics

**Makefile variable**: `ELECTROMAGNETICS := yes`
**Suffix**: `em`
**Direct dependencies**: none

#### Physics

The `electromagnetics` module solves Maxwell's equations in the frequency and
time domains. It handles:

- **Electrostatics**: Gauss's law `div(epsilon * E) = rho_free`, where `E = -grad phi`,
  solved for the electric potential `phi`.
- **Magnetostatics**: `curl H = J_free` with `B = mu H`, solved with vector
  potentials.
- **Full-wave electromagnetics (time-harmonic)**: The vector Helmholtz equation
  for the electric field `E`:

  ```
  curl curl E - k^2 E = -i omega mu J_source
  ```

  where `k^2 = omega^2 mu epsilon - i omega mu sigma`. In MOOSE notation,
  the weak form uses Nedelec edge elements (H(curl)-conforming) for the
  electric field.

- **Transient wave equation**: Second-order time derivative form,

  ```
  epsilon * d^2 E/dt^2  +  sigma * dE/dt  +  curl( (1/mu) * curl E )  =  J_source
  ```

- **Joule heating coupling**: The Joule heat source `Q_J = sigma |E|^2`
  (or `J . E`) can be passed to a heat conduction solve via
  `EMJouleHeatingSource` / `ADJouleHeatingSource` (also in `heat_transfer`).

#### Key Classes

**Kernels**
- `CurlCurlField` — `curl((1/mu) * curl E)` term (weak form residual)
- `ADMatWaveReaction` — `-k^2 * E` wave reaction term with material coefficient
- `ADConductionCurrent` — `sigma * dE/dt` or `i omega sigma E` (frequency domain)
- `VectorSecondTimeDerivative` — `epsilon * d^2 E/dt^2` for transient wave
- `VectorCurrentSource` — source current term `J_source`
- `EMJouleHeatingSource` — deposits Joule heat `sigma |E|^2` into a temperature
  equation

**Boundary Conditions**
- `EMRobinBC` — first-order absorbing / Robin boundary condition for scalar
  problems
- `VectorEMRobinBC` — first-order absorbing (Mur) BC for vector EM problems
- `VectorTransientAbsorbingBC` — transient perfectly-matched-layer-type BC

**Interface Kernels**
- `ElectrostaticContactCondition` — continuity of normal current density
  across a resistive interface
- `ParallelElectricFieldInterface` — tangential E-field continuity
- `PerpendicularElectricFieldInterface` — normal D-field continuity

**AuxKernels**
- `CurrentDensity` — computes `J = sigma * E` or `J = -sigma * grad phi`
- `PotentialToFieldAux` — converts scalar potential gradient to field vector
- `EMJouleHeatingHeatGeneratedAux` — evaluates the Joule heating density
  `sigma |E|^2` as an aux variable for visualization
- `AzimuthMagneticTimeDerivRZ` — azimuthal component of `dB/dt` for
  axisymmetric induction problems
- `SourceCurrentHeating` — heating rate from external source current

**Materials**
- `WaveEquationCoefficient` — material that packages `mu`, `epsilon`, `sigma`
  and derived wave number `k` for the Helmholtz equation

**Postprocessors**
- `ReflectionCoefficient` — computes the power reflection coefficient from a
  boundary flux

**Utilities**
- `ElectromagneticConstants` — physical constants (epsilon_0, mu_0, c, etc.)
- `ElectromagneticEnums` — enumeration types for field components and modes

#### Module Dependencies

None. `electromagnetics` is a leaf module. For Joule-heating coupled problems,
also enable `heat_transfer`.

#### Minimal Example (Electrostatics + Joule Heating)

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 10
[]

[Variables]
  [phi]   []  # electric potential
  [T]     []  # temperature
[]

[Kernels]
  [electric_div]
    type = ADDiffusion
    variable = phi
  []
  [joule_heat]
    type = ADJouleHeatingSource
    variable = T
    elec = phi
    electrical_conductivity = sigma
  []
  [heat_diff]
    type = ADDiffusion
    variable = T
  []
[]

[BCs]
  [phi_left]
    type = DirichletBC
    variable = phi
    boundary = left
    value = 1.0
  []
  [phi_right]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0.0
  []
[]

[Materials]
  [conductivity]
    type = ADGenericConstantMaterial
    prop_names  = 'electrical_conductivity'
    prop_values = '1.0'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
```

#### When to Use

Use `electromagnetics` when electromagnetic fields are the primary physics
(Maxwell equations, induction, microwave heating, electrostatics). For
purely resistive current flow without wave effects, a simple `Diffusion`
kernel on `phi` plus `JouleHeatingSource` from `heat_transfer` may suffice.
Use `electromagnetics` for the full curl-curl operator, frequency-domain
problems, or when boundary conditions specific to EM (absorbing, Robin) are
needed.

---

### 4.5 external_petsc_solver

**Makefile variable**: `EXTERNAL_PETSC_SOLVER := yes`
**Suffix**: `eps`
**Direct dependencies**: none

#### Physics

A thin wrapper that allows an external PETSc solver (e.g., a finite-difference
code using PETSc DMDA structures) to be embedded inside a MOOSE simulation as
a sub-app or as the primary solver. This is an escape hatch for legacy codes
that already use PETSc directly.

#### Key Classes

- `PETScDMDAMesh` — creates a MOOSE mesh from a PETSc `DMDA` structure
- `ExternalPETScProblem` — `FEProblemBase` subclass that delegates solve calls
  to an external PETSc solver callback
- `PETScDiffusionFDM` — example external solver (2D diffusion using PETSc FDM)
  demonstrating the interface
- `ExternalPetscTimeStepper` — time stepper that synchronizes MOOSE and an
  external PETSc adaptive time integrator

#### When to Use

Use this module when you have an existing PETSc-based finite-difference or
finite-volume code that you want to couple to MOOSE physics through the
MultiApp / Transfer system without rewriting the legacy code as MOOSE kernels.

---

### 4.6 fluid_properties

**Makefile variable**: `FLUID_PROPERTIES := yes`
**Suffix**: `fp`
**Direct dependencies**: `misc`

#### Physics

`fluid_properties` is a pure **equation-of-state library** — it does not
contain flow kernels. It provides a consistent C++ interface for computing
thermodynamic and transport properties of fluids: density, specific enthalpy,
viscosity, thermal conductivity, speed of sound, etc., as functions of state
variables such as pressure and temperature `(p, T)`, or specific volume and
specific internal energy `(v, e)`.

The two primary base interfaces are:
- `SinglePhaseFluidProperties` — properties for a single-phase fluid as a
  function of `(p, T)` or `(v, e)`.
- `TwoPhaseFluidProperties` — extends single-phase with phase equilibrium
  (saturation temperature/pressure, latent heat, surface tension).

All property methods provide both the value and its partial derivatives, so
that AD (automatic differentiation) Jacobians remain exact.

#### Implemented Fluids

| Class | Fluid | Notes |
|---|---|---|
| `IdealGasFluidProperties` | Ideal gas | `p = rho R T` |
| `StiffenedGasFluidProperties` | Stiffened gas | for liquids |
| `Water97FluidProperties` | Water (IAPWS-IF97) | steam tables |
| `CO2FluidProperties` | Carbon dioxide | Span-Wagner EOS |
| `NitrogenFluidProperties` | Nitrogen | Span-Wagner EOS |
| `MethaneFluidProperties` | Methane | Setzmann-Wagner EOS |
| `HydrogenFluidProperties` | Hydrogen | NIST correlations |
| `HeliumFluidProperties` | Helium | NIST correlations |
| `SodiumProperties` | Liquid sodium | ANL correlations |
| `SodiumSaturationFluidProperties` | Sodium (two-phase) | saturation curve |
| `LeadFluidProperties` | Liquid lead | |
| `LeadBismuthFluidProperties` | Lead-bismuth eutectic | |
| `LeadLithiumFluidProperties` | Lead-lithium | fusion blanket |
| `FlibeFluidProperties` | FLiBe molten salt | |
| `FlinakFluidProperties` | FLiNaK molten salt | |
| `NaKFluidProperties` | Sodium-potassium | |
| `BrineFluidProperties` | Brine (H2O + NaCl) | multi-component |
| `SimpleFluidProperties` | Simplified liquid | constant properties |
| `TabulatedBicubicFluidProperties` | Tabulated (bicubic) | from external data |
| `TabulatedFluidProperties` | Tabulated (linear) | from external data |
| `TemperaturePressureFunctionFluidProperties` | User-defined functions | |
| `HEMFluidProperties` | Homogeneous equilibrium mixture | two-phase |
| `CaloricallyImperfectGas` | Real gas (NASA polynomials) | |
| `IdealGasMixtureFluidProperties` | Ideal gas mixture | |

**Two-phase interfaces**:
- `TwoPhaseFluidPropertiesIndependent` — two phases with independent fluid objects
- `StiffenedGasTwoPhaseFluidProperties` — stiffened gas two-phase
- `TwoPhaseNCGFluidProperties` — two-phase with non-condensable gas

#### Key Classes

**Actions**
- `AddFluidPropertiesAction` — reads `[FluidProperties]` block and registers
  fluid property objects

**AuxKernels**
- `FluidDensityAux`, `TemperatureAux`, `PressureAux` — evaluate `rho(p,T)`,
  `T(p,h)`, `p(rho,T)` as field variables
- `SpecificEnthalpyAux` — `h(p,T)`
- `StagnationPressureAux`, `StagnationTemperatureAux` — stagnation conditions
  for compressible flow
- `SaturationTemperatureAux` — saturation temperature from two-phase object

**Materials**
- `FluidPropertiesMaterialPT` — evaluates properties as material properties
  given `p` and `T` field variables
- `FluidPropertiesMaterialVE` — given specific volume `v` and energy `e`
- `SodiumPropertiesMaterial` — sodium-specific material wrapper

**UserObjects**
- `FluidPropertiesInterrogator` — command-line utility to print a table of
  fluid properties for a given fluid and range of state points

**Utilities**
- `DimensionlessFlowNumbers` — Reynolds, Prandtl, Nusselt, etc.
- `BrentsMethod` — Brent root finding for property inversions
- `NewtonInversion` — Newton inversion for property conversions

#### Module Dependencies

`fluid_properties` depends on `misc` for physical constants and the
`GravityVectorInterface`.

#### When to Use

Nearly every flow module depends on `fluid_properties`. Use it directly when
you need to query thermodynamic properties in a material or user object without
running a full flow solve. For example, a heat structure material might call
`_fp.rho_from_p_T(p, T, drho_dp, drho_dT)` to compute density derivatives
needed for an energy equation.

---

### 4.7 fsi (Fluid-Structure Interaction)

**Makefile variable**: `FSI := yes`
**Suffix**: `fsi`
**Direct dependencies**: `navier_stokes`, `solid_mechanics`

#### Physics

The `fsi` module couples a fluid flow (incompressible or weakly compressible
Navier-Stokes) with a structural domain (solid_mechanics). Two primary coupling
mechanisms are provided:

1. **Arbitrary Lagrangian-Eulerian (ALE) moving mesh**: The fluid mesh moves
   with the structural displacement. The `ConvectedMesh` kernel accounts for
   the mesh velocity in the advection term:

   ```
   rho * (du/dt + (u - u_mesh) . grad u)  +  grad p  -  div(mu * grad u)  =  0
   ```

   `ConvectedMeshPSPG` adds the PSPG stabilization contribution for the
   mesh-velocity term.

2. **Acoustic-structural coupling**: For acoustic domains (pressure `p_a`
   satisfying the wave equation `rho d^2 p_a/dt^2 = -K div u_s`), the
   `AcousticInertia` kernel and `StructureAcousticInterface` interface kernel
   couple the acoustic pressure to structural acceleration.

3. **Penalty velocity continuity**: `ADPenaltyVelocityContinuity` and its
   Newmark-Beta variant enforce fluid and structure velocity continuity at the
   interface through a penalty method.

#### Key Classes

**Kernels**
- `ConvectedMesh` — corrects the convection term for mesh velocity
- `ConvectedMeshPSPG` — PSPG stabilization for the mesh-velocity correction
- `AcousticInertia` — `rho * d^2 p_a / dt^2` for pressure-acoustic elements

**Interface Kernels**
- `ADPenaltyVelocityContinuity` — penalizes slip between fluid and structure
- `ADPenaltyVelocityContinuityNewmarkBeta` — Newmark-Beta variant for dynamic
  coupling
- `CoupledPenaltyInterfaceDiffusion` — generic penalty diffusion coupling
  (base class for FSI)
- `StructureAcousticInterface` — acoustic pressure — structural acceleration
  coupling

**Boundary Conditions**
- `FluidFreeSurfaceBC` — kinematic free surface boundary condition for the
  fluid height

**AuxKernels**
- `WaveHeightAuxKernel` — computes free-surface wave height from fluid solution

#### Module Dependencies

`fsi` requires both `navier_stokes` (for flow kernels and variables) and
`solid_mechanics` (for structural kernels and displacement variables).

#### When to Use

Use `fsi` for problems where fluid flow deforms or moves a solid boundary: e.g.,
turbine blade vibration, pipe conveying flow, sloshing in a tank, or acoustic
pressure loading on a structure. For one-way coupling (fluid pressure as a load
on a structure with negligible deformation) you may simply use `solid_mechanics`
with a Pressure BC driven by a flow solution via MultiApp transfer, which avoids
the complexity of ALE mesh motion.

---

### 4.8 functional_expansion_tools

**Makefile variable**: `FUNCTIONAL_EXPANSION_TOOLS := yes`
**Suffix**: `fet`
**Direct dependencies**: none

#### Physics

`functional_expansion_tools` is a numerical coupling library rather than a
physics module. It provides a framework for representing and transferring
field quantities between disparate meshes or applications using orthogonal
functional expansions (series representations).

A field `f(x)` is approximated as

```
f(x) ≈  sum_{n=0}^{N} c_n * phi_n(x)
```

where `phi_n(x)` are orthogonal basis functions and `c_n` are expansion
coefficients computed by integration over a domain. The primary use case is
in MultiApp coupling: one app computes a volumetric or boundary integral of
a field to obtain coefficients, then passes those coefficients to another app
which evaluates the series to reconstruct the field. This avoids interpolation
errors on non-matching meshes.

#### Supported Basis Functions

- **Legendre** polynomials — 1D, for use on line segments or Cartesian
  decompositions
- **Zernike** polynomials — 2D, for circular or radial cross-sections
- **Cartesian** composite series — tensor products of 1D Legendre series in x, y, z
- **CylindricalDuo** — composite of radial Zernike and axial Legendre series

#### Key Classes

**Series**
- `FunctionalBasisInterface` — abstract base for all series types
- `SingleSeriesBasisInterface` — base for 1D series
- `CompositeSeriesBasisInterface` — base for multidimensional composites
- `Legendre` — Legendre polynomial evaluation to order N
- `Zernike` — Zernike polynomial evaluation
- `Cartesian` — 3D Cartesian Legendre series
- `CylindricalDuo` — cylindrical composite series

**Functions**
- `FunctionSeries` — MOOSE Function subclass that evaluates a functional
  expansion; coefficients are mutable (updated from UserObjects or Transfers)
- `MemoizedFunctionInterface` — caches function evaluations for efficiency
- `MutableCoefficientsFunctionInterface` — base for functions whose coefficients
  can be set from outside (Transfers)

**UserObjects**
- `FXVolumeUserObject` — integrates a MOOSE variable over a volume to compute
  expansion coefficients
- `FXBoundaryValueUserObject` — integrates over a boundary (value)
- `FXBoundaryFluxUserObject` — integrates a boundary flux

**AuxKernels**
- `FunctionSeriesToAux` — evaluates a `FunctionSeries` at quadrature points and
  stores as an auxiliary variable

**Boundary Conditions**
- `FXValueBC` — Dirichlet BC set from a `FunctionSeries` evaluation
- `FXFluxBC` — Neumann BC (flux) from a `FunctionSeries`
- `FXValuePenaltyBC` — penalty Dirichlet using the series value

**Transfers**
- `MultiAppFXTransfer` — transfers expansion coefficients between apps via
  the Reporter system, enabling FX-based coupling in MultiApp setups

#### Module Dependencies

None. `functional_expansion_tools` is a leaf module.

#### When to Use

Use `functional_expansion_tools` when you need to couple two multiphysics
applications on non-conforming meshes and want a globally smooth representation
of the coupling field rather than a nodal/element-wise transfer. It is
particularly useful in nuclear applications where neutronics codes (on assembly
meshes) are coupled to thermal-hydraulics codes (on pin-level meshes) via
functional expansions. For direct nodal transfers on matching meshes, use
the standard MOOSE Transfer system instead.

---

### 4.9 geochemistry

**Makefile variable**: `GEOCHEMISTRY := yes`
**Suffix**: `gc`
**Direct dependencies**: none

#### Physics

`geochemistry` provides a complete geochemical modeling capability rooted in
equilibrium thermodynamics. It reads from thermodynamic databases (the
Geochemist's Workbench format is supported) and solves the nonlinear system of
mass-action equations to find equilibrium speciation, mineral saturation,
redox reactions, and gas equilibria. Key capabilities include:

- **Equilibrium speciation**: For each secondary species `j`, the activity
  product must satisfy

  ```
  log10(Q_j / K_j(T))  =  0
  ```

  where `K_j(T)` is the temperature-dependent equilibrium constant and
  `Q_j` is the ion activity product computed from primary species activities.

- **Kinetic reactions**: Transition-state theory rates for mineral
  precipitation/dissolution, redox reactions, and biodegradation.

- **Activity coefficients**: Debye-Huckel B-dot model, Davies model, and
  no-activity-coefficient options.

- **Temperature dependence**: Polynomial or analytical fits for `log10(K(T))`.

- **Reactive transport in space**: The `GeochemistrySpatialReactor` enables
  a spatially distributed geochemical solve, with advection and diffusion of
  total dissolved concentrations coupled to the geochemical equilibrium.

- **Database reading and validation**: `GeochemicalDatabaseReader` parses the
  thermodynamic database; `GeochemicalDatabaseValidator` checks consistency.

#### Key Classes

**Actions**
- `AddTimeDependentReactionSolverAction` — sets up a time-evolving batch reactor
- `AddTimeIndependentReactionSolverAction` — equilibrium batch solve
- `AddSpatialReactionSolverAction` — spatially distributed reactive transport
- `AddGeochemistrySolverAction` — general solver setup
- `AddGeochemicalModelInterrogatorAction` — interactive database interrogation

**UserObjects**
- `GeochemicalModelDefinition` — reads the database and defines the chemical
  system (which species/minerals/gases are in the model)
- `GeochemistryTimeIndependentReactor` — solves a single equilibrium problem
- `GeochemistryTimeDependentReactor` — tracks geochemical evolution over time
  with fluid additions, extractions, and kinetics
- `GeochemistrySpatialReactor` — nodal geochemical reactor for spatially
  distributed problems

**Kernels**
- `GeochemistryTimeDerivative` — time derivative for spatial reactive transport
- `GeochemistryDispersion` — diffusion/dispersion of total concentrations

**AuxKernels**
- `GeochemistryQuantityAux` — extracts quantities (pH, pe, activity, mass,
  saturation index, etc.) from the reactor for output
- `NodalVoidVolumeAux` — computes void volume at nodes for saturation
  calculations

**Utilities (Utils)**
- `GeochemicalDatabaseReader` — parses GWB-format database files
- `GeochemicalSystem` — the core data structure holding all species, minerals,
  gases, and their current concentrations
- `GeochemicalSolver` — Newton solve of the nonlinear speciation equations
- `PertinentGeochemicalSystem` — reduces the full database to only the species
  relevant for the current problem
- `GeochemistryActivityCoefficientsDebyeHuckel` — Debye-Huckel B-dot activity
  model
- `GeochemistryIonicStrength` — ionic strength calculation
- `GeochemistryKineticRateCalculator` — kinetic rate expressions
- `GeochemistrySpeciesSwapper` — manages basis swaps when a primary species
  reaches zero concentration
- `GeochemistryUnitConverter` — converts between mol/kg, g/kg, etc.

**Outputs**
- `GeochemistryConsoleOutput` — prints human-readable speciation tables to the
  console at each time step
- `GeochemicalModelInterrogator` — interactive output for single-point queries

#### Module Dependencies

None. `geochemistry` is a leaf module.

#### Minimal Example (HCl Batch Reactor)

```
[GeochemicalModelDefinition]
  name = definition
  database_file = '../../database/moose_geochemical_database.json'
  basis_species = 'H2O H+ Cl-'
[]

[TimeDependentReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  charge_balance_species = 'Cl-'
  constraint_species = 'H2O H+ Cl-'
  constraint_value   = '  1  1e-8 1e-3'
  constraint_meaning = 'kg_solvent free_molality free_molality'
  temperature = 25
[]
```

#### When to Use

Use `geochemistry` for zero-dimensional or 1D batch and flow-through reactor
problems where complex aqueous speciation matters: water chemistry,
mineral saturation indices, pH buffering, redox chemistry. For spatial reactive
transport, add `GeochemistrySpatialReactor` (which is within `geochemistry`
itself) or couple `geochemistry` with `porous_flow` via operator splitting in a
MultiApp setup.

---

### 4.10 heat_transfer

**Makefile variable**: `HEAT_TRANSFER := yes`
**Suffix**: `ht`
**Direct dependencies**: `ray_tracing`

#### Physics

`heat_transfer` covers all mechanisms of thermal energy transport in solids and
at interfaces:

- **Fourier heat conduction**: The diffusion equation

  ```
  rho * Cp * dT/dt  =  div(k * grad T)  +  Q
  ```

  where `rho` is density, `Cp` specific heat, `k` thermal conductivity, and
  `Q` a volumetric heat source.

- **Convective heat transfer**: Newton's law of cooling at a boundary,

  ```
  q_n  =  h * (T - T_infinity)
  ```

  with heat transfer coefficient `h` and far-field temperature `T_infinity`.

- **Radiative heat transfer**: Gray diffuse radiation between surfaces via
  view factors (computed by the `ray_tracing` module),

  ```
  q_rad  =  sigma_SB * epsilon * ( T^4 - T_infinity^4 )
  ```

  For P1 radiation (optically thick media), the diffusion approximation gives

  ```
  - div( (1 / (3 * kappa)) * grad G )  +  kappa * G  =  4 * kappa * sigma_SB * T^4
  ```

- **Gap heat transfer**: Conduction and radiation across a thin gap between two
  surfaces (contact conductance model).

- **Joule heating**: Source term `Q_J = sigma |grad phi|^2` from electrical
  conduction.

- **Fin effects**: `FinEfficiencyFunctorMaterial` and
  `FinEnhancementFactorFunctorMaterial` for homogenized fin-enhanced surfaces.

#### Key Classes

**Kernels (FE)**
- `HeatConduction` / `ADHeatConduction` — `-div(k * grad T)`
- `HeatConductionTimeDerivative` / `ADHeatConductionTimeDerivative` —
  `rho * Cp * dT/dt`
- `HeatCapacityConductionTimeDerivative` — variant using volumetric heat
  capacity `rho * Cp` as a material property
- `SpecificHeatConductionTimeDerivative` — `rho * Cp * dT/dt` with separate
  density and specific heat properties
- `HeatSource` / `ADMatHeatSource` — volumetric heat source `Q`
- `JouleHeatingSource` / `ADJouleHeatingSource` — `sigma * |grad phi|^2`
- `AnisoHeatConduction` — anisotropic conductivity tensor `k_{ij}`
- `HomogenizedHeatConduction` / `AnisoHomogenizedHeatConduction` —
  homogenized conductivity for periodic microstructures
- `TrussHeatConduction` — 1D truss element conduction (structural heat)

**FV Kernels**
- `FVHeatConductionTimeDerivative` / `FVFunctorHeatConductionTimeDerivative`
- `FVThermalRadiationSourceSink` — P1 radiation source/sink in the energy equation

**Boundary Conditions (FE)**
- `ConvectiveHeatFluxBC` / `ADConvectiveHeatFluxBC` — Newton cooling
- `CoupledConvectiveFlux` / `CoupledConvectiveHeatFluxBC` — h from a coupled
  variable
- `ConvectiveFluxFunction` — h(T) from a function
- `GapHeatTransfer` — gap conductance BC
- `FunctionRadiativeBC` — Stefan-Boltzmann radiation to a prescribed temperature
- `GrayLambertNeumannBC` — gray diffuse radiation Neumann BC
- `InfiniteCylinderRadiativeBC` — enclosure radiation for a cylinder
- `GaussianEnergyFluxBC` — Gaussian beam energy input
- `DirectionalFluxBC` — directional (laser/beam) flux BC
- `HeatConductionBC` — gradient (Fourier) flux BC

**FV Boundary Conditions**
- `FVFunctorConvectiveHeatFluxBC` — FV convective heat flux
- `FVFunctorRadiativeBC` — FV radiative BC
- `FVThermalResistanceBC` — layered thermal resistance BC
- `FVMarshakRadiativeBC` — Marshak radiative BC for P1 equation
- `FVGaussianEnergyFluxBC` — FV Gaussian beam

**Constraints**
- `GapConductanceConstraint` — mortar-based gap conductance
- `ModularGapConductanceConstraint` — modular version with plug-in gap physics
- `ADInterfaceJouleHeatingConstraint` — Joule heating at an electrical interface

**Interface Kernels**
- `ConjugateHeatTransfer` — fluid-solid conjugate heat transfer on an
  internal interface
- `SideSetHeatTransferKernel` — thin layer resistance heat transfer
- `ThinLayerHeatTransfer` — heat transfer through a thin layer

**Materials**
- `HeatConductionMaterial` — simple isotropic `k`, `Cp`, `rho`
- `AnisoHeatConductionMaterial` — anisotropic conductivity tensor
- `GapConductance` / `GapConductanceConstant` — gap conductance models
- `ElectricalConductivity` — `sigma(T)` for Joule heating
- `FunctionPathEllipsoidHeatSource` — moving ellipsoidal heat source
  (e.g., welding arc)

**Functor Materials**
- `ChurchillChuHTCFunctorMaterial` — Churchill-Chu natural convection HTC
- `ConvectionHeatFluxFunctorMaterial` — general convective heat flux
- `CylindricalGapHeatFluxFunctorMaterial` — gap heat flux in cylindrical
  geometry
- `FinEfficiencyFunctorMaterial` — fin efficiency factor
- `FinEnhancementFactorFunctorMaterial` — total fin enhancement factor
- `RadiativeP1DiffusionCoefficientMaterial` — P1 radiation diffusion coefficient

**Physics Actions**
- `HeatConductionCG` — `Physics` block action for FE heat conduction
- `HeatConductionFV` — `Physics` block action for FV heat conduction

**Actions**
- `ThermalContactAction` — sets up gap heat transfer objects from the
  `[ThermalContact]` block
- `MortarGapHeatTransferAction` — mortar-based gap heat transfer setup
- `RadiationTransferAction` — reads the `[RadiativeTransfer]` block and
  creates view-factor and radiation exchange objects

#### Module Dependencies

`heat_transfer` requires `ray_tracing` (used internally for computing radiative
view factors by shooting rays between surfaces).

#### Minimal Example (Transient Heat Conduction with Convective BC)

```
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
[]

[Variables]
  [T]
    initial_condition = 200.0
  []
[]

[Kernels]
  [heat_dt]
    type = ADHeatConductionTimeDerivative
    variable = T
  []
  [heat_cond]
    type = ADHeatConduction
    variable = T
  []
[]

[BCs]
  [convective_right]
    type = ADConvectiveHeatFluxBC
    variable = T
    boundary = right
    T_infinity = 100.0
    heat_transfer_coefficient = 10.0
  []
[]

[Materials]
  [props]
    type = HeatConductionMaterial
    thermal_conductivity = 1.0
    specific_heat = 1.0
  []
  [density]
    type = Density
    density = 1.0
  []
[]

[Executioner]
  type = Transient
  num_steps = 10
  dt = 1.0
  solve_type = NEWTON
[]
```

#### When to Use

Use `heat_transfer` for any problem involving thermal conduction, convective
cooling, radiation, or gap heat transfer. It is a very commonly enabled module
and is required by `solid_properties`, `navier_stokes` (for energy equations),
and `thermal_hydraulics`.

---

### 4.11 level_set

**Makefile variable**: `LEVEL_SET := yes`
**Suffix**: `ls`
**Direct dependencies**: none

#### Physics

The `level_set` module implements the level-set method for tracking moving
interfaces. A signed distance function `phi` is advected with the local
velocity field `u`:

```
dphi/dt  +  u . grad phi  =  0
```

Periodically, `phi` is reinitialized to restore the signed-distance property
`|grad phi| = 1` using the Olsson-Kreiss reinitialization equation:

```
dphi/dtau  +  div( phi(1 - phi) * n_hat )  -  epsilon * div(grad phi)  =  0
```

where `tau` is pseudo-time and `n_hat = grad phi / |grad phi|`.

The module includes SUPG (Streamline Upwind Petrov-Galerkin) stabilization
for the advection equation to suppress spurious oscillations.

#### Key Classes

**Kernels**
- `LevelSetAdvection` — `u . grad phi` advection term
- `LevelSetAdvectionSUPG` — SUPG-stabilized advection
- `LevelSetTimeDerivativeSUPG` — SUPG-modified time derivative
- `LevelSetForcingFunctionSUPG` — SUPG-stabilized forcing
- `LevelSetOlssonReinitialization` — reinitialization equation kernel

**Functions**
- `LevelSetOlssonBubble` — initial condition for a spherical bubble
- `LevelSetOlssonPlane` — initial condition for a planar interface
- `LevelSetOlssonVortex` — velocity field for the Olsson vortex benchmark

**Postprocessors**
- `LevelSetCFLCondition` — computes the CFL-limited time step for advection
- `LevelSetVolume` — computes the volume enclosed by the zero level set

**Problems**
- `LevelSetProblem` — base problem class for level-set simulations
- `LevelSetReinitializationProblem` — problem class for the reinitialization
  sub-problem

**MultiApps / Transfers**
- `LevelSetReinitializationMultiApp` — sub-app that performs reinitialization
  on a periodic schedule
- `LevelSetMeshRefinementTransfer` — transfers the level-set field to drive
  AMR refinement near the interface

**UserObjects**
- `LevelSetOlssonTerminator` — terminates pseudo-time reinitialization when
  steady state is reached

#### Module Dependencies

None. `level_set` is a leaf module.

#### Minimal Example (1D Advection)

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
[]

[Variables]
  [phi]
    family = LAGRANGE
    order = FIRST
  []
[]

[ICs]
  [phi_ic]
    type = FunctionIC
    variable = phi
    function = phi0
  []
[]

[Kernels]
  [phi_dt]
    type = TimeDerivative
    variable = phi
  []
  [phi_advection]
    type = LevelSetAdvectionSUPG
    variable = phi
    velocity = vel
  []
[]

[AuxVariables]
  [vel]
    family = LAGRANGE_VEC
  []
[]

[Executioner]
  type = Transient
  scheme = crank-nicolson
  num_steps = 100
  dt = 0.01
[]
```

#### When to Use

Use `level_set` for interface-tracking problems where you need to follow a
moving phase boundary (liquid-gas, solid-melt, etc.) through a fixed Eulerian
mesh. For diffuse interface problems (phase-field) use `phase_field` instead.
`level_set` is often used with `navier_stokes` for two-phase flow with sharp
interfaces.

#### Quickstart Example

**Case 73** (`quickstart-runs/case73-level-set-bubble`) demonstrates
SUPG-stabilized level-set advection using the solid-body rotation benchmark:
a circular bubble is advected by a divergence-free rotating velocity field and
completes one full revolution in unit time. `LevelSetAdvectionSUPG` and
`LevelSetTimeDerivativeSUPG` provide the stabilized kernels; `LevelSetCFLCondition`
drives the adaptive time step; `LevelSetVolume` verifies approximate mass
conservation over the rotation.

---

### 4.12 misc

**Makefile variable**: `MISC := yes`
**Suffix**: `misc`
**Direct dependencies**: none

#### Physics

`misc` is a small utility module providing objects and interfaces that do not
belong naturally to any single physics module but are needed by many. It does
not define any particular governing equations.

#### Key Classes

**Kernels**
- `CoefDiffusion` — diffusion with a constant coefficient: `- c * div(grad u)`
- `ThermoDiffusion` / `ADThermoDiffusion` — Soret (thermophoresis) term:
  `- D_T / T * div(T * grad c)`, the mass flux driven by a temperature gradient

**Materials**
- `ArrheniusMaterialProperty` — Arrhenius temperature dependence:
  `A * exp(-E_a / (R * T))`
- `Density` — constant or temperature-dependent density material property

**AuxKernels**
- `CoupledDirectionalMeshHeightInterpolation` — height interpolation along a
  mesh direction for gravitational head computations

**Postprocessors**
- `InternalVolume` — computes the enclosed volume bounded by a surface
- `GeneralSensorPostprocessor` — generalized sensor model postprocessor
- `ThermocoupleSensorPostprocessor` — thermocouple sensor with response lag

**Interfaces**
- `GravityVectorInterface` — provides a consistent `gravity` input parameter
  across objects that need the gravitational acceleration vector

**Utilities**
- `PhysicalConstants` — fundamental constants (Boltzmann, Avogadro, R, etc.)
  as C++ constexpr values

#### Module Dependencies

None. `misc` is a leaf module. It is depended upon by `fluid_properties`,
`thermal_hydraulics`, `scalar_transport`, and indirectly by all modules that
use those.

---

### 4.13 navier_stokes

**Makefile variable**: `NAVIER_STOKES := yes`
**Suffix**: `ns`
**Direct dependencies**: `fluid_properties`, `rdg`, `heat_transfer`

#### Physics

`navier_stokes` is one of the largest MOOSE modules. It provides kernels and
actions for both compressible and incompressible flow in multiple formulations:

**1. Incompressible Navier-Stokes (FE, INSAD)**

Solves the incompressible continuity and momentum equations,

```
div u = 0

rho * (du/dt + u . grad u) + grad p - div(mu * (grad u + (grad u)^T)) = f
```

using the INSAD (Incompressible NS via Automatic Differentiation) kernels with
SUPG/PSPG stabilization.

**2. Incompressible Navier-Stokes (FV, INSFV)**

Finite-volume formulation using cell-centered pressure-velocity coupling with
the SIMPLE or PIMPLE algorithm. Supports the `Physics/NavierStokes/Flow` and
`Physics/NavierStokes/Turbulence` action blocks. The PIMPLE executioner
handles the pressure-velocity segregated loop.

**3. Weakly Compressible Navier-Stokes (WCNS, FV)**

Density varies with temperature but not pressure (Boussinesq-like or low-Mach).

**4. Compressible Navier-Stokes / Euler Equations (CNS)**

Conservative form using primitive or conservative variables `(rho, rho*u, rho*E)`.

```
d(rho)/dt  +  div(rho u)  =  0

d(rho u)/dt  +  div(rho u ⊗ u + p I)  -  div(tau)  =  f

d(rho E)/dt  +  div((rho E + p) u)  -  div(tau . u)  -  div(k grad T)  =  Q
```

The compressible solver is used with `rdg` for stabilization on structured
or unstructured meshes.

**5. Turbulence Models**

- Mixing-length model via `INSFVMixingLengthTurbulentViscosityAux`
- k-epsilon model via `kEpsilonViscosityAux`
- RANS y+ wall functions

**6. Porous Media Flow (Volumetric Averaging)**

The `Physics/NavierStokes/Flow` block with porosity supports volume-averaged
(porous media) incompressible flow.

**7. HDG (Hybridized Discontinuous Galerkin)**

`NavierStokesLHDGVelocityDirichletBC`, `NavierStokesLHDGOutflowBC`, and
related HDG kernels implement a high-order hybridized DG formulation for
the Stokes and Navier-Stokes equations.

#### Key Classes

**Actions**
- `INSAction` / `INSFVAction` — set up incompressible NS (FE and FV)
- `CNSAction` — set up compressible NS
- Physics block: `Physics/NavierStokes/Flow`, `Physics/NavierStokes/Turbulence`,
  `Physics/NavierStokes/ScalarTransport`, `Physics/NavierStokes/Energy`

**FV Kernels (selected)**
- `INSFVMomentumAdvection` — convection in momentum equation
- `INSFVMomentumDiffusion` — viscous diffusion
- `INSFVMomentumPressure` — pressure gradient
- `INSFVMomentumBodyForce` — body force (gravity)
- `INSFVMomentumGravity` — gravity contribution
- `INSFVMassAdvection` — continuity equation advection
- `PINSFVMomentumAdvection` — porous media momentum convection
- `WCNSFVMomentumTimeDerivative` — weakly compressible time derivative
- `NSFVPhaseChangeMomentumFlux` — phase change contribution to momentum

**FE Kernels (selected)**
- `INSADMomentumAdvection` — INSAD convection
- `INSADMomentumViscous` — INSAD viscous diffusion
- `INSADMomentumPressure` — INSAD pressure gradient
- `INSADMass` — INSAD continuity
- `INSADMomentumSUPGCorrection` — SUPG stabilization

**Compressible (NS) Kernels**
- `NSMassInviscidFlux` — inviscid mass flux
- `NSMomentumInviscidFlux` — inviscid momentum flux
- `NSEnergyInviscidFlux` — inviscid energy flux
- `NSMomentumViscousFlux` — viscous momentum flux
- `NSEnergyViscousFlux` — viscous energy flux

**Executioners**
- `PIMPLE` — pressure-implicit with splitting of operators, segregated NS loop
- `LinearAssemblySegregatedSolve` — linear assembly solver for PIMPLE

**AuxKernels**
- `NSMachAux`, `NSTemperatureAux`, `NSPressureAux`, `NSVelocityAux` —
  derived flow quantities for compressible flow
- `INSFVMixingLengthTurbulentViscosityAux` — mixing-length turbulent viscosity
- `RANSYPlusAux`, `WallFunctionWallShearStressAux`, `WallFunctionYPlusAux` —
  turbulent wall functions

#### Module Dependencies

`navier_stokes` depends on:
- `fluid_properties` — equation of state
- `rdg` — slope reconstruction and limiting for compressible flows
- `heat_transfer` — energy equation coupling

#### Minimal Example (Lid-Driven Cavity, FE INSAD)

```
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 32
    ny = 32
    elem_type = QUAD9
  []
[]

[Modules/IncompressibleNavierStokes]
  equation_type = transient
  velocity_boundary = 'bottom right top             left'
  velocity_function = '0 0    0 0   lid_function 0  0 0'
  pressure_pinned_node = 0
  density_name = rho
  dynamic_viscosity_name = mu
  use_ad = true
  laplace = true
  family = LAGRANGE
  order = SECOND
[]

[Materials]
  [const]
    type = ADGenericConstantMaterial
    prop_names  = 'rho  mu'
    prop_values = '1.0  0.01'
  []
[]

[Functions]
  [lid_function]
    type = ParsedFunction
    expression = '4*x*(1-x)'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  num_steps = 20
  dt = 0.5
[]
```

#### When to Use

Use `navier_stokes` for single-phase fluid flow problems: incompressible
isothermal or non-isothermal flow (conjugate heat transfer), weakly
compressible flow, fully compressible / supersonic flows, and turbulent RANS
simulations. For multi-phase porous media flow use `porous_flow`. For
1D network thermal-hydraulics use `thermal_hydraulics`. For nuclear subchannel
flow use `subchannel`.

---

### 4.14 optimization

**Makefile variable**: `OPTIMIZATION := yes`
**Suffix**: `opt`
**Direct dependencies**: none

#### Physics

`optimization` provides design optimization and PDE-constrained inverse problem
capabilities. Two main problem classes are supported:

1. **PDE-constrained inverse problems**: Given noisy measurements of a physical
   field, find source terms, material properties, or boundary conditions that
   best reproduce the observations. The gradient of the objective function with
   respect to parameters is computed using the **adjoint method**: one forward
   solve and one adjoint solve yield the gradient at the cost of two PDE solves,
   independent of the number of parameters.

   ```
   min_{theta} J(u, theta) = 1/2 sum_i (u(x_i) - d_i)^2 + regularization
   subject to: L(u; theta) = 0
   ```

   The adjoint equation is `L^*(lambda) = dJ/du`, and the gradient is
   `dJ/dtheta = integral lambda * dL/dtheta`.

2. **Topology optimization**: Density-based SIMP (Solid Isotropic Material with
   Penalization) method for structural or thermal optimization. Material density
   `rho_e` is the design variable and stiffness/conductivity is penalized:
   `E_e = E_min + rho_e^p * (E_0 - E_min)`.

External optimizers are accessed through PETSc's TAO library. Supported
solvers include L-BFGS, BFGS, conjugate gradient (BNCG), and bound-constrained
variants.

#### Key Classes

**Executioners**
- `Optimize` — top-level executioner that drives the TAO optimization loop
- `OptimizeSolve` — inner solve manager (forward + adjoint cycle per iteration)
- `AdjointSolve` — steady adjoint solve
- `AdjointTransientSolve` — transient adjoint solve (reverse in time)
- `SteadyAndAdjoint` — combined steady + adjoint in one executioner
- `TransientAndAdjoint` — combined transient + adjoint

**OptimizationReporters**
- `OptimizationReporter` — base class: manages parameter vector, objective,
  and gradient communication with TAO
- `GeneralOptimization` — general-purpose optimization reporter
- `ParameterMeshOptimization` — parameters defined on a mesh (distributed
  field inversion)

**Functions**
- `OptimizationFunction` — base class for functions whose parameters are
  design variables
- `NearestReporterCoordinatesFunction` — maps parameter values to spatial
  locations (for source reconstruction)
- `ParameterMeshFunction` — spatially distributed parameter field on a mesh
- `ParsedOptimizationFunction` — parsed expression where named coefficients
  are optimization parameters

**DiracKernels**
- `ReporterTimePointSource` — point source whose magnitude comes from an
  optimization reporter (used for source term identification)

**VectorPostprocessors (inner product evaluations for gradient)**
- `ElementOptimizationSourceFunctionInnerProduct` — `int lambda * dQ/dtheta dV`
- `ElementOptimizationDiffusionCoefFunctionInnerProduct` — material property
  inner product
- `ElementOptimizationFunctionInnerProduct` — general function inner product
- `SideOptimizationFunctionInnerProduct` — boundary inner product
- `AdjointStrainSymmetricStressGradInnerProduct` — structural optimization
  gradient

**Materials**
- `CostSensitivity` — topology optimization density sensitivity
- `ReporterOffsetFunctionMaterial` — material whose property is shifted by
  optimization reporter values

**UserObjects**
- `DensityUpdate` — SIMP density update for topology optimization
- `DensityUpdateTwoConstraints` — density update with two constraints
- `SensitivityFilter` — density filter for topology optimization
- `AdjointSolutionUserObject` — stores and retrieves adjoint solution fields

#### Module Dependencies

None. `optimization` is a leaf module. For topology optimization of structures,
also enable `solid_mechanics`; for thermal topology optimization, enable
`heat_transfer`.

#### Minimal Example (Parameter Identification)

```
[Optimization]
[]

[OptimizationReporter]
  type = GeneralOptimization
  parameter_names = 'q_source'
  num_values = 1
  initial_condition = '100'
  lower_bounds = '0'
  upper_bounds = '1000'
[]

[Executioner]
  type = Optimize
  tao_solver = taobncg
  petsc_options_iname = '-tao_gatol'
  petsc_options_value  = '1e-6'
[]
```

#### When to Use

Use `optimization` for **inverse problems** (find sources/properties from
measurements), **design optimization** (maximize stiffness or minimize
temperature), and **topology optimization** (find optimal material layout). For
purely stochastic uncertainty quantification without optimization, use
`stochastic_tools` instead.

---

### 4.15 peridynamics

**Makefile variable**: `PERIDYNAMICS := yes`
**Suffix**: `pd`
**Direct dependencies**: `solid_mechanics`

#### Physics

Peridynamics is a nonlocal continuum theory that replaces the classical
differential operators with integral operators, making it naturally suited for
modeling fracture and damage without special treatment at crack surfaces. The
equation of motion is

```
rho * d^2 u(x,t)/dt^2  =  integral_{H(x)} f(u', u, x', x) dV'  +  b(x,t)
```

where `H(x)` is the **horizon** (neighborhood) of point `x`, `f` is the
pairwise force function, and the integral is over all material points within
the horizon. The horizon radius `delta` is a length scale parameter.

Two main formulations are implemented:

1. **Bond-based peridynamics (BPD)**: The simplest model, with pairwise forces
   along the bond direction. Poisson's ratio is fixed at 1/4 (2D) or 1/3 (3D)
   for isotropic materials.

2. **Non-ordinary state-based peridynamics (NOSPDPD)**: A more general model
   that allows arbitrary constitutive relations (Poisson's ratio is free).
   The force state is computed from the deformation state, and the bond-level
   force depends on deformation of all bonds in the neighborhood.

3. **Ordinary state-based peridynamics (OSPD)**: Intermediate formulation.

The module also supports:
- Small and finite strain formulations
- Generalized plane strain
- Thermal conduction in the peridynamic framework via BPD heat conduction
- Failure / damage via bond breaking (stretch-based or stress-based criteria)

#### Key Classes

**Actions**
- `MechanicsActionPD` — sets up peridynamic mechanics kernels from a high-level
  block
- `GeneralizedPlaneStrainActionPD` — sets up GPS for peridynamic models

**Kernels**
- `MechanicsBPD` — bond-based mechanics kernel
- `MechanicsOSPD` — ordinary state-based mechanics kernel
- `ForceStabilizedSmallStrainMechanicsNOSPD` — force-stabilized NOS small strain
- `HorizonStabilizedFormISmallStrainMechanicsNOSPD` — horizon-stabilized NOS
  (Form I, small strain)
- `HorizonStabilizedFormIISmallStrainMechanicsNOSPD` — horizon-stabilized NOS
  (Form II, small strain)
- `HorizonStabilizedFormIFiniteStrainMechanicsNOSPD` — finite strain variant
- `HorizonStabilizedFormIIFiniteStrainMechanicsNOSPD` — finite strain variant
- `WeakPlaneStressNOSPD` — weak plane-stress constraint for NOS
- `GeneralizedPlaneStrainOffDiagNOSPD` — off-diagonal GPS for NOS
- `GeneralizedPlaneStrainOffDiagOSPD` — off-diagonal GPS for OS
- `HeatConductionBPD` — bond-based peridynamic heat conduction
- `HeatSourceBPD` — heat source for BPD thermal model

**Materials**
- `ComputeSmallStrainConstantHorizonMaterialBPD` — BPD small strain, constant
  horizon
- `ComputeSmallStrainVariableHorizonMaterialBPD` — BPD, variable horizon
- `ComputeSmallStrainConstantHorizonMaterialOSPD` — OSPD small strain
- `ComputeSmallStrainVariableHorizonMaterialOSPD` — OSPD variable horizon
- `ComputeSmallStrainNOSPD` — NOSPD small strain
- `ComputeFiniteStrainNOSPD` — NOSPD finite strain
- `ComputePlaneSmallStrainNOSPD` — plane strain NOS
- `ComputePlaneFiniteStrainNOSPD` — plane finite strain NOS

**AuxKernels**
- `BondStatusBasePD` — base for bond status (intact / broken)
- `StretchBasedFailureCriterionPD` — bond failure by stretch threshold
- `RankTwoBasedFailureCriteriaNOSPD` — NOS failure by rank-two stress/strain
- `NodalRankTwoPD` — nodal rank-two tensor output
- `NodalVolumePD` — peridynamic nodal volume output
- `BoundaryOffsetPD` — boundary correction for partial neighborhoods

#### Module Dependencies

`peridynamics` depends on `solid_mechanics` for constitutive model classes
(elasticity tensors, strain and stress materials from solid_mechanics can be
reused in NOS peridynamics).

#### When to Use

Use `peridynamics` when you need to simulate fracture, fragmentation, or
damage evolution without explicitly tracking crack surfaces — particularly when
crack paths are not known a priori. For fatigue crack growth along a known
path, `xfem` or cohesive zone models within `solid_mechanics` may be more
appropriate.

---

### 4.16 phase_field

**Makefile variable**: `PHASE_FIELD := yes`
**Suffix**: `pf`
**Direct dependencies**: `solid_mechanics`

#### Physics

The `phase_field` module implements diffuse-interface models for microstructure
evolution. The two main equation types are:

**Allen-Cahn equation** (non-conserved order parameter `eta`):

```
tau * d(eta)/dt  =  - dF/d(eta)  +  kappa * laplacian(eta)
```

where `F` is a free energy density and `kappa` is the gradient energy coefficient.

**Cahn-Hilliard equation** (conserved composition `c`):

```
dc/dt  =  div( M * grad( dF/dc - kappa * laplacian(c) ) )
```

where `M` is the mobility. The chemical potential `mu = dF/dc - kappa * laplacian(c)`.

The free energy `F` is typically provided as a `DerivativeParsedMaterial`
so that MOOSE's AD system can compute all required derivatives automatically.

The module supports:
- Single and multi-phase systems
- Polycrystal grain growth (N-order-parameter model) via
  `PolycrystalKernelAction` and `GrainGrowthAction`
- KKS (Kim-Steinbach-Kobayashi-Suzuki) multi-phase model
- Phase-field crystal (PFC) models
- Coupling to mechanics via eigenstrain (coherency strains)
- Discrete nucleation via `DiscreteNucleationInserter`
- EBSD data reader for initializing realistic microstructures

#### Key Classes

**Actions**
- `GrainGrowthAction` — sets up N-order-parameter grain growth with the
  `Modules/PhaseField/GrainGrowth` shorthand
- `ConservedAction` — sets up Cahn-Hilliard with `Modules/PhaseField/Conserved`
- `NonconservedAction` — sets up Allen-Cahn with `Modules/PhaseField/Nonconserved`
- `KKSAction` — Kim-Steinbach-Kobayashi-Suzuki multi-phase
- `GrandPotentialKernelAction` — grand potential model setup
- `PolycrystalKernelAction` / `PolycrystalVariablesAction` — polycrystal setup
- `PolycrystalColoringICAction` — IC for multi-grain Voronoi tessellation

**Kernels**
- `AllenCahn` / `ADAllenCahn` — Allen-Cahn bulk driving force `- dF/d(eta)`
- `ACInterface` / `ADACInterface` — Allen-Cahn gradient term `kappa * laplacian(eta)`
- `CahnHilliard` / `CahnHilliardAniso` — Cahn-Hilliard `div(M grad mu)`
- `SplitCHWRes` / `SplitCHMath` — split (Ladyzhenskaya) form CH
- `CHInterface` — CH gradient energy term
- `ACGrGr` — grain-grain interaction driving force for polycrystal
- `ACKKSMultiACBulkF`, `KKSSplitCHCRes`, etc. — KKS model kernels
- `PFFractureBulkRate` — phase-field fracture model

**Materials**
- `DerivativeParsedMaterial` — free energy from a parsed expression; automatic
  differentiation provides `dF/d(eta)`, `d^2F/d(eta)^2`, etc.
- `DerivativeSumMaterial` — sums multiple free energy contributions
- `MathEBFreeEnergy` — double-well free energy `F = (1-c)^2 (1+c)^2`
- `GrandPotentialInterface` — grand potential interfacial contributions
- `PolycrystalElasticDrivingForce` — elastic driving force for grain rotation

**ICs (Initial Conditions)**
- `SmoothCircleIC`, `LatticeSmoothCircleIC`, `PolycrystalCircleGrainIC` —
  circular grain/particle ICs
- `CrossIC`, `ClosePackIC` — two-phase initial conditions
- `EBSDReaderPointDataAux` / `EBSDReaderAvgDataAux` — EBSD microstructure data

**UserObjects**
- `EBSDReader` — reads EBSD data files (TSL OIM format)
- `PolycrystalVoronoi` — Voronoi tessellation for polycrystal IC
- `DiscreteNucleationInserter` — probabilistic nucleation engine
- `FeatureFloodCount` — connected-component labeling for grain tracking

#### Module Dependencies

`phase_field` depends on `solid_mechanics` for coupled mechano-chemical
eigenstrain contributions and stress-assisted phase transformations.

#### Minimal Example (Cahn-Hilliard Spinodal Decomposition)

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 32
  ny = 32
  xmax = 50
  ymax = 50
[]

[Modules/PhaseField/Conserved/cv]
  solve_type = direct
  free_energy = F
  kappa = 2.0
  mobility = 1.0
[]

[ICs]
  [c_ic]
    type = CrossIC
    variable = cv
    x1 = 5.0  y1 = 5.0
    x2 = 45.0 y2 = 45.0
  []
[]

[Materials]
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = F
    coupled_variables = 'cv'
    expression = '(1-cv)^2 * (1+cv)^2'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  dt = 0.7
  num_steps = 50
[]
```

#### When to Use

Use `phase_field` for microstructure evolution problems: spinodal decomposition,
grain growth, solidification, precipitation, phase transformations, and
phase-field fracture. It is the correct choice whenever the interface is
diffuse (width ~ diffuse length scale >> mesh spacing) and you want to avoid
explicit interface tracking. For sharp interfaces use `level_set` or `xfem`.

---

### 4.17 porous_flow

**Makefile variable**: `POROUS_FLOW := yes`
**Suffix**: `pflow`
**Direct dependencies**: `solid_mechanics`, `fluid_properties`, `chemical_reactions`

#### Physics

`porous_flow` is a comprehensive module for multiphase, multicomponent flow
and transport in deformable porous media. It solves coupled THM
(Thermo-Hydro-Mechanical) problems. The governing equations include:

**Mass balance** for fluid component `kappa` in phase `beta`:

```
d/dt ( phi * sum_beta S_beta * rho_beta * X_beta^kappa )
    + div( sum_beta rho_beta * X_beta^kappa * ( - k * k_rel_beta / mu_beta * ( grad p_beta - rho_beta * g ) ) )
    + div( sum_beta rho_beta * D_beta^kappa * grad X_beta^kappa )
    =  q^kappa
```

**Energy balance** (optional):

```
d/dt ( (1 - phi) * rho_rock * Cp_rock * T + phi * sum_beta S_beta * rho_beta * h_beta )
    + div( sum_beta rho_beta * h_beta * q_fluid_beta )
    - div( lambda * grad T )
    =  Q
```

**Mechanical coupling** (optional, via `solid_mechanics`): Biot consolidation,

```
div( sigma - alpha_B * p_f * I )  =  0
```

where `alpha_B` is the Biot coefficient, and the porosity evolves with
volumetric strain.

The Dictator pattern: a `PorousFlowDictator` UserObject is used by all
other PorousFlow objects to coordinate phase structure, number of components,
and variable naming — eliminating any duplication of this information.

#### Key Classes

**Actions (high-level)**
- `PorousFlowBasicTHM` — single-phase fully saturated THM problem
- `PorousFlowFullySaturated` — fully saturated single-phase flow
- `PorousFlowUnsaturated` — unsaturated single-phase (variably saturated)
- `PorousFlowSinglePhaseBase` — base for single-phase actions

**Kernels**
- `PorousFlowMassTimeDerivative` — time rate of change of fluid mass
- `PorousFlowAdvectiveFlux` — Darcy advective flux
- `PorousFlowDispersiveFlux` — diffusion and dispersion
- `PorousFlowFullySaturatedDarcyFlow` — fully saturated Darcy
- `PorousFlowHeatConduction` — Fourier heat conduction
- `PorousFlowHeatAdvection` — advective heat transport
- `PorousFlowEnergyTimeDerivative` — energy time derivative
- `PorousFlowEffectiveStressCoupling` — Biot effective stress coupling
- `PorousFlowHeatMassTransfer` — heat/mass exchange between matrix and fracture
- `FluxLimitedTVDAdvection` — TVD-limited advection for sharp fronts

**DiracKernels (Wells)**
- `PorousFlowPeacemanBorehole` — Peaceman well model for production/injection
- `PorousFlowPolyLineSink` — polyline sink along a well trajectory
- `PorousFlowLineSink` — single line sink
- `PorousFlowPointSourceFromPostprocessor` — flow rate from a postprocessor
- `PorousFlowPointEnthalpySourceFromPostprocessor` — enthalpy source for wells

**Boundary Conditions**
- `PorousFlowSink` — general sink/source BC
- `PorousFlowPiecewiseLinearSink` — flux varies piecewise linearly with
  pressure
- `PorousFlowHalfGaussianSink` / `PorousFlowHalfCubicSink` — smooth flux
  functions
- `PorousFlowEnthalpySink` — sink with enthalpy for thermal BCs
- `PorousFlowOutflowBC` — free outflow

**Materials** (many): capillary pressure, relative permeability, porosity,
permeability, fluid density and viscosity wrappers, etc.

**FV Kernels**
- `FVPorousFlowAdvectiveFlux`, `FVPorousFlowDispersiveFlux`,
  `FVPorousFlowMassTimeDerivative`, `FVPorousFlowHeatAdvection`,
  `FVPorousFlowHeatConduction`, `FVPorousFlowEnergyTimeDerivative`

#### Module Dependencies

`porous_flow` requires:
- `solid_mechanics` — mechanical deformation coupling (Biot)
- `fluid_properties` — equation of state
- `chemical_reactions` — reactive transport

#### Minimal Example (Single-Phase Darcy Flow)

```
[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 1e7
  []
[]

[AuxVariables]
  [temperature]
    initial_condition = 293
  []
[]

[PorousFlowBasicTHM]
  porepressure = porepressure
  temperature = temperature
  coupling_type = Hydro
  gravity = '0 0 -9.81'
  fp = simple_fluid
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2e9
    viscosity = 1e-3
    density0 = 1000
  []
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0  0 1e-13 0  0 0 1e-13'
  []
[]

[Executioner]
  type = Transient
  dt = 100
  num_steps = 10
  solve_type = NEWTON
[]
```

#### When to Use

Use `porous_flow` for subsurface flow applications: groundwater hydrology, CO2
geological storage, oil/gas reservoir simulation, geothermal energy extraction,
nuclear waste repository modeling, and hydro-mechanically coupled problems
(consolidation, land subsidence). For simpler single-phase unsaturated flow
without chemistry, the older `richards` module is also available but `porous_flow`
is more feature-complete and actively developed.

---

### 4.18 ray_tracing

**Makefile variable**: `RAY_TRACING := yes`
**Suffix**: `ray`
**Direct dependencies**: none

#### Physics

`ray_tracing` provides infrastructure for geometric ray tracing through a MOOSE
mesh. Rays are traced element by element, intersecting element faces using
robust geometric algorithms. This enables:

- Line integrals of field quantities along arbitrary ray paths
- View-factor computation for radiation heat transfer
- Neutron/particle transport line-of-sight calculations
- Laser beam path integration in metal additive manufacturing

A ray carries user-defined data (scalars or arrays) that is accumulated or
modified as the ray traverses each element. Ray boundary conditions control
what happens when a ray hits a boundary: it can be killed, reflected,
or trigger a user-defined action.

#### Key Classes

**Base**
- `RayTracingObject` — base class for all ray tracing objects

**RayKernels** (applied to each element a ray crosses)
- `IntegralRayKernel` — accumulates an integral along the ray
- `MaterialIntegralRayKernel` — integrates a material property
- `VariableIntegralRayKernel` — integrates a field variable
- `FunctionIntegralRayKernel` — integrates a function
- `LineSourceRayKernel` — deposits a source term proportional to path length
  (e.g., volumetric laser heating)
- `AuxRayKernel` — applies an aux variable contribution along the ray
- `RayDistanceAux` — stores the distance traveled by the ray at each element
- `KillRayKernel` — terminates the ray on a condition

**RayBCs** (applied when a ray hits a boundary)
- `KillRayBC` — absorb the ray (terminate)
- `ReflectRayBC` — specular reflection
- `GeneralRayBC` — user-defined BC logic

**Postprocessors**
- `RayIntegralValue` — reports the accumulated integral on a named ray
- `RayDataValue` — reports a specific data field on a ray
- `RayTracingStudyResult` — reports study-level statistics

**Study Classes**
- `ParallelRayStudy` — manages parallel execution of ray tracing across
  distributed meshes, handling ray communication between MPI ranks

**Outputs**
- `RayTracingExodus` — visualizes ray paths in Exodus format
- `RayTracingMeshOutput` / `RayTracingNemesis` — alternative ray path outputs

#### Module Dependencies

None. `ray_tracing` is a leaf module. It is used internally by `heat_transfer`
for view-factor calculations.

#### When to Use

Use `ray_tracing` as a dependency when you need geometric ray tracing in your
physics, e.g., for radiation view factors or laser path integration. For users,
this module is usually enabled indirectly through `heat_transfer`.

---

### 4.19 rdg (Reconstructed Discontinuous Galerkin)

**Makefile variable**: `RDG := yes`
**Suffix**: `rdg`
**Direct dependencies**: none

#### Physics

`rdg` provides slope reconstruction and flux limiting infrastructure for
cell-centered finite-volume / DG methods applied to hyperbolic systems of
conservation laws. It implements the AEFV (Average Effective Flux Value)
approach, where a piecewise linear reconstruction of the solution is computed
within each cell and limited to prevent oscillations (slope limiting), then
an upwind or approximate Riemann solver is used to compute inter-cell fluxes.

The module does not contain flow physics — it provides the numerical building
blocks that `navier_stokes` and `thermal_hydraulics` use for their
compressible-flow solvers.

#### Key Classes

**UserObjects (Reconstruction and Limiting)**
- `SlopeReconstructionOneD` / `SlopeReconstructionMultiD` — cell-slope
  reconstruction from neighbors
- `SlopeLimitingBase` / `SlopeLimitingOneD` / `SlopeLimitingMultiDBase` —
  TVD slope limiters (minmod, van Leer, etc.)
- `AEFVSlopeLimitingOneD` — 1D AEFV-specific slope limiting

**UserObjects (Flux Functions)**
- `InternalSideFluxBase` — abstract internal face flux
- `AEFVUpwindInternalSideFlux` — first-order upwind internal flux
- `BoundaryFluxBase` / `BCUserObject` — boundary flux objects
- `AEFVFreeOutflowBoundaryFlux` — free outflow (zero-gradient) boundary flux
- `RDGFluxBase` — base for all RDG flux computations

**DG Kernels**
- `AEFVKernel` — applies the inter-cell flux as a DG residual contribution

**Boundary Conditions**
- `AEFVBC` — applies boundary flux in the DG formulation

**Materials**
- `AEFVMaterial` — reconstructed state material (provides extrapolated values
  from the limited reconstruction for flux evaluation)

**Postprocessors**
- `RDGBoundaryFluxPostprocessor` — integrates a boundary flux

#### Module Dependencies

None. `rdg` is a leaf module. It is used by `navier_stokes` and
`thermal_hydraulics`.

#### When to Use

`rdg` is almost never used directly by end users. It is a numerical
infrastructure module that is automatically enabled when `navier_stokes` or
`thermal_hydraulics` is enabled. Use it directly only if you are implementing
a new hyperbolic conservation law solver that needs slope reconstruction and
limiting on an unstructured mesh.

---

### 4.20 reactor

**Makefile variable**: `REACTOR := yes`
**Suffix**: `rct`
**Direct dependencies**: none

#### Physics

`reactor` is a **mesh generation** module for nuclear reactor geometry. It does
not contain physics objects (kernels, materials, BCs). Instead, it provides
specialized `MeshGenerator` objects that build pin-level, assembly-level, and
core-level meshes for typical reactor geometries (hexagonal or square lattices,
concentric circles, control drums, etc.).

#### Key Classes

**MeshGenerators**
- `PinMeshGenerator` — builds a fuel pin mesh (concentric annular regions)
- `AssemblyMeshGenerator` — arranges pins in a lattice to build an assembly
- `CoreMeshGenerator` — arranges assemblies into a core layout
- `PolygonConcentricCircleMeshGenerator` — concentric circles inside a polygon
  (hex or square) with user-controlled region IDs and mesh refinement
- `HexagonConcentricCircleAdaptiveBoundaryMeshGenerator` — adaptive boundary
  hexagon for smooth hex-hex interfaces
- `CartesianConcentricCircleAdaptiveBoundaryMeshGenerator` — Cartesian variant
- `PatternedHexMeshGenerator` / `PatternedCartesianMeshGenerator` — arranges
  sub-meshes in a hex or Cartesian pattern
- `PatternedHexPeripheralModifier` / `PatternedCartesianPeripheralModifier` —
  modifies the periphery of a patterned mesh
- `SimpleHexagonGenerator` — simple hexagonal mesh
- `HexIDPatternedMeshGenerator` / `CartesianIDPatternedMeshGenerator` — pattern
  with assigned IDs
- `HexagonMeshTrimmer` / `CartesianMeshTrimmer` — trim reactor meshes to
  reflect symmetry planes
- `ControlDrumMeshGenerator` — control drum geometry
- `PeripheralRingMeshGenerator` / `PeripheralTriangleMeshGenerator` —
  radial peripheral regions
- `RevolveGenerator` — revolves a 2D mesh to create 3D geometry
- `ReactorMeshParams` — shared mesh parameter object (pitch, ring counts, etc.)
- `DepletionIDGenerator` — assigns depletion zone IDs for burnup calculations
- `ExtraElementIDCopyGenerator` — propagates element IDs from one mesh to
  another after spatial mapping
- `SubdomainExtraElementIDGenerator` — assigns extra element IDs by subdomain
- `CoarseMeshExtraElementIDGenerator` — maps IDs from a coarse mesh
- `GapLineMeshGenerator` — creates 1D line meshes in gaps between assemblies
- `AzimuthalBlockSplitGenerator` — azimuthally subdivides a hex assembly
- `FlexiblePatternGenerator` — pattern mesh with flexible arrangement
- `TriPinHexAssemblyGenerator` — three-pin hexagonal assembly layout
- `AdvancedConcentricCircleGenerator` — advanced concentric circle mesh

**Meshing Divisions**
- `HexagonalGridDivision` — maps spatial positions to hex grid bins

**Positions**
- `HexagonalGridPositions` — generates positions of pins/assemblies on a
  hexagonal grid
- `CartesianGridPositions` — Cartesian grid positions

**Functions**
- `MultiControlDrumFunction` — evaluates control drum position as a function
  of time (for rotation/insertion models)

#### Module Dependencies

None. `reactor` is a leaf module. It is depended upon by `subchannel`.

#### When to Use

Use `reactor` whenever you need to build reactor-geometry meshes programmatically
without CAD tools. The mesh generators handle pin lattice patterns (square PWR,
hexagonal SFR/MSR), automatic element ID assignment for multi-physics coupling,
and geometric symmetry trimming. The resulting meshes are standard MOOSE meshes
compatible with all physics modules.

---

### 4.21 richards

**Makefile variable**: `RICHARDS := yes`
**Suffix**: `rich`
**Direct dependencies**: none

#### Physics

`richards` implements the Richards equation for variably-saturated subsurface
flow (unsaturated zone):

```
d(phi * rho * S) / dt  +  div( rho * k * k_rel(S) / mu * (grad p - rho g) )  =  q
```

where `S` is saturation, `k_rel(S)` is the relative permeability (a function
of saturation), and `p` is pore pressure. Saturation and capillary pressure
are related by the van Genuchten, Brooks-Corey, or other models.

The module also supports:
- **Two-phase (Q2P) formulation**: Gas and liquid phases with independent
  pressure variables
- **Desorption** from a solid matrix (methane/CO2 in coal seams)
- **Peaceman borehole** models for injection/production wells

Note: `richards` is an older module and is partially superseded by
`porous_flow`, which is more general and better maintained. `porous_flow`
supports Richards-type problems through its unsaturated single-phase action.

#### Key Classes

**Actions**
- `Q2PAction` — sets up the two-phase Q2P problem

**Kernels**
- `DarcyFlux` — Darcy advective flux `div(rho * k * k_rel / mu * (grad p - rho g))`
- `PoroFullSatTimeDerivative` — fully saturated time derivative
- `Q2PNodalMass` / `Q2PNegativeNodalMassOld` — Q2P mass balance kernels

**DiracKernels**
- `PeacemanBorehole` — Peaceman well model
- `RichardsBorehole` — Richards-specific borehole
- `Q2PBorehole` — two-phase borehole
- `RichardsPolyLineSink` — polyline sink

**Boundary Conditions**
- `RichardsExcav` — excavation BC (suddenly changing pressure)
- `RichardsPiecewiseLinearSink` / `RichardsHalfGaussianSink` — flux BCs
- `Q2PPiecewiseLinearSink` — Q2P flux BC

**AuxKernels**
- `RichardsSatAux`, `RichardsDensityAux`, `RichardsRelPermAux` — saturation,
  density, relative permeability output
- `DarcyFluxComponent` — Darcy flux vector component

**UserObjects (constitutive relations)**
- `RichardsDensity*` — density models (ideal gas, constant, linear)
- `RichardsRelPerm*` — relative permeability (van Genuchten, Brooks-Corey,
  Corey, power law)
- `RichardsCapPres*` — capillary pressure models
- `RichardsSat*` / `RichardsSeff*` — saturation / effective saturation models

#### Module Dependencies

None. `richards` is a leaf module.

#### When to Use

Use `richards` for simple single-phase or two-phase unsaturated porous flow
problems. For new development with multi-component, multi-phase flow or
hydromechanical coupling, prefer `porous_flow`.

---

### 4.22 scalar_transport

**Makefile variable**: `SCALAR_TRANSPORT := yes`
**Suffix**: `st`
**Direct dependencies**: `chemical_reactions`, `navier_stokes`, `thermal_hydraulics`,
`fluid_properties`, `heat_transfer`, `rdg`, `ray_tracing`, `solid_properties`, `misc`

#### Physics

`scalar_transport` extends the `chemical_reactions` framework with Lagrange
multiplier (LM) discretizations suitable for tightly coupled advection-
diffusion-reaction in the context of `thermal_hydraulics` and `navier_stokes`
flows. It also provides multi-species diffusion physics.

The LM kernels implement the governing equation

```
d(rho * c) / dt  +  div(rho * u * c - rho * D * grad c)  +  R(c)  =  f
```

in a Lagrange multiplier (discontinuous) sense, which is advantageous for
hyperbolic advection-dominated transport.

#### Key Classes

**Kernels (LM)**
- `LMKernel` — base Lagrange multiplier kernel
- `LMTimeKernel` — `dc/dt` in LM form
- `LMDiffusion` — diffusion in LM form
- `TimeDerivativeLM` — time derivative for LM variables
- `BodyForceLM` — body force in LM form
- `CoupledForceLM` — coupling force in LM form

**Boundary Conditions**
- `BinaryRecombinationBC` — flux BC for surface recombination:
  `J = K_r * c_1 * c_2` (e.g., hydrogen isotope recombination)
- `DissociationFluxBC` — dissociation flux at a surface

**Physics**
- `MultiSpeciesDiffusionCG` / `MultiSpeciesDiffusionPhysicsBase` — multi-species
  diffusion via the `Physics` block system

#### Module Dependencies

`scalar_transport` has the widest dependency tree of any module, requiring
`chemical_reactions`, `navier_stokes`, `thermal_hydraulics`, `fluid_properties`,
`heat_transfer`, `rdg`, `ray_tracing`, `solid_properties`, and `misc`.

#### When to Use

Use `scalar_transport` for species transport problems tightly coupled to
thermal hydraulics flows defined by `thermal_hydraulics`, or for multi-species
diffusion (hydrogen permeation, tritium transport). For simpler advection-
diffusion problems not coupled to `thermal_hydraulics`, use `chemical_reactions`
directly or standard MOOSE `Convection` / `Diffusion` kernels.

---

### 4.23 solid_mechanics (formerly tensor_mechanics)

**Makefile variable**: `SOLID_MECHANICS := yes`
**Suffix**: `sm`
**Direct dependencies**: none

#### Physics

`solid_mechanics` (renamed from `tensor_mechanics` in 2023) is the main
structural mechanics module. It covers a broad range of continuum solid
mechanics:

- **Small strain linear elasticity**:
  ```
  epsilon  =  (grad u + (grad u)^T) / 2
  sigma    =  C : epsilon
  div sigma = -f
  ```

- **Finite strain nonlinear mechanics** (including large deformations):
  Multiplicative decomposition `F = F_e * F_p`, with full geometric nonlinearity
  via updated Lagrangian or total Lagrangian formulations.

- **Inelastic material models**:
  - Isotropic and kinematic hardening plasticity (J2, Drucker-Prager)
  - Power-law creep and Norton creep
  - Damage models (Mazars, isotropic scalar damage)
  - Combined creep-plasticity

- **Eigenstrain contributions**: Thermal expansion, swelling, transformation
  strain, etc. can be added as eigenstrains `epsilon^* = alpha * Delta T * I`.

- **Structural elements**: Beam elements (Timoshenko and Euler-Bernoulli),
  shell elements, truss elements.

- **Dynamics**: Newmark-Beta and HHT time integration for structural dynamics
  (including seismic analysis).

- **Fracture mechanics**: J-integral, stress intensity factors (KI, KII, KIII),
  crack tip enrichment (via XFEM), cohesive zone models.

- **Generalized plane strain** and **plane stress** reduced-order models.

- **Cosserat (micropolar) elasticity** for granular media.

The primary action block is:

```
[Physics/SolidMechanics/QuasiStatic]
```

which automates kernel, material, and variable setup.

#### Key Classes (selected)

**Actions**
- `TensorMechanicsAction` (alias for `QuasiStaticSolidMechanicsPhysics`) —
  main action for quasi-static solid mechanics
- `CommonSolidMechanicsAction` — shared parameters across blocks
- `DomainIntegralAction` — fracture domain integral
- `CohesiveZoneAction` — cohesive zone model setup
- `LineElementAction` — beam/truss element setup
- `GeneralizedPlaneStrainAction` — GPS constraints
- `CavityPressureAction` — internal pressurized cavity
- `GlobalStrainAction` — periodic RVE global strain constraints

**Kernels**
- `StressDivergenceTensors` / `ADStressDivergenceTensors` — `div sigma + f = 0`
- `InertialForce` / `ADInertialForce` — `rho * d^2 u / dt^2`
- `DynamicStressDivergenceTensors` — dynamics with Rayleigh damping
- `StressDivergenceRZTensors` — axisymmetric

**Materials (Kinematics)**
- `ComputeSmallStrain` / `ADComputeSmallStrain` — infinitesimal strain tensor
- `ComputeIncrementalSmallStrain` — incremental small strain
- `ComputeFiniteStrain` / `ADComputeFiniteStrain` — finite strain (incremental)
- `ComputeAxisymmetricRZSmallStrain` — axisymmetric small strain
- `ComputePlaneSmallStrain` — plane strain
- `ComputeEigenstrainFromTemperature` — thermal expansion eigenstrain

**Materials (Elasticity Tensors)**
- `ComputeIsotropicElasticityTensor` / `ADComputeIsotropicElasticityTensor` —
  isotropic `C(E, nu)`
- `ComputeElasticityTensor` — general 6×6 anisotropic tensor
- `ComputeLayeredCosseratElasticityTensor` — Cosserat
- `ComputeRotatedElasticityTensor` — rotated crystal or fiber directions

**Materials (Stress)**
- `ComputeLinearElasticStress` / `ADComputeLinearElasticStress` — `sigma = C:epsilon`
- `ComputeFiniteStrainElasticStress` / `ADComputeFiniteStrainElasticStress`
- `ComputeMultipleInelasticStress` — general inelastic stress with multiple
  return-mapping models
- `ADComputeMultipleInelasticStress`
- `ComputeDamageStress` — continuum damage mechanics

**AuxKernels**
- `RankTwoAux` — extracts components of stress, strain, or any rank-2 tensor
- `RankTwoScalarAux` — scalar invariants (von Mises stress, hydrostatic, etc.)
- `ElasticEnergyAux` — elastic strain energy density
- `NewmarkAccelAux` / `NewmarkVelAux` — time-integration acceleration/velocity

#### Module Dependencies

None. `solid_mechanics` is a leaf module and one of the most important
dependencies in the ecosystem.

#### Minimal Example (3D Finite Strain Elasticity)

```
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4  ny = 4  nz = 4
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    strain = FINITE
    generate_output = 'stress_xx stress_yy vonmises_stress'
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e10
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []
[]

[BCs]
  [fix_bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_left]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_back]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [pull_top]
    type = ADDirichletBC
    variable = disp_y
    boundary = top
    value = 0.01
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
```

#### When to Use

Use `solid_mechanics` for any solid structural analysis: linear elastic statics
or dynamics, plasticity, creep, thermal stress, fracture, composites, shells,
beams, or coupled thermo-mechanical problems. It is the central module depended
upon by `contact`, `fsi`, `porous_flow`, `peridynamics`, `phase_field`, and
`xfem`.

---

### 4.24 solid_properties

**Makefile variable**: `SOLID_PROPERTIES := yes`
**Suffix**: `sp`
**Direct dependencies**: `heat_transfer`

#### Physics

`solid_properties` is a material property library for solid-phase thermal
properties, analogous to `fluid_properties` for fluids. It provides
temperature-dependent specific heat capacity `Cp(T)`, density `rho(T)`, and
thermal conductivity `k(T)` for engineering materials commonly found in
nuclear and energy systems.

#### Implemented Materials

| Class | Material | Notes |
|---|---|---|
| `ThermalGraphiteProperties` | Nuclear graphite | SGL Group data |
| `ThermalMonolithicSiCProperties` | Monolithic SiC | |
| `ThermalCompositeSiCProperties` | SiC/SiC composite | |
| `ThermalUCProperties` | Uranium carbide (UC) | |
| `ThermalSS316Properties` | 316 stainless steel | |
| `ThermalFunctionSolidProperties` | User-defined functions | |
| `ThermalSolidProperties` | Abstract base class | |

#### Key Classes

**Actions**
- `AddSolidPropertiesAction` — reads `[SolidProperties]` block and registers
  solid property objects

**Materials**
- `ThermalSolidPropertiesMaterial` — evaluates `k(T)`, `Cp(T)`, `rho(T)` as
  material properties from a `ThermalSolidProperties` object
- `ConstantDensityThermalSolidPropertiesMaterial` — variant with constant `rho`
- `TungstenThermalPropertiesMaterial` — tungsten-specific material (PFC armor)

**Functor Materials**
- `ThermalSolidPropertiesFunctorMaterial` — functor version for FV problems

**Postprocessors**
- `ThermalSolidPropertiesPostprocessor` — evaluates `k`, `Cp`, or `rho` at a
  given temperature for verification

#### Module Dependencies

`solid_properties` depends on `heat_transfer` because it uses `HeatConductionMaterial`
conventions and the `HeatTransferApp` registry.

#### When to Use

Use `solid_properties` when you need validated, literature-based thermal
properties for nuclear or energy materials rather than hard-coding constants.
For custom materials with arbitrary functional forms, use
`ThermalFunctionSolidProperties` which wraps three MOOSE `Function` objects
for `k(T)`, `Cp(T)`, and `rho(T)`.

---

### 4.25 stochastic_tools

**Makefile variable**: `STOCHASTIC_TOOLS := yes`
**Suffix**: `st`
**Direct dependencies**: none

#### Physics

`stochastic_tools` provides uncertainty quantification (UQ), sensitivity
analysis (SA), surrogate modeling, and Bayesian calibration capabilities.
It operates in a MultiApp framework: a master application runs the sampling
and analysis while sub-applications (the actual physics simulations) are
launched for each sample point.

Key capabilities:

- **Monte Carlo and quasi-Monte Carlo sampling**: Sobol sequences, Latin
  hypercube sampling, Cartesian product grids
- **Polynomial Chaos Expansion (PCE)**: Stochastic collocation or regression
- **Gaussian Process (GP) surrogate models** with covariance kernels
  (squared exponential, Matérn, LMC multi-output)
- **Proper Orthogonal Decomposition (POD)** plus GP for field-level surrogates
- **Active learning / Bayesian optimization**: acquisition functions for
  efficient experimental design
- **Markov Chain Monte Carlo (MCMC)**: for Bayesian posterior sampling
- **Sensitivity analysis**: Sobol indices (first-order, total effect)
- **Statistical postprocessors**: mean, standard deviation, percentiles,
  confidence intervals, CVaR, etc.

#### Key Classes

**Distributions**
- `Normal`, `TruncatedNormal`, `Lognormal` — Gaussian-based distributions
- `Uniform`, `Beta`, `Gamma`, `Weibull` — standard statistical distributions
- `StudentT`, `FDistribution` — t and F distributions
- `Logistic`, `JohnsonSB` — specialized distributions
- `KernelDensity1D` — non-parametric kernel density estimate

**Samplers** (not shown in headers but registered via actions)
- Monte Carlo sampler, Latin hypercube sampler, Sobol sampler,
  Cartesian product sampler, adaptive sampler, MCMC sampler

**Surrogates**
- `PolynomialChaos` — PCE surrogate
- `GaussianProcess` — GP regression surrogate
- `NearestPointSurrogate` — nearest-neighbor interpolation
- `LibtorchANNSurrogate` — neural network surrogate (requires LibTorch)

**Covariance Functions**
- `SquaredExponentialCovariance` — `exp(-r^2 / (2 l^2))`
- `ExponentialCovariance` — `exp(-r / l)`
- `MaternHalfIntCovariance` — Matérn 1/2, 3/2, 5/2
- `LMC` — Linear Model of Coregionalization for multi-output GP

**Acquisition Functions** (for Bayesian optimization / active learning)
- `ExpectedImprovement` — maximizes expected reduction of objective
- `ProbabilityofImprovement` — probability of improving over current best
- `UpperConfidenceBound` — UCB exploration-exploitation
- `CoefficientOfVariation` — maximizes prediction uncertainty
- `UFunction` — U-function for reliability analysis
- `BayesianPosteriorTargeted` — targeted Bayesian sampling

**Controls**
- `MultiAppSamplerControl` — passes sampler values to sub-app parameters
- `SamplerReceiver` — receives control values in the sub-app

**Actions**
- `ParameterStudyAction` — high-level parameter study (replaces manual setup
  for common sampling patterns)
- `StochasticToolsAction` — enables the stochastic tools framework

#### Module Dependencies

None. `stochastic_tools` is a leaf module.

#### Minimal Example (Monte Carlo Parameter Study)

```
[StochasticTools]
[]

[Distributions]
  [uniform_E]
    type = Uniform
    lower_bound = 1e10
    upper_bound = 3e10
  []
[]

[Samplers]
  [mc]
    type = MonteCarlo
    distributions = 'uniform_E'
    num_rows = 100
    seed = 42
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = sub.i
  []
[]

[Transfers]
  [param]
    type = SamplerParameterTransfer
    to_multi_app = sub
    sampler = mc
    parameters = 'Materials/elasticity/youngs_modulus'
  []
  [result]
    type = SamplerReporterTransfer
    from_multi_app = sub
    sampler = mc
    from_reporter = 'max_disp/value'
    to_reporter = 'results/max_disp'
  []
[]

[Reporters]
  [results]
    type = StochasticReporter
  []
[]

[Executioner]
  type = Steady
[]
```

#### When to Use

Use `stochastic_tools` for UQ, SA, Bayesian calibration, or surrogate
modeling workflows. It is well-suited to problems where the physics is
encapsulated in a MOOSE sub-application that can be launched repeatedly.
For deterministic design optimization and inverse problems, use `optimization`
instead.

---

### 4.26 subchannel

**Makefile variable**: `SUBCHANNEL := yes`
**Suffix**: `sc`
**Direct dependencies**: `fluid_properties`, `heat_transfer`, `reactor`

#### Physics

`subchannel` implements subchannel thermal-hydraulics for nuclear fuel rod
bundles. Subchannel analysis is a method that decomposes the rod bundle cross-
section into hydraulic subchannels (the interstitial regions between rods),
and solves simplified 1D-per-subchannel mass, momentum, and energy equations
coupled by lateral crossflow between adjacent subchannels.

Two lattice geometries are supported:
- **Quadrilateral (square lattice)**: PWR-type rod bundles
- **Triangular (hexagonal lattice)**: SFR/CANDU-type bundles

The governing equations per subchannel are:

```
Mass:    d(rho A)/dt  +  d(rho A V)/dz  =  sum_j w_j  (crossflow)

Momentum: d(rho A V)/dt  +  d(rho A V^2)/dz  +  A * dp/dz
          + F_friction + F_spacer  =  sum_j ( w_j * V_j )

Energy:  d(rho A h)/dt  +  d(rho A V h)/dz
         =  q_rod * P_heated / A  +  sum_j w_j h_j
```

where `V` is axial velocity, `h` is specific enthalpy, `w_j` is crossflow
between adjacent subchannels `j`, and `q_rod` is the rod power.

Closure correlations provide friction factors, crossflow resistance, and
heat transfer coefficients (Dittus-Boelter, Gnielinski, etc.).

#### Key Classes

**Mesh Classes**
- `QuadSubChannelMesh` — square lattice subchannel mesh (1D axial + 2D map)
- `TriSubChannelMesh` — triangular lattice subchannel mesh
- `SubChannelMesh` — base class

**Mesh Generators**
- `SCMQuadSubChannelMeshGenerator` — builds a square lattice mesh
- `SCMQuadPinMeshGenerator` — builds the pin mesh for a square lattice
- `SCMTriSubChannelMeshGenerator` — builds a triangular lattice mesh
- `SCMTriPinMeshGenerator` — builds the pin mesh for a triangular lattice
- `SCMDetailedQuadSubChannelMeshGenerator` — detailed (2D cross-section) quad
- `SCMDetailedTriSubChannelMeshGenerator` — detailed triangular
- `SCMQuadDuctMeshGenerator` — assembly duct mesh

**Actions**
- `QuadSubChannelBuildMeshAction` — high-level quad mesh builder
- `TriSubChannelBuildMeshAction` — high-level tri mesh builder
- `SubChannelAddVariablesAction` — registers subchannel variables
  (mass flux, enthalpy, pressure, velocity)
- `SubChannelCreateProblemAction` — creates the subchannel problem
- `AddSCMClosureAction` — adds closure relations

**AuxKernels**
- `SCMMassFlowRateAux` — mass flow rate per subchannel
- `SCMFlatMassFlowRateAux` — flattened mass flow rate
- `SCMQuadPowerAux` / `SCMTriPowerAux` — rod power contributions

**Initial Conditions**
- `SCMQuadFlowAreaIC` / `SCMTriFlowAreaIC` — flow area from geometry
- `SCMQuadWettedPerimIC` / `SCMTriWettedPerimIC` — wetted perimeter
- `SCMMassFlowRateIC` — initial mass flow rate
- `SCMQuadPowerIC` / `SCMTriPowerIC` — initial power distribution

**Fluid Properties**
- `PBSodiumFluidProperties` — sodium properties (ANL correlation) for
  liquid metal fast reactor subchannel analysis

#### Module Dependencies

`subchannel` requires:
- `fluid_properties` — coolant equation of state
- `heat_transfer` — rod-to-coolant heat transfer
- `reactor` — hexagonal grid positions and pin mesh building

#### When to Use

Use `subchannel` for nuclear fuel assembly thermal-hydraulics analysis where
a full 3D CFD (Navier-Stokes) solution is too expensive. It gives rapid axial
temperature and enthalpy profiles for licensing and design studies. For 3D
CFD of coolant in rod bundles, use `navier_stokes` with the `reactor`
geometries.

---

### 4.27 tensor_mechanics (Alias)

`tensor_mechanics` is **not a separate module directory**. It was the original
name of the solid mechanics module. The directory does not exist; instead,
the build system provides a compatibility alias:

```makefile
# In modules/modules.mk:
ifeq ($(TENSOR_MECHANICS),yes)
  SOLID_MECHANICS := yes
  $(warning The tensor mechanics module was renamed to the solid mechanics \
    module. Please update your Makefile and replace TENSOR_MECHANICS with \
    SOLID_MECHANICS)
endif
```

In source code, `TensorMechanicsApp` is a typedef:

```cpp
// solid_mechanics/include/base/TensorMechanicsApp.h
#include "SolidMechanicsApp.h"
typedef SolidMechanicsApp TensorMechanicsApp;
```

Similarly, `TensorMechanicsAction` is a typedef for `QuasiStaticSolidMechanicsPhysics`.

**Action required**: Update any Makefile with `TENSOR_MECHANICS := yes` to
`SOLID_MECHANICS := yes`, and any `include "TensorMechanicsApp.h"` directives to
`include "SolidMechanicsApp.h"` in downstream application code. Input files
using `[Modules/TensorMechanics/...]` blocks still parse correctly because
the action aliases are maintained for backward compatibility.

---

### 4.28 thermal_hydraulics

**Makefile variable**: `THERMAL_HYDRAULICS := yes`
**Suffix**: `th`
**Direct dependencies**: `navier_stokes`, `fluid_properties`, `heat_transfer`, `rdg`,
`ray_tracing`, `solid_properties`, `misc`

#### Physics

`thermal_hydraulics` (THM) solves 1D single-phase compressible flow through
pipe networks. The governing equations are the 1D compressible Euler equations
with wall friction and heat transfer source terms:

```
d(rho A)/dt  +  d(rho A V)/dx  =  0

d(rho A V)/dt  +  d((rho V^2 + p) A)/dx  =  p * dA/dx  -  f * rho A V|V| / (2 D_h)
                                             +  rho A g_x

d(rho A E)/dt  +  d((rho E + p) A V)/dx  =  q_wall * P_htf  +  rho A V g_x
```

where `f` is the Darcy-Weisbach friction factor, `D_h` is hydraulic diameter,
`q_wall` is wall heat flux, and `P_htf` is the heated perimeter.

The module uses a **component-based input syntax**: the user defines a network
of components (`FlowChannel1Phase`, `HeatStructureCylindrical`,
`Inlet*`, `Outlet*`, junctions, pumps, turbines, etc.) and the framework
automatically builds the mesh, variables, kernels, and boundary conditions.

#### Key Components

**Flow Components**
- `FlowChannel1Phase` — a straight pipe segment
- `FlowChannelGasMix` — gas mixture flow channel
- `ElbowPipe1Phase` — pipe with directional change
- `Junction*`, `VolumeJunction1Phase` — pipe junctions

**Boundary Components**
- `InletStagnationPressureTemperature1Phase` — inlet stagnation conditions
- `InletVelocityTemperature1Phase` — inlet velocity + temperature
- `InletMassFlowRate1Phase` — inlet mass flow rate
- `InletDensityVelocity1Phase` — inlet density and velocity
- `Outlet1Phase` — pressure outlet
- `SolidWall1Phase` — solid wall (no flow) boundary
- `FreeBoundary1Phase` — free outflow

**Heat Structure Components**
- `HeatStructureCylindrical` — cylindrical heat structure (fuel rod wrapper)
- `HeatStructurePlate` — flat plate heat structure
- `HeatStructureFromFile3D` — 3D heat structure from an Exodus file

**Heat Transfer Components**
- `HeatTransferFromHeatStructure1Phase` — 1D-1D coupling to a heat structure
- `HeatTransferFromSpecifiedTemperature1Phase` — prescribed wall temperature
- `HeatTransferFromHeatFlux1Phase` — prescribed wall heat flux
- `HeatTransferFromExternalApp1Phase` — coupling via MultiApp transfer

**Rotating Machine Components**
- `ShaftConnectedPump1Phase` — centrifugal pump on a shaft
- `ShaftConnectedTurbine1Phase` — turbine on a shaft
- `ShaftConnectedCompressor1Phase` — compressor on a shaft
- `Shaft` — rotating shaft connecting machines
- `ShaftConnectedMotor` — electric motor
- `SimpleTurbine1Phase` — simplified turbine model

**Well Components**
- `InjectionWell` / `ProductionWell` — geothermal/petroleum well models

**Closures** (friction and heat transfer correlations)
- `Closures1PhaseSimple` — simple friction and HTC correlations
- `Closures1PhaseVoigtDittusBoelter` — Dittus-Boelter heat transfer

#### Key Classes (selected)

**Base**
- `FlowModel` / `FlowModel1PhaseBase` — flow model definition
- `THMActionComponent` — base for all THM components

**Actions**
- `AddComponentAction` — dispatches component creation
- `THMBuildMeshAction` — builds the network mesh from components
- `CoupledHeatTransferAction` — sets up fluid-structure coupling

**AuxKernels**
- `MachNumberAux` — Mach number `M = V / c`
- `ReynoldsNumberAux` — Reynolds number
- `PrandtlNumberAux` — Prandtl number
- `ADConvectiveHeatFlux1PhaseAux` — convective heat flux

#### Module Dependencies

`thermal_hydraulics` depends on:
- `navier_stokes` — compressible flow infrastructure
- `fluid_properties` — equation of state
- `heat_transfer` — heat structure physics
- `rdg` — slope reconstruction for compressible flow
- `ray_tracing` — radiation in heat structures
- `solid_properties` — heat structure material properties
- `misc` — utility classes

#### Minimal Example (Simple Pipe Flow)

```
[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [pipe]
    type = FlowChannel1Phase
    orientation = '1 0 0'
    position = '0 0 0'
    length = 1.0
    n_elems = 10
    A = 0.01
    initial_T = 300
    initial_p = 1e5
    initial_vel = 1.0
    fp = fp
    closures = simple_closures
    f = 0.01
  []
  [inlet]
    type = InletVelocityTemperature1Phase
    input = 'pipe:in'
    vel = 1.0
    T = 300
  []
  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 1e5
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  dt = 0.001
  num_steps = 100
[]
```

#### When to Use

Use `thermal_hydraulics` for 1D pipe network simulations of nuclear reactor
primary and secondary loops, steam generators, feedwater systems, and
geothermal wells. For higher-dimensional flow (2D/3D CFD), use `navier_stokes`.
For multi-phase porous media flow, use `porous_flow`.

#### Quickstart Example

**Case 72** (`quickstart-runs/case72-thm-pipe-flow`) demonstrates the Component
DSL by modeling a single horizontal pipe with a velocity inlet, pressure outlet,
and constant Darcy-Weisbach wall friction using `IdealGasFluidProperties` and
`Closures1PhaseSimple`. Postprocessors report steady-state pressure drop, Mach
number, and wall friction force for comparison with the analytic Darcy-Weisbach
formula.

---

### 4.29 xfem

**Makefile variable**: `XFEM := yes`
**Suffix**: `xfem`
**Direct dependencies**: `solid_mechanics`

#### Physics

`xfem` (Extended Finite Element Method) implements crack simulation by
enriching the standard finite element approximation with discontinuous
and near-crack-tip singular basis functions. Unlike a standard FE mesh, XFEM
does not require remeshing as the crack grows — instead, the mesh is cut
geometrically by a crack surface description.

The displacement field near a crack is decomposed as:

```
u(x) = sum_i N_i(x) u_i                   (standard FE part)
      + sum_j N_j(x) H(x) a_j             (Heaviside enrichment, cut elements)
      + sum_k N_k(x) sum_l phi_l(x) b_kl  (crack-tip enrichment functions)
```

where `H(x)` is the Heaviside function (±1 on either side of the crack) and
`phi_l` are the four Westergaard near-tip functions capturing the `r^{1/2}`
singularity.

Key capabilities:
- 2D and 3D crack propagation
- Stationary cracks with accurate SIF computation (via domain integrals from
  `solid_mechanics`)
- Crack growth driven by Paris law, stress corrosion cracking, or user-defined
  rates
- Level-set based crack representation
- Cohesive zone interface model with XFEM

#### Key Classes

**Base**
- `XFEM` — the main XFEM infrastructure object that manages crack cut data,
  element subdivision, and enrichment DOFs
- `XFEMCutElem` / `XFEMCutElem2D` / `XFEMCutElem3D` — element cutting geometry
- `XFEMFuncs` — geometric utility functions

**UserObjects (Crack Geometry)**
- `CircleCutUserObject` — circular (penny-shaped) crack
- `CrackMeshCut3DUserObject` — 3D crack surface defined by a mesh
- `ComboCutUserObject` — combination of multiple cut geometries
- `GeometricCut2DUserObject` (base) — 2D crack segment
- `GeometricCutUserObject` (base) — general cut geometry interface
- `CutElementSubdomainModifier` — marks cut element subdomains for XFEM

**Actions**
- `XFEMAction` — reads the `[XFEM]` block and sets up enrichment variables,
  the XFEM object, and crack growth reporters

**Kernels**
- `CrackTipEnrichmentStressDivergenceTensors` — stress divergence with
  near-tip enrichment DOFs

**Materials**
- `ComputeCrackTipEnrichmentSmallStrain` — strain with crack-tip enrichment
- `LevelSetBiMaterialReal` / `LevelSetBiMaterialRankTwo` /
  `LevelSetBiMaterialRankFour` — bi-material properties switching across the
  crack (useful for bi-material interface cracks)
- `XFEMCutSwitchingMaterial` — material that switches value across the cut

**Constraints**
- `XFEMSingleVariableConstraint` — Nitsche-type constraint across the cut
- `XFEMEqualValueAtInterface` — equality constraint at the XFEM interface

**AuxKernels**
- `XFEMVolFracAux` — element volume fraction on each side of the cut
- `XFEMMarkerAux` — marks cut vs. non-cut elements
- `XFEMCutPlaneAux` — cut plane normal and intercept
- `CutSubdomainIDAux` — subdomain ID based on which side of the cut
- `MeshCutLevelSetAux` — level set value defining the crack surface

**Boundary Conditions**
- `CrackTipEnrichmentCutOffBC` — cuts off enrichment functions outside the
  active zone

**Reporters (Crack Growth)**
- `ParisLaw` — Paris-law fatigue crack growth `da/dN = C * DeltaK^m`
- `StressCorrosionCrackingExponential` — SCC exponential growth law

**Outputs**
- `XFEMCutMeshOutput` — outputs the crack mesh for visualization

#### Module Dependencies

`xfem` depends on `solid_mechanics` for displacement variables and stress
computation.

#### Minimal Example (Stationary Crack with Heaviside Enrichment)

```
[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [crack_cut]
    type = CircleCutUserObject
    cut_data = '0.5 0.5 0   0.5 0.0 0   0.5 0.5 0.25'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    strain = SMALL
  []
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
```

#### When to Use

Use `xfem` for fracture mechanics problems where crack paths are not known in
advance, cracks need to propagate through a fixed mesh, or you need to capture
the near-tip singularity accurately without a highly refined mesh around the
crack tip. For predefined stationary cracks, standard `solid_mechanics` with
a conforming mesh and domain integrals may be simpler. For diffuse damage
(not a sharp crack), use `phase_field` fracture models.

#### Quickstart Example

**Case 71** (`quickstart-runs/case71-xfem-heat-crack`) demonstrates XFEM
applied to heat conduction across a stationary insulating crack. A
`LineSegmentCutUserObject` defines the crack geometry; the `[XFEM]` block
activates Heaviside enrichment for cut elements; `XFEMMarkerAux` and
`XFEMVolFracAux` provide visualization of which elements are cut. The
temperature jump across the crack face is clearly visible in the Exodus output.

---

## 5. Cross-Module Coupling Patterns

### 5.1 Thermo-Mechanical Coupling

Couple `heat_transfer` and `solid_mechanics` via thermal eigenstrains:

```
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    strain = SMALL
    eigenstrain_names = 'thermal_eigenstrain'
  []
[]

[Kernels]
  [heat_cond]
    type = ADHeatConduction
    variable = T
  []
  [heat_dt]
    type = ADHeatConductionTimeDerivative
    variable = T
  []
[]

[Materials]
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 1.2e-5
    temperature = T
    eigenstrain_name = 'thermal_eigenstrain'
  []
[]
```

### 5.2 Fluid-Structure Coupling via MultiApp

For loosely coupled FSI (fluid provides pressure on solid surface, solid
provides displacement to fluid mesh):

**Master app** (solid): uses `solid_mechanics` + `MeshDisplacedProblem`
**Sub-app** (fluid): uses `navier_stokes`

Transfers:
- `MultiAppMeshFunctionTransfer` — sends fluid pressure to solid surface BC
- `MultiAppMeshFunctionTransfer` — sends solid displacement to fluid mesh
  ALE (via `ConvectedMesh` kernels)

### 5.3 Porous Flow + Geochemistry via Operator Splitting

`porous_flow` solves transport; `geochemistry` solves equilibrium at each
time step. Operator splitting via MultiApp:

```
[MultiApps]
  [geochemistry_sub]
    type = TransientMultiApp
    input_files = geochemistry.i
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [concentrations_to_geochemistry]
    type = MultiAppMeshFunctionTransfer
    to_multi_app = geochemistry_sub
    source_variable = 'c_Ca c_CO3'
    variable = 'transported_Ca transported_CO3'
  []
  [minerals_from_geochemistry]
    type = MultiAppMeshFunctionTransfer
    from_multi_app = geochemistry_sub
    source_variable = 'calcite_vol_frac'
    variable = 'calcite_vol_frac'
  []
[]
```

### 5.4 Stochastic UQ with Any Physics Sub-App

`stochastic_tools` wraps any MOOSE sub-application:

```
[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = physics.i   # any physics module
  []
[]

[Transfers]
  [param_transfer]
    type = SamplerParameterTransfer
    to_multi_app = sub
    sampler = mc
    parameters = 'Materials/steel/youngs_modulus'
  []
[]
```

### 5.5 Reactor Geometry + Neutronics + THM

A typical nuclear multi-physics workflow:

1. **Mesh**: `reactor` module builds pin/assembly/core mesh
2. **Neutronics**: external app (Griffin, OpenMC) provides power distribution
3. **Heat transfer**: `heat_transfer` + `solid_properties` for fuel and
   cladding temperature
4. **Subchannel**: `subchannel` for coolant temperature and pressure
5. **Mechanics**: `solid_mechanics` for fuel rod deformation (optional)
6. **Coupling**: `stochastic_tools` for uncertainty propagation across the chain

### 5.6 Contact + Thermal Contact

Mechanical and thermal contact can be combined:

```
# Mechanical contact (contact module)
[Contact]
  [mortar]
    primary = bottom_top
    secondary = top_bottom
    formulation = mortar
    model = coulomb
    friction_coefficient = 0.3
  []
[]

# Thermal contact (heat_transfer module)
[ThermalContact]
  [gap_conductance]
    type = GapHeatTransfer
    variable = T
    primary = bottom_top
    secondary = top_bottom
    emissivity_primary = 0.8
    emissivity_secondary = 0.8
  []
[]
```

---

*This reference document was generated by reading the MOOSE source tree at
`modules/` in the `moose-next` repository. All class names, dependency
declarations, and example inputs are drawn directly from the source headers
and test files. For the most current list of input parameters for any class,
run `moose-opt --dump ClassName` or consult the online MOOSE documentation
at https://mooseframework.inl.gov.*
