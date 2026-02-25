# MOOSE Framework Architecture

**Multiphysics Object-Oriented Simulation Environment**

*Audience: C/C++/Java developers new to scientific computing frameworks.*

---

## Table of Contents

1. [The Thirty-Second View](#1-the-thirty-second-view)
2. [The Three-Layer Stack](#2-the-three-layer-stack)
3. [Directory Structure Map](#3-directory-structure-map)
4. [The MooseObject Class Hierarchy](#4-the-mooseobject-class-hierarchy)
5. [The Factory Pattern](#5-the-factory-pattern)
6. [The Action System: From Input to Objects](#6-the-action-system-from-input-to-objects)
7. [The Solve Loop](#7-the-solve-loop)
8. [The Material System](#8-the-material-system)
9. [The MultiApp and Transfer System](#9-the-multiapp-and-transfer-system)
10. [The Finite Volume Subsystem](#10-the-finite-volume-subsystem)
11. [The Mesh Generator Pipeline](#11-the-mesh-generator-pipeline)
12. [The Input Parameter System](#12-the-input-parameter-system)
13. [Parallel Architecture](#13-parallel-architecture)
14. [Communication Flow Diagrams](#14-communication-flow-diagrams)

---

## 1. The Thirty-Second View

MOOSE is a C++ finite-element multiphysics simulation framework that sits on top of two heavyweight scientific libraries: **libMesh** (finite-element mesh and DOF management) and **PETSc** (algebraic solvers). You write physics as small C++ classes called *Kernels*, *BoundaryConditions*, *Materials*, etc., describe your problem in a plain-text HIT input file, and MOOSE assembles and solves the resulting nonlinear system of equations. The framework handles MPI parallelism, mesh partitioning, ghost-element communication, thread-level parallelism for element loops, restart/checkpointing, and output to Exodus/VTK — so your physics code never touches MPI directly.

```
+------------------+
|   Input File     |  (HIT syntax: blocks, subblocks, key=value pairs)
|  problem.i       |
+--------+---------+
         | hit::parse()
         v
+------------------+         +------------------+
|   Moose::Builder | ------> | ActionWarehouse  |
|  (HIT Walker)    |         |  ordered task    |
|  Builder.h       |         |  graph           |
+------------------+         +--------+---------+
                                       | executeAllTasks()
                                       v
                              +------------------+
                              |   ActionFactory  |  resolves "type=" strings
                              |   (per-task)     |  to Action subclasses
                              +--------+---------+
                                       | Action::act()
                                       v
                              +------------------+
                              |     Factory      |  resolves "type=" strings
                              |   Factory.h      |  to MooseObject subclasses
                              +--------+---------+
                                       | Factory::create<T>()
                                       v
                         +---------------------------+
                         |       MooseObjects        |
                         |  Kernel, BC, Material,    |
                         |  UserObject, Postprocessor|
                         +-------------+-------------+
                                       |
                                       v
                         +---------------------------+
                         |      FEProblemBase        |
                         |  (owns all systems,       |
                         |   warehouses, Assembly)   |
                         +-------------+-------------+
                                       |
                          +------------+------------+
                          |                         |
                          v                         v
              +------------------+      +------------------+
              |   Executioner    |      |   Output System  |
              |  Steady/Transient|      | Exodus, CSV,     |
              |  Executioner.h   |      | Console, VTK...  |
              +--------+---------+      +------------------+
                       |
                       v
              +------------------+
              |  NonlinearSystem |  wraps libMesh::NonlinearImplicitSystem
              |  NonlinearSystem |
              |   .h / .C        |
              +--------+---------+
                       |
                       v
              +------------------+
              |   PETSc SNES     |  Newton / PJFNK nonlinear solver
              |   KSP + PC       |  linear solver + preconditioner
              +--------+---------+
                       |  calls residualFunction(), jacobianFunction()
                       v
              +------------------+
              |  Assembly +      |  element loops, quadrature, shape fns
              |  Kernel loops    |  ThreadedElementLoop
              +------------------+
```

Every component above corresponds to a real class in `framework/include/`. The following sections peel back each layer.

---

## 2. The Three-Layer Stack

MOOSE is best understood as three cooperating layers:

```
+=========================================================================+
|  LAYER 3 — MOOSE                                                        |
|                                                                         |
|  MooseApp           owns Factory, ActionWarehouse, MeshGeneratorSystem, |
|  FEProblemBase      owns NonlinearSystem, AuxiliarySystem, Assembly,    |
|                     MaterialWarehouse, Kernel/BC/UserObject warehouses  |
|  Executioner        Steady, Transient (+ time steppers, integrators)    |
|  Action system      Builder (HIT) -> ActionWarehouse -> Action::act()  |
|  Object system      Registry -> Factory -> MooseObject subclasses       |
|  Output system      OutputWarehouse -> Exodus, CSV, Console, ...        |
|  Parallel helpers   ParallelObject, Relationship Managers, ghosting     |
|                                                                         |
|  KEY CLASSES (framework/include/)                                       |
|    base/MooseApp.h          base/Factory.h          base/Registry.h    |
|    problems/FEProblemBase.h systems/NonlinearSystem.h                   |
|    executioners/Executioner.h, Transient.h, FEProblemSolve.h           |
|    actions/Action.h          actions/ActionWarehouse.h                  |
|    base/Assembly.h           utils/InputParameters.h                    |
+=========================================================================+
         |                            |                          |
         | libMesh::EquationSystems   | libMesh::MeshBase        | PETSc calls
         v                            v                          v
+========================+   +==================+   +========================+
|  LAYER 2 — PETSc       |   |  LAYER 1 — libMesh|   (libMesh wraps PETSc  |
|                        |   |                  |    solver interface)     |
|  SNES                  |   |  MeshBase /      |                          |
|    nonlinear solver    |   |  ReplicatedMesh / |   SNES: nonlinear solve |
|    Newton iterations   |   |  DistributedMesh  |   KSP:  linear solve    |
|                        |   |                  |   PC:   preconditioner  |
|  KSP                   |   |  Elem, Node,     |                          |
|    Krylov linear solver|   |  FEBase          |   Data flow:            |
|    GMRES, CG, ...      |   |                  |   MOOSE fills PETSc     |
|                        |   |  DofMap          |   vectors/matrices;     |
|  PC                    |   |  (global DOF     |   PETSc calls MOOSE     |
|    ILU, AMG, ...       |   |   numbering,     |   residual/Jacobian     |
|                        |   |   sparsity)      |   callbacks.            |
|  Mat, Vec              |   |                  |                          |
|    sparse matrices     |   |  QuadratureRule  |                          |
|    parallel vectors    |   |  (Gauss points)  |                          |
+========================+   +==================+   +========================+
```

### Data Flow Between Layers (Residual Evaluation)

```
FEProblemBase::computeResidualInternal()
        |
        v
NonlinearSystemBase::computeResidual()
        |
        | launches ThreadedElementLoop (OpenMP / TBB threads)
        v
ComputeResidualThread::onElement(const Elem * elem)
        |
        |-- Assembly::reinit(elem)           // libMesh: shape fns, JxW, q-pts
        |-- MaterialWarehouse::reinit(elem)  // fill MaterialProperty arrays
        |-- for each Kernel in warehouse:
        |       Kernel::computeResidual()
        |           -> for each qp: computeQpResidual()
        |               assembles into local dense vector _local_re
        |-- Assembly::addResidualLocal()     // libMesh: scatter to PETSc Vec
        v
libMesh::NumericVector<Number> (= PETSc Vec)
        |
        v
PETSc SNES receives residual, decides next Newton step
```

---

## 3. Directory Structure Map

The framework source is split into `framework/include/` (headers) and `framework/src/` (implementations). They have identical subdirectory names. Below is the full annotated tree, grouped by responsibility.

```
framework/
  include/                       (and mirror: src/)
  |
  +-- CORE INFRASTRUCTURE
  |   +-- base/                  Core classes: MooseApp, MooseBase, MooseObject,
  |   |                          Factory, Registry, Assembly, InputParameters (re-exported),
  |   |                          MeshGeneratorSystem, TheWarehouse, MooseError, MooseInit
  |   +-- parser/                HIT input file parsing: Builder (the HIT Walker),
  |   |                          Parser, Syntax, MooseSyntax, CommandLine
  |   +-- actions/               Action base class + ActionWarehouse + ActionFactory +
  |   |                          ~60 concrete Add*Action classes (one per object type)
  |   +-- utils/                 InputParameters, MooseUtils, MooseTypes, MooseEnum,
  |   |                          ExecFlagEnum, DataFileUtils, conversion helpers, math utils
  |   +-- interfaces/            Mix-in interfaces: Coupleable, BlockRestrictable,
  |   |                          BoundaryRestrictable, SetupInterface, PerfGraphInterface,
  |   |                          PostprocessorInterface, UserObjectInterface, Restartable, ...
  |   +-- warehouses/            MooseObjectWarehouse, MooseObjectTagWarehouse,
  |   |                          ExecuteMooseObjectWarehouse (exec-flag scheduling)
  |
  +-- PROBLEMS & SYSTEMS
  |   +-- problems/              FEProblemBase (central hub), FEProblem (concrete),
  |   |                          SubProblem, DisplacedProblem, EigenProblem, ExternalProblem
  |   +-- systems/               SystemBase, SolverSystem, NonlinearSystemBase,
  |   |                          NonlinearSystem, LinearSystem, AuxiliarySystem,
  |   |                          MooseEigenSystem, DisplacedSystem
  |
  +-- PHYSICS OBJECTS
  |   +-- kernels/               FE volume integrals (strong-form residuals).
  |   |                          KernelBase -> Kernel, ADKernel, ArrayKernel, VectorKernel.
  |   |                          Specific: Diffusion, Reaction, TimeDerivative, BodyForce, ...
  |   +-- bcs/                   Boundary conditions: IntegratedBC (Neumann-type),
  |   |                          NodalBC (Dirichlet), DirichletBC, NeumannBC, RobinBC, ...
  |   +-- dgkernels/             Discontinuous Galerkin interior-face terms
  |   +-- interfacekernels/      Terms on internal interfaces between subdomains
  |   +-- nodalkernels/          Residuals evaluated at nodes (ODEs, lumped)
  |   +-- scalarkernels/         Scalar (0-D) equation kernels
  |   +-- hdgkernels/            Hybridized DG kernels
  |
  +-- FINITE VOLUME
  |   +-- fvkernels/             FVKernel base, FVFluxKernel (face loop, divergence theorem),
  |   |                          FVElementalKernel (cell loop), FVTimeKernel, FVDiffusion, ...
  |   +-- fvbcs/                 FV boundary conditions: FVDirichletBC, FVNeumannBC, ...
  |   +-- fvics/                 FV initial conditions
  |   +-- fviks/                 FV interface kernels
  |   +-- linearfvkernels/       Kernels for the linear FV solver pathway
  |   +-- linearfvbcs/           BCs for the linear FV solver pathway
  |   +-- limiters/              Slope limiters for FV advection (Minmod, VanLeer, ...)
  |
  +-- MATERIALS
  |   +-- materials/             MaterialBase -> Material. declareProperty<T>(),
  |   |                          getMaterialProperty<T>(). Stateful (old/older) support.
  |   |                          DerivativeMaterial*, GenericConstant*, InterfaceMaterial, ...
  |   +-- functormaterials/      FunctorMaterial: lazy evaluation via addFunctorProperty().
  |   |                          FunctorMaterialProperty wraps a lambda called on demand.
  |
  +-- AUXILIARY VARIABLES & ICS
  |   +-- auxkernels/            AuxKernels fill auxiliary variable fields (post-process).
  |   +-- auxscalarkernels/      Scalar auxiliary kernels
  |   +-- ics/                   InitialConditions: set solution vector at t=0
  |
  +-- USER OBJECTS & POSTPROCESSORS
  |   +-- userobjects/           UserObject base hierarchy: GeneralUserObject,
  |   |                          ElementUserObject, SideUserObject, NodalUserObject,
  |   |                          InterfaceUserObject. Execute on EXEC flags.
  |   +-- postprocessors/        Scalar reduction objects (area, volume integrals, norms).
  |   |                          Inherits from GeneralUserObject + Reporter.
  |   +-- vectorpostprocessors/  1-D vectors of values per timestep.
  |   +-- reporters/             Typed key-value reporter data (modern replacement for PP/VPP).
  |   +-- samplers/              Monte Carlo sample generation.
  |
  +-- EXECUTIONERS & SOLVERS
  |   +-- executioners/          Executioner (base), Steady, Transient (-> TransientBase),
  |   |                          FEProblemSolve, FixedPointSolve, PicardSolve, SecantSolve,
  |   |                          SolveObject, NonlinearSolveObject, MultiSystemSolveObject
  |   +-- executors/             Lower-level executor objects for modular solve pipelines
  |   +-- timeintegrators/       ImplicitEuler, CrankNicolson, BDF2, NewmarkBeta, Explicit* ...
  |   +-- timesteppers/          ConstantDT, IterationAdaptiveDT, FunctionDT, PostprocessorDT ...
  |   +-- predictors/            Solution predictors for time stepping
  |   +-- correctors/            Post-step solution correction
  |   +-- linesearches/          Custom PETSc line searches
  |   +-- convergence/           Convergence assessment objects
  |   +-- preconditioners/       MoosePreconditioner wrappers (FieldSplit, ASM, ...)
  |   +-- splits/                FieldSplit preconditioner block definitions
  |
  +-- MESH
  |   +-- mesh/                  MooseMesh (wraps libMesh::MeshBase), GeneratedMesh,
  |   |                          FileMesh, AnnularMesh, FaceInfo (FV face metadata)
  |   +-- meshgenerators/        ~50+ MeshGenerator subclasses (see Section 11)
  |   +-- meshmodifiers/         Post-generation mesh modifications
  |   +-- meshdivisions/         Mesh subdivision strategies for samplers
  |   +-- geomsearch/            Geometric search: PenetrationLocator, NearestNodeLocator
  |   +-- partitioner/           Custom mesh partitioners (Petsc, Parmetis, Hierarch.)
  |
  +-- CONSTRAINTS & CONTACT
  |   +-- constraints/           NodeFaceConstraint, EqualDOF, TiedValueConstraint
  |
  +-- MULTIAPP & TRANSFERS
  |   +-- multiapps/             MultiApp base, TransientMultiApp, FullSolveMultiApp,
  |   |                          CentroidMultiApp, QuadraturePointMultiApp
  |   +-- transfers/             MultiAppTransfer base + ~15 concrete Transfer types
  |   |                          (NearestNode, Projection, UserObject, Reporter, ...)
  |
  +-- CONTROLS & FUNCTIONS
  |   +-- controls/              Control objects: modify parameters at runtime
  |   +-- chaincontrols/         Chain of controls for complex runtime logic
  |   +-- functions/             Function (analytic expressions), ParsedFunction, ...
  |   +-- positions/             Named position sets for MultiApp placement
  |   +-- times/                 Named time sets (e.g., for time-sequence steppers)
  |
  +-- INDICATORS & MARKERS (AMR)
  |   +-- indicators/            Estimate local solution error (for AMR)
  |   +-- markers/               Flag elements for refinement/coarsening
  |   +-- dampers/               Limit Newton update magnitude per DOF
  |   +-- bounds/                Enforce variable bounds
  |
  +-- OUTPUT
  |   +-- outputs/               Output base + Exodus, CSV, Console, Checkpoint,
  |   |                          Nemesis, VTK/Tecplot, JSON, GMV, Gnuplot, ...
  |   +-- restart/               Checkpoint/restart data serialization
  |
  +-- COMPUTE LOOPS (internal)
  |   +-- loops/                 ThreadedElementLoop infrastructure.
  |   |                          ComputeResidualThread, ComputeJacobianThread,
  |   |                          ComputeMaterialsObjectThread, ComputeUserObjectsThread,
  |   |                          ComputeFVFluxThread, ComputeNodalKernelsThread, ...
  |
  +-- DISTRIBUTIONS & SAMPLERS
  |   +-- distributions/         Probability distributions (Normal, Uniform, Weibull, ...)
  |
  +-- RELATIONSHIP MANAGERS
  |   +-- relationshipmanagers/  Tell libMesh which ghost elements are needed for
  |                              algebraic coupling and geometric stencil requirements.
  |
  +-- VARIABLES
  |   +-- variables/             MooseVariable (FE scalar), MooseVariableFE<T>,
  |                              MooseVariableFV<T>, MooseVariableScalar,
  |                              ArrayMooseVariable, MooseVariableData<T>
  |
  +-- OPTIONAL / EXPERIMENTAL
      +-- physics/               PhysicsBase: higher-level physics convenience classes
      +-- hdgkernels/            Hybridized DG (experimental)
      +-- libtorch/              LibTorch ML integration (if enabled)
      +-- kokkos/                Kokkos GPU/CPU portability layer (if enabled)
      +-- neml2/                 NEML2 material model integration (if enabled)
      +-- mfem/                  MFEM backend (experimental)
      +-- csg/                   CSG mesh generation (OpenCASCADE-backed)
```

---

## 4. The MooseObject Class Hierarchy

Every user-facing class in MOOSE ultimately inherits from `MooseBase`, and most from `MooseObject`. The hierarchy is wide rather than deep, using multiple inheritance for mix-in interfaces. Here is the primary structural tree with key virtual methods shown.

```
ConsoleStreamInterface
    |
    +-- MooseBase  (framework/include/base/MooseBase.h)
    |       name(), type(), typeAndName(), parameters()
    |       _app, _type, _name, _pars
    |       |
    |       +-- Action  (framework/include/actions/Action.h)
    |       |       virtual act() = 0
    |       |       (also: ParallelParamObject, MeshMetaDataInterface,
    |       |        PerfGraphInterface, SolutionInvalidInterface)
    |       |
    |       +-- MooseObject  (framework/include/base/MooseObject.h)
    |               virtual enabled() const
    |               (also: ParallelParamObject, SolutionInvalidInterface,
    |                enable_shared_from_this<MooseObject>)
    |               |
    |               +-- PHYSICS RESIDUAL OBJECTS
    |               |   |
    |               |   +-- ResidualObject  (base/ResidualObject.h)
    |               |   |       computeResidual(), computeJacobian()
    |               |   |       computeResidualAndJacobian()
    |               |   |
    |               |   +-- KernelBase  (kernels/ -> ResidualObject)
    |               |   |   +-- Kernel
    |               |   |   |     computeQpResidual() [pure virtual]
    |               |   |   |     computeQpJacobian()
    |               |   |   |     _u, _grad_u, _phi, _test, _qp
    |               |   |   |
    |               |   |   +-- ADKernelTempl<T>   (AD auto-diff path)
    |               |   |   |     computeQpResidual() -> ADReal
    |               |   |   |
    |               |   |   +-- ArrayKernel         (vector of DOFs per node)
    |               |   |   +-- VectorKernel         (Nedelec / H(curl) variables)
    |               |   |
    |               |   +-- BoundaryCondition  (bcs/)
    |               |   |   +-- IntegratedBCBase -> IntegratedBC, ADIntegratedBC
    |               |   |   |     computeQpResidual(), computeQpJacobian()
    |               |   |   +-- NodalBCBase -> NodalBC, ADNodalBC
    |               |   |         computeQpResidual()
    |               |   |
    |               |   +-- DGKernelBase  (dgkernels/)
    |               |   |     computeQpResidual(Moose::DGResidualType)
    |               |   |
    |               |   +-- InterfaceKernelBase  (interfacekernels/)
    |               |   |     computeQpResidual(Moose::InterfaceKernelType)
    |               |   |
    |               |   +-- FVKernel  (fvkernels/)
    |               |       +-- FVFluxKernel
    |               |       |     computeQpResidual() -> ADReal  [face-centered]
    |               |       +-- FVElementalKernel
    |               |             computeQpResidual() -> ADReal  [cell-centered]
    |               |
    |               +-- MaterialBase  (materials/MaterialBase.h)
    |               |   computeProperties() [pure virtual]
    |               |   initStatefulProperties()
    |               |   declareProperty<T>(), declareADProperty<T>()
    |               |   _qp (current quadrature point index)
    |               |   |
    |               |   +-- Material  (materials/Material.h)
    |               |   |   computeQpProperties() [pure virtual]
    |               |   |   getMaterialProperty<T>()  -- consumer side
    |               |   |   getMaterialPropertyOld<T>(), getMaterialPropertyOlder<T>()
    |               |   |   |
    |               |   |   +-- ADMaterial  (AD-only override)
    |               |   |   +-- GenericConstantMaterial
    |               |   |   +-- DerivativeParsedMaterial
    |               |   |   +-- (user-defined materials)
    |               |   |
    |               |   +-- FunctorMaterial  (functormaterials/)
    |               |         computeProperties() {}  [empty, lazy]
    |               |         addFunctorProperty<T>(name, lambda)
    |               |         (evaluated only when consumer calls the functor)
    |               |
    |               +-- UserObjectBase / UserObject
    |               |   |
    |               |   +-- GeneralUserObject
    |               |   |     initialize(), execute(), finalize()  [pure virtual]
    |               |   |     |
    |               |   |     +-- Postprocessor  (also Reporter)
    |               |   |           getValue() -> Real
    |               |   |
    |               |   +-- ElementUserObject
    |               |   |     executeOnElement()  [per element]
    |               |   |
    |               |   +-- SideUserObject
    |               |   |     executeOnBoundary()  [per boundary face]
    |               |   |
    |               |   +-- NodalUserObject
    |               |         executeOnNode()  [per node]
    |               |
    |               +-- AuxKernel  (auxkernels/)
    |               |     computeValue() -> Real  [pure virtual]
    |               |     fills AuxiliarySystem DOFs (not in nonlinear solve)
    |               |
    |               +-- Executioner  (executioners/Executioner.h)
    |               |   virtual init(), execute(), preExecute(), postExecute()
    |               |   |
    |               |   +-- Steady
    |               |   |     execute() -> FEProblemSolve::solve()
    |               |   |
    |               |   +-- TransientBase
    |               |       execute() -> loop: takeStep() while keepGoing()
    |               |       |
    |               |       +-- Transient  (owns FEProblemSolve _feproblem_solve)
    |               |       +-- EigenExecutionerBase -> Eigenvalue, InversePowerMethod
    |               |
    |               +-- MeshGenerator  (meshgenerators/)
    |               |     generate() -> std::unique_ptr<MeshBase>  [pure virtual]
    |               |     (see Section 11)
    |               |
    |               +-- Output  (outputs/Output.h)
    |               |     output() [pure virtual]
    |               |     |
    |               |     +-- FileOutput -> Exodus, Nemesis, CSV, VTK, Checkpoint
    |               |     +-- PetscOutput -> Console (also prints PETSc log)
    |               |     +-- Console
    |               |
    |               +-- MultiApp  (multiapps/MultiApp.h)
    |               |     initialSetup(), solveStep()
    |               |     createApp(), resetApp(), restoreApp()
    |               |
    |               +-- Transfer  (transfers/Transfer.h)
    |               |     execute() [pure virtual]
    |               |     |
    |               |     +-- MultiAppTransfer
    |               |           +-- MultiAppCopyTransfer
    |               |           +-- MultiAppNearestNodeTransfer
    |               |           +-- MultiAppProjectionTransfer
    |               |           +-- MultiAppReporterTransfer
    |               |           +-- (many more...)
    |               |
    |               +-- Function  (functions/)
    |               |     value(Real t, const Point & p) -> Real
    |               |
    |               +-- Distribution  (distributions/)
    |               +-- Sampler       (samplers/)
    |               +-- TimeStepper   (timesteppers/)
    |               +-- TimeIntegrator(timeintegrators/)
    |               +-- Damper        (dampers/)
    |               +-- Indicator     (indicators/)
    |               +-- Marker        (markers/)
    |               +-- Constraint    (constraints/)
    +-- (many mix-in interfaces, not shown in full — see interfaces/ directory)
```

---

## 5. The Factory Pattern

### Problem Solved

An input file says `type = Diffusion` inside a `[Kernels]` block. At compile time, MOOSE does not know which class to instantiate — that is decided by user code. The Factory pattern maps string names to constructor calls, allowing new physics to be added without modifying framework source.

### Registration: `registerMooseObject`

In the `.C` file for each physics class, a single macro call registers it before `main()` runs:

```cpp
// framework/src/kernels/Diffusion.C  (example)
registerMooseObject("MooseApp", Diffusion);
```

This macro expands (see `framework/include/base/Registry.h`) to:

```cpp
static char dummyvar_for_registering_obj_Diffusion_<counter> =
    Registry::add<Diffusion>({"MooseApp", "Diffusion", "", "", __FILE__, __LINE__, "", ""});
```

The `Registry::add<T>()` template method:
1. Creates a `RegistryEntry<T>` object (a concrete subclass of `RegistryEntryBase`).
2. Stores it in `Registry::_per_label_objects["MooseApp"]`.
3. Records the mapping `typeid(T).name() -> "Diffusion"` for reverse lookup.

The `RegistryEntry<T>` struct holds three virtual factory methods:

```cpp
// Simplified from Registry.h
template <typename T>
struct RegistryEntry : public RegistryEntryBase {
    std::unique_ptr<MooseObject> build(const InputParameters & params) override {
        return std::make_unique<T>(params);   // T::T(const InputParameters &)
    }
    InputParameters buildParameters() override {
        return T::validParams();              // calls static method on T
    }
    // ...
};
```

### Loading Into the Factory

When `MooseApp` is constructed, it calls:

```cpp
Registry::registerObjectsTo(factory, {"MooseApp"});
```

This iterates `_per_label_objects["MooseApp"]` and calls `Factory::reg()` for each entry, populating `Factory::_name_to_object` (a `std::map<string, shared_ptr<RegistryEntryBase>>`).

### Object Creation

When an Action processes `type = Diffusion`, it calls:

```cpp
// Typical Action::act() pattern
auto params = _factory.getValidParams("Diffusion");
// ... fill params from input file ...
auto kernel = _factory.create<KernelBase>("Diffusion", "diff", params, tid);
```

`Factory::create<T>()` (see `framework/include/base/Factory.h`):

```
Factory::create<T>(obj_name, name, params, tid)
    |
    +-- look up _name_to_object[obj_name]  -> RegistryEntry<Diffusion>
    +-- Factory::initialize(type, name, params, tid)
    |       creates InputParameters in InputParameterWarehouse
    +-- _currently_constructing.push_back(&params)
    +-- entry->buildShared(params)
    |       -> std::make_shared<Diffusion>(params)
    +-- _currently_constructing.pop_back()
    +-- Factory::finalize(type, object)
    |       sanity-checks params
    +-- createCheckObjectType<T>(object, obj_name)
    +-- return std::static_pointer_cast<T>(object)
```

### Aliases and Deprecation

The registry supports object aliasing and deprecation:

```cpp
registerMooseObjectAliased("MooseApp", MyClass, "OldName")
registerMooseObjectDeprecated("MooseApp", OldClass, "06/01/2026 00:00")
registerMooseObjectReplaced("MooseApp", OldClass, "06/01/2026 00:00", NewClass)
```

### Actions Use the Same Pattern

`registerMooseAction("MooseApp", AddKernelAction, "add_kernel")` registers an Action class associated with the task name `"add_kernel"`. The `ActionFactory` mirrors `Factory` for Actions.

### Registration Flow Diagram

```
Application static initialization
        |
        | (before main())
        v
Registry::add<Diffusion>(...)    <- registerMooseObject("MooseApp", Diffusion)
Registry::add<ADDiffusion>(...)  <- registerMooseObject("MooseApp", ADDiffusion)
Registry::add<Reaction>(...)     ...
...
        |
        v
MooseApp constructor
        |
        +-- Registry::registerObjectsTo(factory, labels)
        |       iterates _per_label_objects, calls Factory::reg()
        |
        +-- Registry::registerActionsTo(action_factory, labels)
                iterates _per_label_actions, calls ActionFactory::reg()
        |
        v
Factory::_name_to_object populated:
    "Diffusion"       -> RegistryEntry<Diffusion>
    "ADDiffusion"     -> RegistryEntry<ADDiffusion>
    "Reaction"        -> RegistryEntry<Reaction>
    ...
        |
        | (during input file processing)
        v
AddKernelAction::act()
    -> Factory::create<KernelBase>("Diffusion", "diff", params)
    -> RegistryEntry<Diffusion>::buildShared(params)
    -> new Diffusion(params)
```

---

## 6. The Action System: From Input to Objects

### What Is an Action?

An Action is a task object that *sets up* the simulation — it does not participate in the solve. Think of it like a build script: it reads parameters from the input file and calls framework APIs to create MooseObjects, register systems, set up meshes, etc.

Every concrete Add*Action calls `_problem->addKernel(...)`, `_problem->addBC(...)`, etc. which stores the object in the appropriate `MooseObjectWarehouse` inside `FEProblemBase`.

### HIT Input Blocks Map to Actions

```
Input File                     Action triggered
-------------------------------+----------------------------------
[Mesh]                         SetupMeshAction, AddMeshGeneratorAction
[Variables]                    AddVariableAction
[AuxVariables]                 AddAuxVariableAction
[Kernels]                      AddKernelAction       (task: "add_kernel")
  [./diff]
    type = Diffusion
  [../]
[BCs]                          AddBCAction           (task: "add_bc")
[Materials]                    AddMaterialAction     (task: "add_material")
[Executioner]                  CreateExecutionerAction
[Outputs]                      AddOutputAction
[MultiApps]                    AddMultiAppAction
[Transfers]                    AddTransferAction
```

### Task Ordering in ActionWarehouse

The ActionWarehouse maintains an ordered list of task names. Tasks are declared in `Moose::registerAll()` via calls like:

```
registerTask("setup_mesh",          /*is_required=*/false);
registerTask("add_mesh_generator",  false);
registerTask("create_mesh",         true);
registerTask("add_variable",        false);
registerTask("add_aux_variable",    false);
registerTask("add_kernel",          false);
registerTask("add_bc",              false);
registerTask("add_material",        false);
registerTask("init_problem",        true);
registerTask("setup_executioner",   false);
...
```

`ActionWarehouse::executeAllTasks()` iterates tasks in dependency-resolved order, calling `action->timedAct()` on every Action registered for each task.

### Example: `[Kernels]/[diff]/type=Diffusion`

```
  Input File Parser
  (Moose::Builder, builder.h)
        |
        | walk() called for each HIT node
        | sees block path "[Kernels]/[diff]", param "type = Diffusion"
        v
  ActionFactory::create("AddKernelAction", ...)
        |
        v
  ActionWarehouse::addActionBlock(shared_ptr<AddKernelAction>)
        |
        | ... later, during executeAllTasks("add_kernel") ...
        v
  AddKernelAction::act()                        (actions/AddKernelAction.C)
        |
        +-- type_str = getParam<std::string>("type")  // "Diffusion"
        +-- params = _factory.getValidParams(type_str)
        +-- _moose_object_pars.applyInputFileParameters(params, input_node)
        +-- _problem->addKernel(type_str, name, params, tid)
                |
                v
          FEProblemBase::addKernel(type, name, params, tid)
                |
                +-- kernel = _factory.create<KernelBase>(type, name, params, tid)
                |       // Factory resolves "Diffusion" -> new Diffusion(params)
                +-- _nl->addKernel(kernel, tid)
                        // stores in NonlinearSystemBase::_kernels warehouse

+----------------------------------------------------------+
|  Diagram: Input to Object                                |
|                                                          |
|  problem.i                                               |
|    [Kernels]                                             |
|      [./diff]                                            |
|        type = Diffusion  ----+                           |
|      [../]                   |                           |
|                              v                           |
|                       Builder::walk()                    |
|                        (HIT walker)                      |
|                              |                           |
|                              v                           |
|                    ActionFactory::create                 |
|                    ("AddKernelAction")                   |
|                              |                           |
|                              v                           |
|                    ActionWarehouse stores                 |
|                    AddKernelAction instance              |
|                              |                           |
|                (task "add_kernel" fires)                 |
|                              |                           |
|                              v                           |
|                   AddKernelAction::act()                 |
|                              |                           |
|                    Factory::create<KernelBase>           |
|                    ("Diffusion", "diff", params)         |
|                              |                           |
|                              v                           |
|                   RegistryEntry<Diffusion>               |
|                   ::buildShared(params)                  |
|                              |                           |
|                              v                           |
|                   new Diffusion(params)                  |
|                              |                           |
|                              v                           |
|                   FEProblemBase stores in                |
|                   NonlinearSystemBase::_kernels          |
+----------------------------------------------------------+
```

### Full Task Sequence (Abbreviated)

```
setup_mesh
  -> SetupMeshAction: constructs MooseMesh
add_mesh_generator
  -> AddMeshGeneratorAction: stores generator params
create_mesh
  -> MeshGeneratorSystem::executeMeshGenerators(): runs DAG, builds mesh
setup_aux_variables
add_variable
  -> AddVariableAction: calls _problem->addVariable()
add_aux_variable
add_kernel
  -> AddKernelAction: (see above)
add_bc
add_material
add_user_object
add_postprocessor
add_transfer
add_multi_app
add_output
init_problem
  -> FEProblemBase::init(): sets up DOF map, allocates PETSc vectors
setup_executioner
  -> CreateExecutionerAction: instantiates Steady or Transient
execute
  -> ExecuteMeshGenerators, then Executioner::execute()
```

---

## 7. The Solve Loop

This section traces the exact call chain from `Executioner::execute()` all the way down to `computeQpResidual()` for a `Transient` simulation.

### Top-Level Loop (Transient)

```
MooseApp::run()
    |
    v
Executioner::execute()          (framework/include/executioners/Executioner.h)
    |  [overridden by TransientBase]
    v
TransientBase::execute()        (executioners/TransientBase.h)
    |
    +-- preExecute()            (output initial state, etc.)
    +-- while (keepGoing()):
    |       |
    |       +-- computeDT()     (call TimeStepper::computeDT())
    |       +-- takeStep()
    |       |       |
    |       |       v
    |       |   TransientBase::takeStep()
    |       |       |
    |       |       +-- _fixed_point_solve->solve()    [FixedPointSolve / PicardSolve]
    |       |       |       |
    |       |       |       v
    |       |       |   FEProblemSolve::solve()        (executioners/FEProblemSolve.h)
    |       |       |       |
    |       |       |       +-- [optional grid sequencing outer loop]
    |       |       |       +-- _problem.solve()
    |       |       |               |
    |       |       |               v
    |       |       |           FEProblemBase::solve()  (problems/FEProblemBase.h)
    |       |       |               |
    |       |       |               +-- for each NonlinearSystem _nl:
    |       |       |                       _nl->solve()
    |       |       |                           |
    |       |       |                           v
    |       |       |                  NonlinearSystem::solve()
    |       |       |                           |
    |       |       |                           v
    |       |       |                  libMesh::NonlinearImplicitSystem::solve()
    |       |       |                           |
    |       |       |                           v
    |       |       |                  PETSc SNES solve
    |       |       |                  (SNESSolve or SNESSolve_Private)
    |       |       |
    |       +-- endStep()       (accept/reject, advance time)
    |       +-- incrementStepOrReject()
    |
    +-- postExecute()           (finalize outputs)
```

### PETSc Callbacks Into MOOSE

PETSc SNES calls back into MOOSE for residual and Jacobian evaluations. The callback is wired in `NonlinearSystem`:

```cpp
// Simplified from NonlinearSystem.C
_nl_implicit_sys.nonlinear_solver->residual =
    [this](const NumericVector<Number> & X, NumericVector<Number> & R, NonlinearImplicitSystem &) {
        _nl_residual_functor(X, R, sys);  // ComputeResidualFunctor
    };
_nl_implicit_sys.nonlinear_solver->jacobian =
    [this](...) { ... };  // ComputeJacobianFunctor
```

### Residual Evaluation Path

```
PETSc SNES calls residualFunction(SNES, Vec X, Vec F, void* ctx)
    |
    v
ComputeResidualFunctor::operator()(X, R, sys)
    |
    v
FEProblemBase::computeResidualInternal(X, R, tags)
    |
    +-- updateSolution(X)          // scatter X to local ghosted vector
    +-- computeResidualTags(tags)
            |
            v
    NonlinearSystemBase::computeResidual(R, tags)
            |
            +-- _fe_problem.reinitMaterials()    // reset material storage
            +-- Parallel_reduce over mesh using ComputeResidualThread
                    |
                    v
            [For each local Elem e on thread T]:
            ComputeResidualThread::onElement(e)
                    |
                    +-- _fe_problem.prepare(e, tid)
                    |       Assembly::reinit(e)
                    |           // libMesh: compute shape fns phi, dphi, JxW
                    |           // at all quadrature points
                    |
                    +-- _fe_problem.reinitMaterialsOnElement(e, tid)
                    |       // MaterialWarehouse: call computeProperties()
                    |       // on all materials active on this subdomain
                    |
                    +-- for each Kernel k active on e:
                    |       k->computeResidual()
                    |           |
                    |           +-- for qp in 0..nQP:
                    |                   _re(i) += _JxW[qp] * _coord[qp]
                    |                             * _test[i][qp]
                    |                             * k->computeQpResidual()
                    |                                 // pure virtual, user overrides
                    |
                    +-- Assembly::addResidualLocal(tid)
                            // scatter _re to global PETSc Vec F

            [For each boundary face on boundary BC]:
            ComputeResidualThread::onBoundary(e, side, bnd_id)
                    // similar: BC::computeResidual() -> computeQpResidual()
```

### Jacobian Evaluation Path

```
PETSc SNES calls jacobianFunction(SNES, Vec X, Mat A, Mat B, void* ctx)
    |
    v
FEProblemBase::computeJacobian(X, A, tags)
    |
    v
ComputeJacobianThread::onElement(e)
    |
    +-- Assembly::reinit(e)            // same as residual
    +-- reinitMaterialsOnElement(e)    // same as residual
    +-- for each Kernel k:
    |       k->computeJacobian()
    |           // for i in test DOFs:
    |           //   for j in trial DOFs:
    |           //     _ke(i,j) += _JxW[qp]*_test[i]*computeQpJacobian()*_phi[j]
    +-- Assembly::addJacobianLocal(tid)
            // scatter _ke to global PETSc Mat A
```

### AD (Automatic Differentiation) Path

For `ADKernel` subclasses, MOOSE uses forward-mode AD via `DualNumber<Real, DualNumberSupportedVector>`. The Jacobian is derived automatically:

```
ADKernelTempl::computeResidual()
    |
    +-- seed DOF values with dual numbers: _u[qp] = DualNumber(val, seed)
    +-- computeQpResidual() -> ADReal  (carries derivatives w.r.t. all DOFs)
    +-- extract .value() for residual
    +-- extract .derivatives() for Jacobian entries
    |   (no separate computeQpJacobian() needed)
```

### Element Loop Threading

`ThreadedElementLoop` uses `libMesh::Threads::parallel_reduce()`:

```
[Thread 0]  elements 0..N/4     ComputeResidualThread (copy per thread)
[Thread 1]  elements N/4..N/2   ComputeResidualThread (copy per thread)
[Thread 2]  elements N/2..3N/4  ComputeResidualThread (copy per thread)
[Thread 3]  elements 3N/4..N    ComputeResidualThread (copy per thread)
    |
    v (join phase)
All local dense vectors assembled into global PETSc vectors
(thread-safe: each thread accumulates into thread-local storage,
 final scatter is serialized by libMesh)
```

---

## 8. The Material System

### Overview

Materials compute physical properties (conductivity, viscosity, stress, etc.) at every quadrature point of every element before the Kernel loops need those values. The system follows a **producer/consumer** pattern:

- **Producer**: A `Material` subclass calls `declareProperty<T>("thermal_conductivity")` — this registers the property and returns a writable reference.
- **Consumer**: A `Kernel` or other object calls `getMaterialProperty<T>("thermal_conductivity")` — this registers a read dependency and returns a const reference.

MOOSE resolves the dependency graph automatically and ensures materials are computed in the right order.

### Class Hierarchy (Material side)

```
MooseObject
    |
    +-- MaterialBase  (materials/MaterialBase.h)
    |   Mixin interfaces: BlockRestrictable, BoundaryRestrictable,
    |   SetupInterface, Coupleable, FunctionInterface, Restartable, ...
    |
    |   computeProperties()         [pure virtual]
    |   initStatefulProperties()    [virtual]
    |   declareProperty<T>(name)    -> MaterialProperty<T>&   (writable)
    |   declareADProperty<T>(name)  -> ADMaterialProperty<T>& (AD-aware)
    |   _qp                         current quadrature point index
    |
    +-- Material  (materials/Material.h)
    |   Adds: Coupleable, MaterialPropertyInterface
    |
    |   computeQpProperties()       [pure virtual - user overrides THIS]
    |   computeProperties()         calls computeQpProperties() for each _qp
    |   getMaterialProperty<T>(name, state=0)   -> const MaterialProperty<T>&
    |   getMaterialPropertyOld<T>(name)          state=1
    |   getMaterialPropertyOlder<T>(name)        state=2
    |
    +-- FunctorMaterial  (functormaterials/FunctorMaterial.h)
        computeProperties() {}           [intentionally empty]
        addFunctorProperty<T>(name, lambda)
        (lazy: property evaluated only when consumer calls the functor)
```

### Stateful Properties

The material system maintains up to three time-level "copies" of each property:

```
Material Property Storage (MaterialPropertyStorage)

              timestep n-2         timestep n-1         timestep n (current)
              (older)              (old)                (current)
+-----------+-----------+--------+-----------+--------+-----------+
| prop name | prop[qp]  |        | prop[qp]  |        | prop[qp]  |
+-----------+-----------+--------+-----------+--------+-----------+
| "damage"  |  0.12     |        |  0.14     |        |  0.0 (?)  |
| "stress"  | [tensor]  |        | [tensor]  |        |  [...]    |
+-----------+-----------+--------+-----------+--------+-----------+
        ^                               ^                    ^
getMaterialPropertyOlder()     getMaterialPropertyOld()  getMaterialProperty()
```

Registration of stateful access is done at construction time, causing `MaterialPropertyStorage` to allocate the extra time levels. State is advanced by `FEProblemBase::advanceState()` at the end of each successful timestep.

### FunctorMaterial: Lazy Evaluation

Traditional `Material::computeQpProperties()` is **eager**: it computes all declared properties for every element, every quadrature point, every residual evaluation — even if most properties are never consumed on that element.

`FunctorMaterial` is **lazy**: the property value is a `FunctorMaterialProperty` holding a C++ lambda. The lambda runs only when a consumer actually calls the functor:

```
// Producer (FunctorMaterial subclass):
addFunctorProperty<Real>("mu",
    [this](const auto & r, const auto & t) -> ADReal {
        return _mu_0 * std::exp(-_activation_energy / _gas_constant / r.getElem().temperature);
    });

// Consumer (FVFluxKernel or similar):
const auto & mu = getFunctor<ADReal>("mu");
// ...
ADReal val = mu(face_arg, time_arg);  // lambda runs here, on demand
```

The functor argument types (`ElemArg`, `FaceArg`, `ElemSideQpArg`, etc.) are defined in `MooseFunctorArguments.h` and carry enough context (element pointer, side index, quadrature point) for the lambda to compute the right value.

### Material Compute Sequence Per Element

```
(Inside ComputeMaterialsObjectThread::onElement)
1. MaterialWarehouse::residualSetup()    [first residual call each timestep]
2. for each Material m active on subdomain:
       m->reinit(elem, qrule)            [update _qp_point, _JxW, etc.]
3. Dependency-ordered loop:
       for each Material m (ordered by DependencyResolver):
           m->computeProperties()
               for _qp = 0..nQP:
                   m->computeQpProperties()  // fills _prop_name[_qp]
4. Kernel loops can now read MaterialProperty values via const refs
```

---

## 9. The MultiApp and Transfer System

### Concept

MultiApp enables multiscale or multiphysics coupling by running sub-applications (sub-apps) as child `MooseApp` instances owned by the parent app. Each sub-app has its own `FEProblemBase`, mesh, variables, and executioner. Data moves between parent and sub-apps through `Transfer` objects.

### Parallel Decomposition

```
MPI World (16 ranks total)
|
+-- Parent App (ranks 0..15)
|   FEProblemBase, NonlinearSystem, ...
|   MPI_Comm = MPI_COMM_WORLD
|
+-- MultiApp "micro" (4 sub-app instances)
    |
    +-- sub-app[0]  on ranks  0..3   (MPI_Comm = sub_comm_0)
    +-- sub-app[1]  on ranks  4..7   (MPI_Comm = sub_comm_1)
    +-- sub-app[2]  on ranks  8..11  (MPI_Comm = sub_comm_2)
    +-- sub-app[3]  on ranks 12..15  (MPI_Comm = sub_comm_3)
```

Each sub-app gets a disjoint communicator carved out by `MultiApp::init()`. If there are more sub-apps than ranks, the remaining sub-apps are run serially on the same ranks (batch mode).

### MultiApp Class Hierarchy

```
MooseObject
    |
    +-- MultiApp  (multiapps/MultiApp.h)
    |   createApp(i, start_time)   -- constructs i-th sub MooseApp
    |   solveStep(dt, target_time) -- advance all sub-apps one step
    |   resetApp(i, time)          -- restore sub-app to checkpoint state
    |   getLocalApp(i)             -- returns local MooseApp pointer
    |
    +-- TransientMultiApp          -- sub-apps use Transient executioner
    +-- FullSolveMultiApp          -- sub-apps use Steady executioner per call
    +-- CentroidMultiApp           -- one sub-app per element centroid
    +-- QuadraturePointMultiApp    -- one sub-app per quadrature point
```

### Transfer Class Hierarchy

```
MooseObject
    |
    +-- Transfer
    |   execute() [pure virtual]
    |   direction: FROM_MULTIAPP or TO_MULTIAPP (or BETWEEN)
    |
    +-- MultiAppTransfer
        _from_multi_app, _to_multi_app (shared_ptr<MultiApp>)
        |
        +-- MultiAppCopyTransfer         exact DOF copy (same mesh)
        +-- MultiAppNearestNodeTransfer  interpolate to nearest nodes
        +-- MultiAppProjectionTransfer   L2 projection onto target mesh
        +-- MultiAppShapeEvaluationTransfer  shape function evaluation
        +-- MultiAppUserObjectTransfer   query a UserObject value
        +-- MultiAppPostprocessorTransfer  scalar postprocessor value
        +-- MultiAppReporterTransfer     Reporter key-value data
        +-- MultiAppDofCopyTransfer      raw DOF vector copy
        +-- MultiAppGeneralFieldTransfer base for geometric field transfers
```

### Picard (Fixed-Point) Iteration

When a MultiApp simulation requires iterative coupling between parent and sub-apps within a single timestep, `PicardSolve` manages the loop:

```
PicardSolve::solve()
    |
    +-- iteration = 0
    +-- while not converged:
    |       |
    |       +-- Execute "TO_MULTIAPP" Transfers
    |       |     parent field -> sub-app field
    |       |
    |       +-- MultiApp::solveStep(dt, target_time)
    |       |     sub-apps solve their physics
    |       |
    |       +-- Execute "FROM_MULTIAPP" Transfers
    |       |     sub-app field -> parent field
    |       |
    |       +-- parent FEProblemBase::solve()
    |       |     (parent solves with updated sub-app data)
    |       |
    |       +-- check convergence (postprocessor change, etc.)
    |       +-- iteration++
    |
    +-- return converged
```

### Parent-Sub App Communication Diagram

```
+------------------------------------+
|  Parent App                        |
|  FEProblemBase                     |
|    NonlinearSystem                 |
|      solution vector U_parent      |
|    AuxiliarySystem                 |
|      aux field T_parent            |
|    MultiApp "thermal_micro"        |
+------------------------------------+
          |          ^
          |          |
 Transfer: TO        | Transfer: FROM
 (T_parent           | (T_micro
  -> T_micro)        |  -> T_parent)
          |          |
          v          |
+------------------------------------+
|  Sub-App[0]                        |
|  FEProblemBase                     |
|    NonlinearSystem                 |
|      solution vector T_micro       |
|    AuxiliarySystem                 |
|      output field flux_micro       |
+------------------------------------+
+------------------------------------+
|  Sub-App[1]   (similar)            |
+------------------------------------+
+------------------------------------+
|  Sub-App[2]   (similar)            |
+------------------------------------+
```

---

## 10. The Finite Volume Subsystem

### FE vs FV: Conceptual Difference

| Aspect               | Finite Element (FE)                  | Finite Volume (FV)                     |
|----------------------|--------------------------------------|----------------------------------------|
| DOF location         | Nodes (Lagrange) or element interior  | Cell centroid (constant per cell)      |
| Integration domain   | Element volume integral               | Face flux + cell volume                |
| Residual form        | Weak form: integral of test * residual | Strong form: div theorem on cell faces |
| Class base           | `KernelBase`                         | `FVKernel`                             |
| Loop type            | `ComputeResidualThread` (elem loop)  | `ComputeFVFluxThread` (face loop)      |
| Variable type        | `MooseVariableFE<T>`                 | `MooseVariableFV<T>`                   |

### FV Kernel Hierarchy

```
ResidualObject
    |
    +-- FVKernel  (fvkernels/FVKernel.h)
    |   BlockRestrictable, ADFunctorInterface, FVRelationshipManagerInterface
    |
    +-- FVFluxKernel  (fvkernels/FVFluxKernel.h)
    |   NeighborCoupleable, TwoMaterialPropertyInterface
    |   computeQpResidual() -> ADReal   [pure virtual]
    |   Represents: surface integral term (divergence theorem applied)
    |   Loop: over all internal faces + boundary faces
    |   Inputs: _u_elem, _u_neighbor, _face_info, _normal
    |   Examples: FVDiffusion, FVAdvection, FVMatAdvection
    |
    +-- FVElementalKernel  (fvkernels/FVElementalKernel.h)
        computeQpResidual() -> ADReal   [pure virtual]
        Represents: volume integral term (source, reaction, time derivative)
        Loop: over all elements (cell-centered)
        Examples: FVTimeKernel, FVBodyForce, FVReaction
```

### Face Loop vs Element Loop

```
FE Residual Loop (ComputeResidualThread):
    for each element e in local mesh partition:
        Assembly::reinit(e)               // shape fns at quad pts
        for each active Kernel k:
            k->computeResidual()          // integrate over element volume

FV Residual Loop (ComputeFVFluxThread):
    // Flux kernels: face-centered loop
    for each FaceInfo fi in local mesh:
        // fi knows: elem side, neighbor side, normal, area, centroids
        for each active FVFluxKernel k:
            k->computeResidual(fi)
                // evaluates at fi._qp=0 (single face quadrature pt)
                // uses reconstructed cell-centered values on both sides
                flux = k->computeQpResidual()   // ADReal
                // scatter: residual[elem_dof] += flux * area
                //          residual[neighbor_dof] -= flux * area (conservation)

    // Elemental kernels: cell-centered loop (normal element loop)
    for each element e:
        for each active FVElementalKernel k:
            k->computeResidual()
```

### Ghost Cells for FV

For a parallel-distributed mesh, the `FVFluxKernel` at a partition boundary needs the DOF value from the neighboring cell, which lives on a different MPI rank. MOOSE uses `RelationshipManager` to tell libMesh which elements must be ghosted:

```
FVKernel::FVRelationshipManagerInterface
    |
    +-- requests AlgebraicGhostingFunctor
            (ensures neighboring cell values are communicated
             via libMesh's ghosting mechanism before each solve)

Ghost layer:
  Rank 0 owns elements 0..499
  Rank 1 owns elements 500..999
  At the partition boundary:
    Rank 0 holds ghost copy of element 500 (from rank 1)
    Rank 1 holds ghost copy of element 499 (from rank 0)
  FVFluxKernel on rank 0 can read U[500] without MPI call —
  the ghost has already been communicated during vector scatter.
```

---

## 11. The Mesh Generator Pipeline

### What Are Mesh Generators?

A `MeshGenerator` is a `MooseObject` with a single pure-virtual method:

```cpp
virtual std::unique_ptr<MeshBase> generate() = 0;
```

Each generator receives zero or more input meshes (by name, via the `input` or `inputs` parameter) and returns a new (possibly modified) mesh. They form a **directed acyclic graph (DAG)** driven by their `input` parameter references.

### DAG Execution

```
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator  # leaf: no inputs
    dim = 2
    nx = 10, ny = 10
  []
  [refine]
    type = RefineBlockGenerator     # depends on "gmg"
    input = gmg
    block = 0
    refinement = 2
  []
  [sideset]
    type = ParsedGenerateSideset    # depends on "refine"
    input = refine
    combinatorial_geometry = 'x > 0.9'
    new_sideset_name = right
  []
[]

DAG:
  GeneratedMeshGenerator("gmg")
      |
      v
  RefineBlockGenerator("refine")   input="gmg"
      |
      v
  ParsedGenerateSideset("sideset") input="refine"
      |
      v
  [final mesh used by FEProblemBase]
```

`MeshGeneratorSystem::executeMeshGenerators()` (owned by `MooseApp`):
1. Topologically sorts generators using the input dependency graph.
2. Executes them in order, passing the output of each to its consumers.
3. The final generator's output becomes `MooseMesh::_mesh`.

### MeshGenerator Subclass Categories

```
PRIMITIVE GENERATORS (create mesh from scratch):
    GeneratedMeshGenerator    regular quad/hex grid
    AnnularMeshGenerator      annular (ring) regions
    CartesianMeshGenerator    arbitrary Cartesian blocks
    FileMeshGenerator         load from Exodus/Gmsh/etc.
    DistributedRectilinearMeshGenerator  parallel generation

MODIFICATION GENERATORS (take input, return modified copy):
    RefineBlockGenerator      uniform h-refinement
    CoarsenBlockGenerator     coarsen selected blocks
    BlockDeletionGenerator    remove elements
    BoundaryLayerSubdomainGenerator
    AdvancedExtruderGenerator 2D -> 3D extrusion with twisting
    BreakMeshByBlockGenerator create element interfaces at subdomain boundaries
    CombinerGenerator         merge multiple meshes

SIDESET / NODESET GENERATORS:
    AllSideSetsByNormalsGenerator
    BoundingBoxNodeSetGenerator
    ParsedGenerateSideset
    ExtraNodesetGenerator

TRANSFORMATION GENERATORS:
    TransformGenerator        translate/rotate/scale
    SymmetryTransformGenerator mirror
    FlipSidesetGenerator

UTILITY:
    AddMetaDataGenerator      attach mesh metadata
    LowerDBlockFromSidesetGenerator  create lower-dimensional elements on sidesets
    MeshRepairGenerator       fix mesh quality issues
    RenameBlockGenerator, RenameBoundaryGenerator
```

---

## 12. The Input Parameter System

### What Is `InputParameters`?

`InputParameters` (declared in `framework/include/utils/InputParameters.h`, inheriting from `libMesh::Parameters`) is a typed dictionary: every parameter has a name, a C++ type, a description, optional default, and metadata about how it was set (by user, by default, via command line).

It is the universal interface between:
- The HIT input file parser (which populates parameters from text)
- The Factory (which passes parameters to constructors)
- MooseObjects (which read parameters in their constructors via `getParam<T>()`)

### `validParams()` Chain

Every MOOSE class that can appear in an input file implements the static method:

```cpp
static InputParameters validParams();
```

This method builds and returns the complete parameter specification for the class, including all parameters inherited from base classes:

```cpp
// Simplified from Kernel.C
InputParameters
Kernel::validParams()
{
    InputParameters params = KernelBase::validParams();  // call parent
    // add this class's own parameters:
    params.addParam<bool>("use_displaced_mesh", false, "Whether to use displaced mesh");
    return params;
}

InputParameters
KernelBase::validParams()
{
    InputParameters params = ResidualObject::validParams();
    params.addParam<std::vector<AuxVariableName>>("save_in", {}, "...");
    params.addParam<std::vector<AuxVariableName>>("diag_save_in", {}, "...");
    return params;
}

InputParameters
ResidualObject::validParams()
{
    InputParameters params = MooseObject::validParams();
    params.addRequiredParam<NonlinearVariableName>("variable", "The variable this kernel operates on");
    return params;
}

InputParameters
MooseObject::validParams()
{
    InputParameters params = MooseBase::validParams();
    params.addParam<bool>("enable", true, "Enable/disable this object");
    params.addParam<std::vector<std::string>>("control_tags", {}, "...");
    return params;
}
```

The full chain: `MooseBase::validParams()` -> `MooseObject::validParams()` -> `ResidualObject::validParams()` -> `KernelBase::validParams()` -> `Kernel::validParams()` -> `Diffusion::validParams()`.

### Parameter Types

`InputParameters` supports a rich set of types via `addParam<T>()`:

```
Scalars:   bool, int, unsigned int, Real, std::string, MooseEnum
Vectors:   std::vector<Real>, std::vector<int>, std::vector<std::string>
MOOSE types: NonlinearVariableName, AuxVariableName, MaterialPropertyName,
             FunctionName, UserObjectName, PostprocessorName, MeshGeneratorName,
             BoundaryName, SubdomainName, FileName, DataFileName
Special:   Point, RealVectorValue, RealTensorValue
```

Special parameter types enable the framework to validate cross-references (e.g., a `UserObjectName` must match a registered UserObject) and to provide meaningful error messages.

### How Parameters Reach Member Variables

```
// In problem.i:
[Kernels]
  [./diff]
    type = Diffusion
    variable = u           # <- this string is stored in InputParameters
  [../]
[]

// In Diffusion.C:
Diffusion::Diffusion(const InputParameters & params)
  : Kernel(params)          // passes params up the chain
{
  // Kernel constructor:
  //   _var = _sys.getVariable<MooseVariableFE<Real>>(
  //               params.get<NonlinearVariableName>("variable"));
  //   _u = &_var.sln();
  //   _grad_u = &_var.gradSln();
}

Real Diffusion::computeQpResidual() {
    return _grad_u[_qp] * _grad_test[_qp];  // -div(grad(u)) weak form
}
```

### Command-Line Override

Parameters flagged with `addParam<T>("key", val, "desc", CommandLineMetadata)` can be overridden at runtime:

```
./myapp-opt -i problem.i Mesh/nx=20 Executioner/num_steps=100
```

### InputParameterWarehouse

The `Factory` stores all `InputParameters` objects created for live MooseObjects in an `InputParameterWarehouse`. This ensures that:
1. Exactly one `InputParameters` object lives for each MooseObject instance.
2. The `Factory::currentlyConstructing()` method can return parameters to child objects during construction (used for sub-object creation patterns).

---

## 13. Parallel Architecture

### MPI: Mesh Partitioning

MOOSE uses libMesh's `DistributedMesh` (or `ReplicatedMesh` for small problems). For parallel runs, the mesh is partitioned — typically with Parmetis or a space-filling curve — so that each MPI rank owns a contiguous subdomain:

```
Global Mesh (1000 elements)         After partitioning (4 ranks):
                                    Rank 0: elements   0..249  (+ ghost layer)
+----+----+----+----+               Rank 1: elements 250..499  (+ ghost layer)
| 0  | 1  | 2  | 3  |               Rank 2: elements 500..749  (+ ghost layer)
+----+----+----+----+   ------>     Rank 3: elements 750..999  (+ ghost layer)
| 4  | 5  | 6  | 7  |
+----+----+----+----+
```

Each rank only holds:
- Its own elements (full data)
- **Ghost elements**: elements from neighboring partitions needed for stencil computations (flux kernels, DG kernels, geometric searches)

### DOF Numbering

libMesh's `DofMap` assigns globally unique DOF indices across ranks. For a problem with variable `u` on LAGRANGE_FIRST:
- Each node has one DOF for `u`.
- DOF indices are contiguous within each rank, allowing efficient PETSc vector assembly.
- The `PetscVector<Number>` is partitioned: each rank owns rows `first..last`.

### Ghost Communication

Before each residual/Jacobian evaluation:
1. `FEProblemBase::updateSolution(X)` scatters the solution from the PETSc global vector to local ghosted vectors.
2. `libMesh::NumericVector::localize_to_one()` or `close()` synchronizes ghost values across MPI ranks.
3. Kernels and FV flux kernels can then read neighbor values from their local ghost elements without MPI calls during the loop.

### RelationshipManagers

`RelationshipManager` objects (in `framework/include/relationshipmanagers/`) tell libMesh how many layers of ghost elements each system needs:

```
ElementSideNeighborLayers(n=1):   ghost 1 layer of face-neighbors
                                  (sufficient for standard FV and DG)
ElementPointNeighborLayers(n=2):  ghost 2 layers via node connectivity
                                  (needed for some reconstruction stencils)
```

These are registered by MooseObjects in their `validParams()` via `addRelationshipManager()`, and processed before mesh distribution.

### Thread Parallelism: THREAD_ID

Within each MPI rank, MOOSE uses OpenMP or Intel TBB (selected at build time via libMesh's threading backend) for element-loop parallelism:

```
#define THREAD_ID unsigned int

// Each thread has its own:
// - Assembly object  (shape function values, local dense vectors)
// - MaterialData     (quadrature-point property storage)
// - Kernel warehouse copy (thread-safe: no shared mutable state in kernels)

// MOOSE allocates numThreads copies of thread-sensitive objects:
std::vector<std::unique_ptr<Assembly>> _assembly;   // indexed by THREAD_ID
```

Thread safety is maintained by:
- Each thread working on a disjoint set of elements.
- Thread-local copies of all mutable per-element state.
- Synchronized reduction at the end of the threaded loop (Assembly scatter phase).

### ParallelObject

Most MOOSE infrastructure objects inherit from `libMesh::ParallelObject`, which carries a `libMesh::Parallel::Communicator`. This communicator may be:
- `MPI_COMM_WORLD` for the main app
- A sub-communicator for a MultiApp sub-app instance

This allows sub-apps to safely call collective MPI operations without colliding with the parent app's communicator.

### Typical Parallel Solve Timing

```
MPI barrier (all ranks sync)
    |
    v
FEProblemBase::computeResidual()
    |
    +-- [rank-local] ThreadedElementLoop (OpenMP threads per rank)
    |       each thread: Assembly::reinit + Material + Kernel loops
    +-- [rank-local] Assembly::addResidualLocal()  -> PETSc Vec (local rows)
    +-- PETSc VecAssemblyBegin/End()     <- MPI collective
    |       gathers boundary contributions across ranks
    v
PETSc SNES evaluates ||F|| (requires VecNorm -> MPI_Allreduce)
    |
    v
PETSc KSP linear solve (GMRES/CG)
    |
    +-- mat-vec product A*x:
    |       MatMult (PETSc) -> MOOSE-filled sparse matrix
    |       requires ghost DOF communication: VecGhostUpdateBegin/End()
    +-- preconditioner apply
    v
Newton update x += dx   (local, no communication needed)
```

---

## 14. Communication Flow Diagrams

### Diagram A: Initialization — Parser to Objects

This diagram shows how the input file travels from text to live C++ objects during application startup.

```
                 INITIALIZATION COMMUNICATION FLOW
                 ===================================

  problem.i (text file)
       |
       | MooseApp::run()
       | -> MooseApp::setupOptions()
       | -> MooseApp::runInputFile()
       |
       v
  Parser::parse(filename)
  Moose::Builder::build(root_node)          // iterates HIT tree
       |
       | For each HIT block node:
       |   Builder::walk(fullpath, nodename, hit::Node*)
       |       |
       |       +-- Syntax::isActionRequired(path)?
       |       |       yes: auto-create action from syntax
       |       |
       |       +-- ActionFactory::create(action_type, name, params)
       |       |       -> new AddKernelAction(params)  [example]
       |       |
       |       +-- ActionWarehouse::addActionBlock(action)
       |
       v
  ActionWarehouse populated:
    task "setup_mesh":       [SetupMeshAction, AddMeshGeneratorAction, ...]
    task "create_mesh":      [CreateMeshAction]
    task "add_variable":     [AddVariableAction("u"), AddVariableAction("v")]
    task "add_kernel":       [AddKernelAction("diff"), AddKernelAction("rxn")]
    task "add_bc":           [AddBCAction("left_bc"), AddBCAction("right_bc")]
    task "add_material":     [AddMaterialAction("mat1")]
    task "init_problem":     [InitProblemAction]
    task "setup_executioner":[CreateExecutionerAction]
       |
       | ActionWarehouse::executeAllTasks()
       v
  For each task in dependency order:
       |
       | task = "add_kernel"
       v
  AddKernelAction::timedAct()
  AddKernelAction::act()
       |
       +-- type = getParam<std::string>("type")  // "Diffusion"
       +-- params = _factory.getValidParams("Diffusion")
       |                  |
       |                  +-- Registry lookup: "Diffusion" -> RegistryEntry<Diffusion>
       |                  +-- RegistryEntry<Diffusion>::buildParameters()
       |                      -> Diffusion::validParams()
       |                      <- InputParameters{variable, enable, ...}
       +-- _moose_object_pars.applyInputParameters(params, hit_node)
       |       // fills: params["variable"] = "u"  [from input file]
       +-- _problem->addKernel("Diffusion", "diff", params, tid=0)
               |
               v
         FEProblemBase::addKernel()
               |
               +-- kernel = _factory.create<KernelBase>("Diffusion", "diff", params)
               |               |
               |               +-- InputParameterWarehouse: store params
               |               +-- RegistryEntry<Diffusion>::buildShared(params)
               |               |       -> std::make_shared<Diffusion>(params)
               |               |              Diffusion constructor:
               |               |              _var = _sys.getVariable<>(_u_var_name)
               |               +-- Factory::finalize(...)
               |               <- std::shared_ptr<Diffusion>
               +-- _nl_sys[0].addKernel(kernel, tid=0)
                       // stored in MooseObjectTagWarehouse<KernelBase>
```

### Diagram B: Solve Loop — Executioner to Quadrature Point

This diagram shows the call chain during a single Newton iteration of a `Transient` solve.

```
            SOLVE LOOP COMMUNICATION FLOW
            ==============================

  TransientBase::execute()
       |  (timestep loop)
       v
  TransientBase::takeStep(dt)
       |
       v
  FixedPointSolve::solve()       [optional: Picard loop for MultiApp]
       |
       v
  FEProblemSolve::solve()        [FEProblemSolve.h]
       |
       v
  FEProblemBase::solve()         [problems/FEProblemBase.h]
       |
       +-- execMultiApps(EXEC_TIMESTEP_BEGIN)   [sub-apps advance first if needed]
       +-- _nl[0]->solve()
       v
  NonlinearSystem::solve()       [systems/NonlinearSystem.h]
       |
       | (delegates to libMesh)
       v
  libMesh::NonlinearImplicitSystem::solve()
       |
       v
  PETSc SNESSolve()
       |
       +-- Newton iteration 0:
       |     SNESComputeFunction()
       |       -> ComputeResidualFunctor::operator()()
       |             -> FEProblemBase::computeResidualInternal()
       |                   |
       |                   +-- updateSolution(X)    [ghost comm: MPI]
       |                   +-- reinitMaterials()
       |                   +-- Parallel_reduce<ComputeResidualThread>
       |                         |
       |                         | [Thread 0..N, each handles subset of elements]
       |                         v
       |                   ComputeResidualThread::onElement(Elem* e)
       |                         |
       |                         +-- Assembly::reinit(e)
       |                         |       libMesh: compute phi, dphi at Gauss pts
       |                         |       update _JxW[], _coord[], _q_point[]
       |                         |
       |                         +-- MaterialWarehouse: computeMaterials(e, tid)
       |                         |   DependencyResolver orders materials
       |                         |   for each Material m:
       |                         |       m->computeProperties()
       |                         |           for qp=0..nQP:
       |                         |               m->computeQpProperties()  <-- USER CODE
       |                         |               // fills _prop_name[qp]
       |                         |
       |                         +-- for each Kernel k active on e:
       |                         |       k->computeResidual()
       |                         |           for i in test_dofs:
       |                         |               for qp=0..nQP:
       |                         |                   _re[i] += _JxW[qp]*_coord[qp]
       |                         |                             *_test[i][qp]
       |                         |                             *k->computeQpResidual()  <-- USER CODE
       |                         |                   // computeQpResidual uses:
       |                         |                   //   _u[qp], _grad_u[qp]   (solution)
       |                         |                   //   _phi[j][qp]           (trial fn)
       |                         |                   //   _JxW[qp]              (det*weight)
       |                         |                   //   MaterialProperty[qp]  (from Material)
       |                         |
       |                         +-- Assembly::addResidualLocal(tid)
       |                               scatter _re -> global PETSc Vec F
       |
       |     SNESComputeJacobian()
       |       -> (similar path via ComputeJacobianThread)
       |
       +-- Newton iteration 1: (repeat until ||F|| < tol)
       v
  NonlinearSystem::converged()
       |
       v
  FEProblemBase::computeAuxiliaryKernels()
  FEProblemBase::computeUserObjects(EXEC_TIMESTEP_END)
  FEProblemBase::computePostprocessors(EXEC_TIMESTEP_END)
       |
       v
  OutputWarehouse::outputStep()    -> Exodus::output(), CSV::output(), ...
```

### Diagram C: MultiApp Coupled Solve

This diagram shows communication during a Picard-iterated parent-sub-app coupled solve.

```
          MULTIAPP COUPLED SOLVE COMMUNICATION FLOW
          ===========================================

  [Parent App - rank 0..N]          [Sub-App 0 - rank 0..k]  [Sub-App 1 - rank k+1..N]

  TransientBase::execute()
       |
       v
  TransientBase::takeStep(dt)
       |
       v
  PicardSolve::solve()              (FixedPointSolve subclass)
       |
       +-- iter=0
       |
       |   STEP 1: Push data TO sub-apps
       |   FixedPointSolve::execMultiAppTransfers(EXEC_MULTIAPP_FIXED_POINT_BEGIN)
       |       |
       |       | MultiAppCopyTransfer::execute()
       |       |   direction = TO_MULTIAPP
       |       |   for each sub-app i:
       |       |       from: parent FEProblemBase variable "T_macro"
       |       |       to:   sub_app[i] FEProblemBase variable "T_in"
       |       |       [MPI communication to sub-app ranks if needed]
       |       |                                |                        |
       |       |                          [receive T_macro]       [receive T_macro]
       |       |                          T_in = T_macro          T_in = T_macro
       |
       |   STEP 2: Solve sub-apps
       |   MultiApp::solveStep(dt, target_time)
       |                                |                        |
       |                         sub_app[0]               sub_app[1]
       |                         Transient::execute()     Transient::execute()
       |                         PETSc SNES solve         PETSc SNES solve
       |                         (uses T_in as BC)        (uses T_in as BC)
       |                         produces flux_micro[0]   produces flux_micro[1]
       |
       |   STEP 3: Pull data FROM sub-apps
       |   FixedPointSolve::execMultiAppTransfers(EXEC_MULTIAPP_FIXED_POINT_END)
       |       |
       |       | MultiAppUserObjectTransfer::execute()
       |       |   direction = FROM_MULTIAPP
       |       |   for each sub-app i:
       |       |       query: sub_app[i] UserObject "avg_flux"
       |       |       ->getValue() [MPI comm back to parent ranks if needed]
       |       |       write to parent AuxVariable "flux_from_micro[i]"
       |       |
       |
       |   STEP 4: Solve parent
       |   FEProblemSolve::solve()
       |       FEProblemBase::solve()
       |           NonlinearSystem::solve()
       |           PETSc SNES -> computeQpResidual() (uses flux_from_micro)
       |
       |   STEP 5: Check convergence
       |   PicardSolve::checkConvergence()
       |       compute ||delta_U|| / ||U||   (between current and previous Picard iter)
       |       if < picard_abs_tol: break
       |
       +-- iter=1
       |   (repeat steps 1-5)
       +-- iter=2 (converged)
       |
       v
  TransientBase::endStep()
       |
       v
  OutputWarehouse::outputStep()
      Exodus output, Postprocessors, etc.
```

---

## Appendix A: Key Data Structures

### MaterialProperty<T>

```
MaterialProperty<T>  (framework/include/materials/MaterialProperty.h)
    std::vector<T> _values;   // indexed by quadrature point [0..nQP-1]

Access pattern:
    _prop[_qp]     // in Kernel or Material: read/write at current qp
    _prop.size()   // == nQP at current element
```

### FaceInfo

Used by FV kernels to describe a mesh face:

```
FaceInfo  (framework/include/mesh/FaceInfo.h)
    const Elem * _elem;          // element on elem side
    const Elem * _neighbor;      // element on neighbor side (null at boundary)
    Point _face_centroid;        // face center point
    RealVectorValue _normal;     // unit normal pointing from elem to neighbor
    Real _face_area;             // area of the face
    Real _elem_volume;           // volume of elem
    Real _neighbor_volume;       // volume of neighbor
    Point _elem_centroid;
    Point _neighbor_centroid;
```

### TheWarehouse

`TheWarehouse` (framework/include/base/TheWarehouse.h) is a tag-based object store used to query MooseObjects by attributes:

```cpp
// Example query: get all Kernels active on subdomain 0, thread 0,
//                that should execute on EXEC_LINEAR
auto & kernels = _warehouse.query()
    .condition<AttribSubdomains>(0)
    .condition<AttribThread>(0)
    .condition<AttribExecFlags>(EXEC_LINEAR)
    .queryInto(results);
```

Attributes include: `AttribThread`, `AttribSubdomains`, `AttribBoundaries`, `AttribExecFlags`, `AttribMatrixTags`, `AttribVectorTags`, `AttribSystem`, `AttribName`.

---

## Appendix B: Common Class Lookup Table

| Task | Base Class | Virtual Method to Override | Register Macro |
|------|-----------|---------------------------|----------------|
| Volume PDE term (FE) | `Kernel` or `ADKernel` | `computeQpResidual()` | `registerMooseObject` |
| Boundary condition (Neumann) | `IntegratedBC` | `computeQpResidual()` | `registerMooseObject` |
| Boundary condition (Dirichlet) | `DirichletBC` | `computeQpValue()` | `registerMooseObject` |
| Face flux (FV) | `FVFluxKernel` | `computeQpResidual()` | `registerMooseObject` |
| Cell source (FV) | `FVElementalKernel` | `computeQpResidual()` | `registerMooseObject` |
| Material property | `Material` | `computeQpProperties()` | `registerMooseObject` |
| Lazy material property | `FunctorMaterial` | `addFunctorProperty<T>()` | `registerMooseObject` |
| Scalar reduction | `GeneralUserObject` | `initialize()/execute()/finalize()` | `registerMooseObject` |
| Element integral | `ElementIntegralUserObject` | `computeQpIntegral()` | `registerMooseObject` |
| Postprocessor value | subclass of `Postprocessor` | `getValue()` | `registerMooseObject` |
| Auxiliary field fill | `AuxKernel` | `computeValue()` | `registerMooseObject` |
| Mesh construction | `MeshGenerator` | `generate()` | `registerMooseObject` |
| Simulation driver | `Executioner` | `execute()` | `registerMooseObject` |
| Output writer | `Output` (or `FileOutput`) | `output()` | `registerMooseObject` |
| Sub-app coupling | `MultiApp` | `solveStep()` | `registerMooseObject` |
| Data transfer | `MultiAppTransfer` | `execute()` | `registerMooseObject` |
| Setup task | `Action` | `act()` | `registerMooseAction` |

---

## Appendix C: Glossary

| Term | Meaning |
|------|---------|
| **HIT** | Hierarchical Input Text — MOOSE's input file format. Blocks denoted `[Name]...[../]`, parameters as `key = value`. |
| **DOF** | Degree of freedom — one scalar unknown in the linear algebra system. A LAGRANGE_FIRST variable has one DOF per node. |
| **qp / `_qp`** | Quadrature point index. Gauss integration points within an element. Material properties and solution values are indexed by `_qp`. |
| **Residual** | The vector F(u) that is zero at the solution. Newton's method iterates to drive ||F|| to zero. |
| **Jacobian** | The matrix dF/du. Used by Newton iterations and the linear solver. |
| **SNES** | Scalable Nonlinear Equation Solver — PETSc's nonlinear solver framework. |
| **KSP** | Krylov Subspace Method solver — PETSc's linear solver framework (GMRES, CG, ...). |
| **PC** | Preconditioner — transforms the linear system to improve solver convergence. |
| **AD** | Automatic Differentiation — MOOSE computes Jacobian entries automatically from AD-enabled kernels using `DualNumber` forward-mode AD. |
| **EXEC flag** | Execution flag: `EXEC_INITIAL`, `EXEC_TIMESTEP_BEGIN`, `EXEC_LINEAR`, `EXEC_NONLINEAR`, `EXEC_TIMESTEP_END`, etc. Controls when UserObjects/outputs fire. |
| **Warehouse** | A `MooseObjectWarehouse` — typed container for MooseObjects, queryable by subdomain, boundary, thread, and exec flag. |
| **Ghosting** | The process of replicating elements from neighboring MPI ranks so that stencil computations can proceed without explicit MPI calls during loops. |
| **Stateful property** | A `MaterialProperty` that retains values from previous timesteps (`_old`, `_older`). Requires registration at construction so storage can be allocated. |
| **Functor** | A `Moose::FunctorBase<T>` — callable object that evaluates a field at a given spatial/temporal argument. Used by FunctorMaterial and FV systems. |
| **Picard iteration** | Fixed-point iteration coupling parent app and sub-apps. Alternates between solving sub-apps and parent until convergence. |
| **Task** | A named setup step in the ActionWarehouse ordering system (e.g., `"add_kernel"`, `"init_problem"`). |
| **Label** | The string identifier for an application or module in the Registry (e.g., `"MooseApp"`, `"SolidMechanicsApp"`). |
| **THREAD_ID** | `unsigned int` index for the current OpenMP/TBB thread. Used to index thread-local copies of Assembly, MaterialData, etc. |
