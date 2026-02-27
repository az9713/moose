# Zero to Hero: A Self-Study Plan for MOOSE Framework Development

**Audience:** Self-directed learners with C/C++/Java experience who are new to scientific
computing. No prior finite-element, PDE, or solver knowledge assumed.

**Goal:** Bring you from zero MOOSE knowledge to the point where you can design, implement,
test, and document your own physics module — understanding both the theory behind what you
are writing and the framework machinery that runs it.

**How to use this plan:** Each phase has daily reading assignments, theory explanations, and
hands-on exercises. Do the theory and the code reading on the same day. Do not rush — the
goal is deep understanding, not speed. Every success check is a gate: do not proceed until
you can pass it.

---

## Table of Contents

- [Phase 0: Prerequisites — Setting the Stage](#phase-0-prerequisites--setting-the-stage)
- [Phase 1: User Foundations — Running Simulations](#phase-1-user-foundations--running-simulations)
- [Phase 2: Understanding the Framework — How MOOSE Works](#phase-2-understanding-the-framework--how-moose-works)
- [Phase 3: First C++ Object — Writing Your First Kernel](#phase-3-first-c-object--writing-your-first-kernel)
- [Phase 4: Materials and BCs — Building a Complete App](#phase-4-materials-and-bcs--building-a-complete-app)
- [Phase 5: Advanced Topics — Mastering the Framework](#phase-5-advanced-topics--mastering-the-framework)
- [Phase 6: Building a Physics Module — Your Own Module](#phase-6-building-a-physics-module--your-own-module)
- [Phase 7: Production Skills — Going Pro](#phase-7-production-skills--going-pro)
- [Appendix: Quick Reference Card](#appendix-quick-reference-card)

---

## Phase 0: Prerequisites — Setting the Stage

**Objective:** Arrive at Day 1 of Phase 1 with the C++ skills and mathematical vocabulary
you will need. Install the framework and confirm it builds.

---

### C++ Concepts to Review

MOOSE is a framework of base classes. You write new physics by subclassing and overriding.
If any of the following patterns are unfamiliar, spend a few hours with a C++ reference
before reading on.

#### Inheritance and Virtual Functions

```cpp
// MOOSE declares the interface in a base class with pure virtuals:
class Kernel {
  virtual Real computeQpResidual() = 0;  // = 0 means pure virtual — you MUST override
  virtual Real computeQpJacobian() { return 0; }  // default: optional override
};

// You provide the physics:
class MyKernel : public Kernel {
  Real computeQpResidual() override;  // override keyword catches signature mismatches
};
```

MOOSE calls your overrides during the solve; you never call them yourself. The pure virtual
`= 0` pattern is how MOOSE enforces the "plugin" contract.

**File to read:** `framework/include/kernels/Kernel.h` — notice `virtual Real computeQpResidual() = 0` on line 46.

#### Templates

Templates appear everywhere in MOOSE's type-generic code. You will encounter them most
often when declaring material properties:

```cpp
MaterialProperty<Real> & _conductivity = declareProperty<Real>("thermal_conductivity");
MaterialProperty<RealTensorValue> & _stress = declareProperty<RealTensorValue>("stress");
```

You do not need to write your own templates for most MOOSE development, but you must be
comfortable reading template syntax to understand error messages.

The `GenericKernel<is_ad>` pattern (used in `framework/include/kernels/Reaction.h`) is
the main example you will encounter: a single class body parameterized over whether it
uses automatic differentiation.

#### Smart Pointers

`std::unique_ptr<T>` is used heavily in `MeshGenerator::generate()`. When a method
signature returns `unique_ptr`, you must move it, not copy it:

```cpp
std::unique_ptr<MeshBase> generate() override {
  auto mesh = buildMeshBaseObject();
  return mesh;  // transfers ownership — no copy
}
```

`std::shared_ptr<T>` appears in MultiApp ownership inside MOOSE internals. You rarely
create these yourself.

**File to read:** `framework/include/meshgenerators/GeneratedMeshGenerator.h` — the
`generate()` override returns `std::unique_ptr<MeshBase>`.

#### RAII (Resource Acquisition Is Initialization)

MOOSE objects live in warehouses owned by `FEProblemBase`. They are created via
`Factory::create<T>()` and stored in `std::shared_ptr`. When the problem is destroyed,
smart-pointer destructors clean up automatically. You never call `delete` on MOOSE objects.

#### const References

Almost everything MOOSE gives your code arrives as a `const` reference:

```cpp
const VariableValue & _u;         // solution values at quadrature points
const VariableGradient & _grad_u; // gradients at quadrature points
```

The `&` means no copy is made. `const` means you may read but not write. This is a
performance contract enforced by the type system — not just a style preference.

---

### Math Background: PDEs, Weak Forms, and the Galerkin Method

This section explains the mathematics that MOOSE is built on. Read it carefully. Every
concept maps directly to code you will write.

#### What is a PDE?

A Partial Differential Equation (PDE) is an equation relating an unknown function
`u(x, y, z, t)` to its partial derivatives. The three archetypal examples:

**Steady diffusion (Laplace / Poisson equation):**
```
  -∇·(D ∇u) = f    in Ω
```
`u` is temperature (or concentration, pressure, potential), `D` is diffusivity, `f` is a
source. `∇` is the gradient operator, `∇·` is the divergence. This is the equation that
`simple_diffusion.i` solves when `D = 1`, `f = 0`.

**Transient diffusion (heat equation):**
```
  ∂u/∂t - ∇·(D ∇u) = 0    in Ω × [0, T]
```
Add a time derivative. The `TimeDerivative` kernel contributes `∂u/∂t`. Solved with
`type = Transient` in `[Executioner]`.

**Wave equation (second-order in time):**
```
  ∂²u/∂t² - c² ∇²u = 0
```
Hyperbolic PDE. Less common in MOOSE's primary domain (parabolic and elliptic problems),
but solid mechanics uses it in dynamic simulations.

**Nonlinear diffusion (reaction-diffusion):**
```
  -∇·(D(u) ∇u) + R(u) = 0
```
When `D` or `R` depends on `u`, the system is nonlinear. Newton's method handles this in
MOOSE — the Jacobian you compute in `computeQpJacobian()` drives the Newton update.

#### What is a Weak Form? Why Do We Need It?

Computers cannot represent continuous functions. They work with numbers at discrete points.
The weak form is the bridge.

Start with the strong form: find `u` such that `-∇·(∇u) = 0` everywhere in `Ω`. Multiply
both sides by an arbitrary smooth "test function" `ψᵢ` and integrate over the domain:

```
  ∫_Ω (-∇·∇u) ψᵢ dΩ = 0
```

Apply integration by parts (Green's theorem in multiple dimensions). The divergence theorem
converts the volume integral of `∇·(...)` into a surface integral, lowering the order of
the derivative requirement:

```
  ∫_Ω ∇u · ∇ψᵢ dΩ - ∮_∂Ω (∇u · n) ψᵢ dS = 0
```

The boundary integral is where Neumann BCs enter naturally. The **weak form** is:

```
  ∫_Ω ∇u · ∇ψᵢ dΩ = ∮_∂Ω (∇u · n) ψᵢ dS    for all i
```

**Why we need it:**
1. The weak form requires only first derivatives of both `u` and `ψᵢ`, while the strong
   form requires second derivatives of `u`. This means we can use simpler polynomial basis
   functions.
2. Neumann boundary conditions (prescribing the flux `∇u · n`) appear naturally in the
   right-hand side — you do not need to impose them explicitly.
3. The mathematical theory of well-posedness (Lax-Milgram theorem) applies to weak forms.

**Connection to code:** `computeQpResidual()` returns the integrand of the left-hand side
at one quadrature point, for one test function index `_i`. In `Diffusion.C`:

```cpp
Real Diffusion::computeQpResidual() {
  return _grad_u[_qp] * _grad_test[_i][_qp];  // ∇u · ∇ψᵢ  at quadrature point _qp
}
```

MOOSE loops over all `_qp` and all `_i`, multiplies by the weight `_JxW[_qp]`, and
accumulates into the global residual vector. You only write the integrand at one point.

#### What is the Galerkin Method? How Does it Connect to Code?

The Galerkin method chooses the approximation space equal to the test function space.
Concretely:

1. **Represent `u` as a linear combination of basis (shape) functions `φⱼ`:**
   ```
   u_h(x) = Σⱼ uⱼ · φⱼ(x)
   ```
   The `uⱼ` are the unknowns — one number per DOF (degree of freedom) per node.

2. **Choose test functions `ψᵢ = φᵢ` (same basis).** This is the "Galerkin" choice.

3. **The weak form then becomes a nonlinear algebraic system `R(u) = 0`:**
   ```
   Rᵢ = ∫_Ω ∇u_h · ∇φᵢ dΩ = Σⱼ uⱼ · ∫_Ω ∇φⱼ · ∇φᵢ dΩ
   ```

4. **Newton's method solves `R(u) = 0` iteratively:**
   ```
   J · δu = -R,    u ← u + δu
   ```
   where `J` is the Jacobian matrix `J_ij = ∂Rᵢ/∂uⱼ`.

**Connection to code:** The Jacobian entry `J_ij` at one quadrature point is what your
`computeQpJacobian()` returns. In `Diffusion.C`:

```cpp
Real Diffusion::computeQpJacobian() {
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];  // ∇φⱼ · ∇ψᵢ
}
```

`_phi[_j][_qp]` is the value of shape function `j` at quadrature point `_qp`. `_grad_phi`
is its gradient. These are provided by libMesh's finite element machinery. You access them
as pre-computed arrays in `Kernel.h`.

#### What is a Quadrature Point? How Does Numerical Integration Work?

On each element, integrals like `∫_element f(x) dΩ` cannot be evaluated exactly. Instead,
MOOSE uses Gaussian quadrature: pick a small set of points `x_qp` inside the element with
weights `w_qp`, then approximate:

```
∫_element f(x) dΩ  ≈  Σ_qp  f(x_qp) · w_qp · |J(x_qp)|
```

where `|J|` is the determinant of the mapping from reference to physical element.
In MOOSE, `_JxW[_qp]` already contains `w_qp · |J(x_qp)|`, so your integrand just needs
to return `f(x_qp)`.

For first-order Lagrange elements on quads (the default in MOOSE), Gaussian quadrature with
2×2 points is exact for polynomials up to degree 3 — which is more than sufficient for
first-order shape functions.

**Key insight:** Your `computeQpResidual()` is called once per quadrature point per test
function index `_i`. MOOSE handles the loops. You write the physics at one point.

---

### Installing MOOSE and Building the Test Application

Follow the official installation guide at `https://mooseframework.inl.gov/getting_started/`.
The recommended path is the conda package manager:

```bash
# After activating the moose conda environment:
cd /path/to/moose-next

# Build the framework library
cd framework && make -j8

# Build the test application (links framework + test objects)
cd ../test && make -j8
```

Key build variables (see `CLAUDE.md` in the repo root):
- `METHOD=opt` (default) — optimized build, no assertions
- `METHOD=dbg` — debug build with assertions and symbols; use this while learning
- `MOOSE_JOBS=N` — number of parallel make jobs

Build with `METHOD=dbg` while learning. Debug builds catch out-of-bounds access, null
pointer dereferences, and assertion failures that the optimized build silently ignores.

```bash
METHOD=dbg cd framework && make -j8
METHOD=dbg cd ../test && make -j8
```

---

### Success Check: Build and Run simple_diffusion

After a successful build, run the canonical first test:

```bash
cd test
./run_tests -j4 --re=simple_diffusion
```

You should see:
```
tests/kernels/simple_diffusion ... OK
```

Now run it manually to see what MOOSE prints:

```bash
./moose_test-dbg -i tests/kernels/simple_diffusion/simple_diffusion.i
```

Read the file `test/tests/kernels/simple_diffusion/simple_diffusion.i`. Every line should
now make sense given the theory above:

```
[Mesh]           GeneratedMesh: produces a 10x10 quad mesh on the unit square
[Variables]      u: the unknown field (first-order Lagrange, by default)
[Kernels]        Diffusion: contributes ∫ ∇u·∇ψᵢ dΩ  (the weak-form LHS)
[BCs]            DirichletBC: u=0 on left, u=1 on right (essential BCs)
[Executioner]    Steady, PJFNK: one Newton step (this problem is linear)
[Outputs]        exodus: write simple_diffusion_out.e, viewable in ParaView
```

**Checklist after Phase 0:**
- [ ] Build succeeds for both `framework` and `test`
- [ ] `./run_tests --re=simple_diffusion` passes
- [ ] Can explain in your own words: PDE, weak form, quadrature point, shape function, test function
- [ ] Know what `_grad_u[_qp]`, `_grad_test[_i][_qp]`, `_JxW[_qp]` represent
- [ ] Understand what `override` does and why MOOSE uses it everywhere

---

## Phase 1: User Foundations — Running Simulations

**Objective:** Become fluent in the MOOSE input file format. Run a variety of simulations.
Write new input files from scratch.

**Primary resources:**
- `docs/user-guide.md` — full input file syntax reference
- `docs/quick-start.md` — 63 worked examples
- `test/tests/kernels/` — real input files for every kernel type

---

### Day 1: The Input File Format (HIT Syntax)

Read the HIT (Hierarchical Input Text) section of `docs/user-guide.md`.

HIT is the input language for all MOOSE simulations. Every input file is a tree of
named blocks. Here is the structure:

```
[BlockName]
  parameter = value
  [SubBlockName]
    parameter = value
  []
[]
```

Key rules:
- Block names are case-sensitive
- Sub-objects are declared as sub-blocks: `[kernel_name]` inside `[Kernels]`
- `type = TypeName` is how MOOSE knows which C++ class to instantiate
- Strings containing spaces need single quotes: `petsc_options_value = 'hypre boomeramg'`
- Inline comments use `#`
- Variable substitution: `AD = ''` at the top, then `type = ${AD}Diffusion` expands to `type = Diffusion`

Read `test/tests/kernels/simple_diffusion/simple_diffusion.i` (45 lines) as your
reference. Every block maps to a C++ system.

**Exercise 1.1:** Copy `simple_diffusion.i` to a working directory. Change the mesh to 3D
(`dim = 3`, add `nz = 5`). Add `[right]` and `[left]` BCs on x faces, and run it. What
changes in the output?

---

### Day 2: Mesh and Variables

Read the `[Mesh]` and `[Variables]` sections in `docs/user-guide.md`.

**Mesh:** `GeneratedMesh` and `GeneratedMeshGenerator` produce structured grids:

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10   ny = 10
  xmin = 0  xmax = 2.0
[]
```

Boundary names for `GeneratedMesh`: `left`, `right`, `top`, `bottom` (2D); additionally
`front`, `back` (3D).

**Variables:** Each `[u]` sub-block declares one unknown field. Default is first-order
Lagrange. Override with:

```
[Variables]
  [T]
    order = SECOND
    family = LAGRANGE
    initial_condition = 300.0
  []
[]
```

Higher-order elements need more quadrature points — MOOSE selects these automatically.

**Exercise 1.2:** Run `test/tests/kernels/2d_diffusion/neumannbc.i`. Inspect the input.
Notice `type = NeumannBC` — this is the natural boundary condition that appears from the
integration-by-parts boundary term. Compare the solution to `simple_diffusion.i`.

**Exercise 1.3:** Run `test/tests/kernels/2d_diffusion/matdiffusion.i`. This uses a
material-dependent diffusivity `D = 0.01 + u²`. Notice `type = MatDiffusion` and a
`[Materials]` block with `DerivativeParsedMaterial`.

---

### Day 3: Kernels

Read the `[Kernels]` section in `docs/user-guide.md`.

Each kernel contributes one term to the weak-form residual for one variable. Multiple
kernels for the same variable are added together — this is how you build a complete PDE.

Standard kernels in `framework/include/kernels/`:

| Type | Weak form contribution | Use case |
|---|---|---|
| `Diffusion` | `∫ ∇u·∇ψᵢ dΩ` | `-∇²u = 0` |
| `TimeDerivative` | `∫ (∂u/∂t) ψᵢ dΩ` | transient term |
| `Reaction` | `∫ λ u ψᵢ dΩ` | `-λu` reaction |
| `BodyForce` | `∫ f ψᵢ dΩ` | source term `f` |
| `ADDiffusion` | same as `Diffusion`, but AD | automatic Jacobian |
| `MatDiffusion` | `∫ D(x) ∇u·∇ψᵢ dΩ` | non-uniform diffusivity |

Study `framework/src/kernels/Diffusion.C` (35 lines). Memorize this pattern — it is the
template for every kernel you will ever write:

```cpp
Real Diffusion::computeQpResidual() {
  return _grad_u[_qp] * _grad_test[_i][_qp];
}

Real Diffusion::computeQpJacobian() {
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
```

Study `framework/src/kernels/Reaction.C` (52 lines). Notice:
- `getParam<Real>("rate")` — how you read parameters from the input file
- `_rate` stored as a `const Real &` — the framework re-reads the parameter value
  when controls change it at runtime

**Exercise 1.4:** Run `test/tests/kernels/ad_transient_diffusion/ad_transient_diffusion.i`.
This solves `∂u/∂t - 0.1∇²u = 0` with `type = Transient`, `num_steps = 20`, `dt = 0.1`.
View the time evolution in ParaView. Add a `Reaction` kernel with `rate = 0.5` and observe
how the solution changes.

---

### Day 4: Boundary Conditions

Read the `[BCs]` section in `docs/user-guide.md`.

**Three BC types:**

1. **Dirichlet (essential):** Prescribe `u = value` on a boundary. Imposed by modifying
   the assembled system, not by a residual contribution.
   ```
   type = DirichletBC    boundary = left    value = 0
   ```

2. **Neumann (natural):** Prescribe the normal flux `∂u/∂n = value`. Arises naturally
   from the integration-by-parts surface term. Implemented as `IntegratedBC`:
   ```
   type = NeumannBC    boundary = right    value = 1.0
   ```
   The residual contribution is `-∫_∂Ω value · ψᵢ dS`.

3. **Robin:** A linear combination: `∂u/∂n + α·u = β`. In MOOSE, see `ADRobinBC`
   in `framework/include/bcs/ADRobinBC.h`.

**Key files to read:**
- `framework/include/bcs/IntegratedBC.h` — notice `_normals` (outward normal at QPs),
  `_phi`, `_test`, `_u`, `_grad_u` — same members as Kernel but on boundary faces
- `framework/include/bcs/NeumannBC.h` — a clean two-method implementation

**Exercise 1.5:** Take the simple_diffusion problem. Replace the right `DirichletBC`
with a `NeumannBC` with `value = 2.0`. How does the solution change? Why?

---

### Day 5: Materials, Executioners, Postprocessors

Read the `[Materials]`, `[Executioner]`, and `[Postprocessors]` sections.

**Materials** compute properties at quadrature points. The framework calls
`computeQpProperties()` on each material for each element before kernels run:

```
[Materials]
  [thermal]
    type = HeatConductionMaterial
    thermal_conductivity = 5.0
    specific_heat = 450.0
  []
[]
```

Properties declared in a material (e.g., `"thermal_conductivity"`) are retrieved by
kernels via `getMaterialProperty<Real>("thermal_conductivity")`.

**Executioners** control the solve strategy:
```
[Executioner]
  type = Steady         # or Transient
  solve_type = PJFNK    # Preconditioned JFNK (matrix-free Newton)
                        # or NEWTON (full Newton, needs SMP preconditioner)
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value  = 'hypre boomeramg'
[]
```

**Postprocessors** compute scalar quantities from the solution:
```
[Postprocessors]
  [average_T]
    type = ElementAverageValue
    variable = T
  []
  [max_flux]
    type = NodalMaxValue
    variable = T
  []
[]
```

Run `test/tests/postprocessors/element_average_value/element_average_value_test.i` to
see postprocessors with CSV output.

---

### Days 6-7: Worked Examples from quick-start.md

Work through `docs/quick-start.md` examples, correlating each to the matching test file.

**Exercise 1.6 (Day 6):** Write an input file from scratch for the 2D heat equation:
```
∂T/∂t - ∇·(k ∇T) = Q    in [0,1]²
T = 0   on left and bottom
∂T/∂n = 0  on right and top
T(x,y,0) = 0
Q = 1 (constant source)
```

Use:
- `type = GeneratedMesh`, `dim = 2`, `nx = 20`, `ny = 20`
- `type = HeatConduction` in `[Kernels]` with a `HeatConductionMaterial`
- `type = HeatConductionTimeDerivative` in `[Kernels]`
- `type = MatHeatSource` with `value = 1.0`
- `type = Transient`, `dt = 0.01`, `num_steps = 100`
- `type = ElementAverageValue` postprocessor for average temperature

Verify that the average temperature increases with time, approaching a steady state.

**Exercise 1.7 (Day 7):** Read `tutorials/darcy_thermo_mech/step02_darcy_pressure/problems/step2.i`
and `tutorials/darcy_thermo_mech/step08_postprocessors/problems/`. Trace how the physics
grows from step to step (Darcy pressure → add material → add heat conduction → add
postprocessors). This tutorial is the best existing walkthrough of MOOSE from a user
perspective.

---

### Success Check: Phase 1

**Gate:** Write a working input file for the steady 2D reaction-diffusion equation:
```
  -∇²u + 5u = f(x,y)     in [0,1]²
  u = 0                    on all boundaries
  f = 10·sin(π x)·sin(π y) (use ParsedFunction)
```

Assemble it entirely from memory using:
- `Diffusion` + `Reaction` (rate = 5) kernels
- `FunctionDirichletBC` with `function = 0`
- `BodyForce` with a `ParsedFunction` for `f`
- `ElementL2Error` postprocessor compared to the analytic solution
  `u* = 10/(π²+5) · sin(πx)·sin(πy)`

**Checklist after Phase 1:**
- [ ] Can write a working Steady and Transient input file from scratch
- [ ] Know what each top-level block (`Mesh`, `Variables`, `Kernels`, `BCs`, `Materials`,
      `Executioner`, `Postprocessors`, `Outputs`) does
- [ ] Understand the difference between Dirichlet and Neumann BCs
- [ ] Know when to use `PJFNK` vs `NEWTON` solve type
- [ ] Can add a postprocessor and have it write to a CSV file

---

## Phase 2: Understanding the Framework — How MOOSE Works

**Objective:** Understand how MOOSE transforms an input file into a running simulation.
Trace the call path from `type = Diffusion` in the input file to `computeQpResidual()`.

**Primary resources:**
- `docs/architecture.md` — complete subsystem diagrams
- `framework/include/base/MooseObject.h` — the root of every MOOSE class
- `framework/include/base/Factory.h` — how objects are created
- `framework/include/actions/Action.h` — how input blocks become objects
- `framework/include/loops/ComputeResidualThread.h` — the solve loop

---

### Day 8: The Three-Layer Stack

Read sections 1-3 of `docs/architecture.md`.

MOOSE sits on two scientific libraries:

```
+----------------------------------------------------------+
|  MOOSE (framework/include/, framework/src/)              |
|  Your objects: Kernels, BCs, Materials, PostProcessors   |
+---------------------------+---------------------------+--+
                            |                           |
                 libMesh (mesh, DOFs, FEM)    PETSc (nonlinear solver)
```

**libMesh** provides:
- `MeshBase` — the mesh (elements, nodes, connectivity)
- `DofMap` — maps unknowns to global equation numbers
- `FEBase` — evaluates shape functions and gradients at quadrature points
- `QuadratureRule` — the quadrature points and weights

**PETSc** provides:
- `SNES` — the nonlinear Newton solver (`solve_type = NEWTON` or `PJFNK`)
- `KSP` — the Krylov linear solver (GMRES, CG, etc.)
- `PC` — the preconditioner (ILU, AMG/hypre, etc.)

MOOSE's job is to fill PETSc's residual vector and Jacobian matrix from your physics objects.

**Exercise 2.1:** Read `framework/include/base/MooseObject.h`. Every single MOOSE plugin
class inherits from here. Notice `static InputParameters validParams()` — this static
method is how MOOSE discovers the parameters your object accepts before constructing it.

**Exercise 2.2:** Read `framework/include/base/Factory.h` (first 80 lines). The `Factory`
stores a registry of constructors keyed by string name. When the input file says
`type = Diffusion`, the Factory looks up the string `"Diffusion"`, calls the registered
constructor, and returns a `std::shared_ptr<MooseObject>`.

---

### Day 9: The Action System — From Input to Objects

Read section 6 of `docs/architecture.md`: "The Action System: From Input to Objects".

The path from input file to objects is:

```
Input file parsed by HIT
  → Builder (HIT walker) creates Actions via ActionFactory
    → ActionWarehouse orders Actions by task name
      → Actions call Factory::create<T>()
        → Factory calls registerMooseObject-registered constructors
          → Objects stored in FEProblemBase warehouses
```

**File to read:** `framework/src/actions/AddKernelAction.C` (52 lines).

This is the action that runs when MOOSE processes a sub-block of `[Kernels]`. The `act()`
method calls `_problem->addKernel(_type, _name, _moose_object_pars)`. The `_type` variable
is the value of `type =` from your input file — e.g., `"Diffusion"`.

`_problem->addKernel` looks in the Factory registry for `"Diffusion"`, calls
`Factory::create<KernelBase>("Diffusion", ...)`, which invokes `Diffusion::Diffusion(...)`.

**The `registerMooseObject` macro** (seen in every `.C` file):
```cpp
registerMooseObject("MooseApp", Diffusion);
```
This runs at program startup (via a static initializer) and registers a constructor for
`"Diffusion"` in the Factory. Without this line, `type = Diffusion` would fail with
"Unknown object type".

**Exercise 2.3:** Read `framework/include/actions/AddKernelAction.h` and
`framework/src/actions/AddKernelAction.C`. Then trace: when the parser encounters the
`[diff]` sub-block inside `[Kernels]` in `simple_diffusion.i`, what sequence of calls
creates the `Diffusion` object?

---

### Day 10: The Kernel Hierarchy

Read section 4 of `docs/architecture.md`: "The MooseObject Class Hierarchy".

Full inheritance chain for `Diffusion`:
```
MooseObject
  └── ResidualObject
        └── KernelBase         (framework/include/kernels/KernelBase.h)
              └── Kernel       (framework/include/kernels/Kernel.h)
                    └── Diffusion (framework/include/kernels/Diffusion.h)
```

Read each of these headers in order from bottom to top (`Diffusion.h` → `Kernel.h` →
`KernelBase.h`).

**From `KernelBase.h`:** Protected members available to every kernel:
```cpp
const Elem * const & _current_elem;  // current element being assembled
unsigned int _qp;                    // current quadrature point index
const MooseArray<Point> & _q_point;  // physical coordinates of quadrature points
const MooseArray<Real> & _JxW;       // quadrature weights * Jacobian det
unsigned int _i;                     // test function index
unsigned int _j;                     // shape function index (Jacobian only)
```

**From `Kernel.h`:** Additional members for scalar-variable kernels:
```cpp
MooseVariable & _var;                       // the variable this kernel acts on
const VariableTestValue & _test;            // test function values  ψᵢ[_qp]
const VariableTestGradient & _grad_test;    // test function gradients ∇ψᵢ[_qp]
const VariablePhiValue & _phi;              // shape function values  φⱼ[_qp]
const VariablePhiGradient & _grad_phi;      // shape function gradients ∇φⱼ[_qp]
const VariableValue & _u;                   // solution values  u_h[_qp]
const VariableGradient & _grad_u;           // solution gradient  ∇u_h[_qp]
```

**Exercise 2.4:** Read `framework/include/kernels/KernelBase.h` in full (77 lines).
Then read `framework/include/kernels/Kernel.h` in full (91 lines). For each data member,
write a one-sentence description of what it represents mathematically and when you would
use it in `computeQpResidual`.

---

### Day 11: The Solve Loop — One Newton Iteration

Read section 7 of `docs/architecture.md`: "The Solve Loop".

PETSc's SNES calls two functions: `residualFunction()` and `jacobianFunction()`. MOOSE
implements these as calls into `FEProblemBase::computeResidualInternal()` and
`FEProblemBase::computeJacobianInternal()`.

**Residual assembly for one Newton iteration:**

```
FEProblemBase::computeResidualInternal()
  └── NonlinearSystemBase::computeResidual()
        └── Threads::parallel_reduce(mesh elements, ComputeResidualThread)
              For each element:
                1. reinit(element)          — libMesh evaluates shape fns at QPs
                2. for each active Kernel:
                     for _i in 0..n_test_functions:
                       for _qp in 0..n_quadrature_points:
                         residual[i] += _JxW[_qp] * _coord[_qp] * computeQpResidual()
                3. Add local residual to global PETSc Vec
              For each boundary element (if BCs):
                4. same pattern for IntegratedBC::computeQpResidual()
```

**File to read:** `framework/include/loops/ComputeResidualThread.h` (44 lines).

**Jacobian assembly:**

Same structure but calls `computeQpJacobian()` for `_i`, `_j` pairs and assembles into
the global PETSc `Mat`. The Jacobian matrix entry `K[i][j]` is:
```
K[global_dof_i][global_dof_j] += _JxW[_qp] * _coord[_qp] * computeQpJacobian()
```

**Exercise 2.5:** Draw a call graph from `Steady::execute()` down to
`Diffusion::computeQpResidual()`. Start with `Steady.h` in `framework/include/executioners/`.
Use the architecture diagram in `docs/architecture.md` for guidance. Write the sequence of
method names with the class that owns each.

---

### Day 12: The Material System

Read section 8 of `docs/architecture.md`: "The Material System".

Materials are called **before** kernels on each element. This is critical: when a kernel
calls `getMaterialProperty<Real>("thermal_conductivity")`, it gets a reference to a
`MaterialProperty<Real>` object that was filled by the material's `computeQpProperties()`.

The material-kernel contract:
- Material declares `MaterialProperty<Real> & _k = declareProperty<Real>("k")` in its header
- Kernel retrieves `const MaterialProperty<Real> & _k = getMaterialProperty<Real>("k")` in its constructor
- On each element, at each QP: material's `computeQpProperties()` fills `_k[_qp]` first,
  then kernel's `computeQpResidual()` reads `_k[_qp]`

**File to read:** `framework/include/materials/Material.h` (319 lines, skim for structure).
Focus on `declareProperty<T>`, `getMaterialProperty<T>`, `computeQpProperties()`.

**Real module example:** Read `modules/heat_transfer/include/materials/HeatConductionMaterial.h`.
This is a template class (`HeatConductionMaterialTempl<bool is_ad>`) that works for both
the non-AD and AD cases. Notice `computeQpProperties()` is the key override.

---

### Day 13: Review and Code-Tracing Exercise

**Exercise 2.6 (Main exercise of Phase 2):** Starting from
`test/tests/kernels/simple_diffusion/simple_diffusion.i`, trace in writing:

1. How `type = Diffusion` is parsed (HIT → Builder → ActionFactory → AddKernelAction)
2. How `AddKernelAction::act()` creates the `Diffusion` object (Factory lookup)
3. How `Diffusion::computeQpResidual()` is called during a Newton iteration
   (Steady → FEProblemBase → NonlinearSystemBase → ComputeResidualThread → Diffusion)
4. What `_grad_u[_qp]` and `_grad_test[_i][_qp]` contain at the moment `computeQpResidual` is called
5. How the returned `Real` ends up in PETSc's residual vector

Write out each step. This exercise is the checkpoint for the whole phase.

---

### Success Check: Phase 2

**Gate:** Explain out loud (or in writing) what happens during one Newton iteration for
`simple_diffusion.i`. Your explanation must cover:
- Which C++ objects are involved and what method each calls
- What data `_grad_u[_qp]` holds and where it comes from
- Why the problem converges in one Newton step (the equation is linear)
- How the Jacobian entry relates to `computeQpJacobian()`

**Checklist after Phase 2:**
- [ ] Can explain the Factory pattern: how `type = Diffusion` maps to a C++ constructor call
- [ ] Know the inheritance chain: `Diffusion → Kernel → KernelBase → ResidualObject → MooseObject`
- [ ] Understand what `registerMooseObject("MooseApp", Diffusion)` does
- [ ] Can describe the solve loop: residual assembly → linear solve → update → repeat
- [ ] Know what `_JxW`, `_test`, `_phi`, `_u`, `_grad_u` are and when they are populated
- [ ] Understand the material-kernel evaluation order

---

## Phase 3: First C++ Object — Writing Your First Kernel

**Objective:** Write, compile, test, and validate a custom kernel. Understand the full
workflow from header file to passing test.

**Primary resources:**
- `docs/developer-guide.md` — sections 3-4 (Creating an App, Writing a Kernel)
- `stork/` — application template
- `tutorials/tutorial01_app_development/step05_kernel_object/` — worked example
- `tutorials/darcy_thermo_mech/step02_darcy_pressure/` — another worked example

---

### Day 14: Creating Your Application with stork

Read `docs/developer-guide.md` sections 1-3.

The fastest way to create a new application is the stork scaffolding script:

```bash
cd ~/projects
$MOOSE_DIR/scripts/stork.sh MyPhysicsApp
```

This creates:
```
MyPhysicsApp/
  Makefile
  src/
    main.C
    base/MyPhysicsAppApp.C
  include/
    base/MyPhysicsAppApp.h
  test/
    include/base/MyPhysicsAppTestApp.h
    src/base/MyPhysicsAppTestApp.C
    tests/
  unit/src/main.C
  unit/src/SampleTest.C
```

Study `tutorials/tutorial01_app_development/step05_kernel_object/src/base/BabblerApp.C`
to understand the `registerAll()` pattern:

```cpp
void BabblerApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax) {
  ModulesApp::registerAllObjects<BabblerApp>(f, af, syntax);  // registers all module objects
  Registry::registerObjectsTo(f, {"BabblerApp"});             // registers your app's objects
  Registry::registerActionsTo(af, {"BabblerApp"});
}
```

Build the scaffolded app to confirm it compiles:
```bash
cd MyPhysicsApp && make -j8 METHOD=dbg
```

---

### Day 15: Writing a Kernel — Non-AD First

Read `docs/developer-guide.md` section 4 (Writing Your First Kernel).

You will implement a **reaction-diffusion kernel**:
```
  -∇²u + κ(x) u = 0
```

where `κ(x)` is a spatially-varying reaction coefficient read from a material property.
The weak form:

```
  ∫_Ω ∇u·∇ψᵢ dΩ + ∫_Ω κ·u·ψᵢ dΩ = 0
```

You could use two separate kernels (`Diffusion` + `Reaction`), but implementing one
combined kernel teaches you more.

**Step 1:** Create the header `include/kernels/ReactionDiffusion.h`:

```cpp
#pragma once
#include "Kernel.h"

class ReactionDiffusion : public Kernel
{
public:
  static InputParameters validParams();
  ReactionDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Material property for the reaction coefficient
  const MaterialProperty<Real> & _kappa;
};
```

**Step 2:** Create the source `src/kernels/ReactionDiffusion.C`:

```cpp
#include "ReactionDiffusion.h"
registerMooseObject("MyPhysicsApp", ReactionDiffusion);

InputParameters ReactionDiffusion::validParams() {
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Implements -∇²u + κu = 0 as a single kernel.");
  params.addRequiredParam<MaterialPropertyName>("kappa",
      "Material property name for the reaction coefficient");
  return params;
}

ReactionDiffusion::ReactionDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _kappa(getMaterialProperty<Real>("kappa"))
{}

Real ReactionDiffusion::computeQpResidual() {
  // Diffusion term: ∇u · ∇ψᵢ
  // Reaction term:  κ · u · ψᵢ
  return _grad_u[_qp] * _grad_test[_i][_qp] + _kappa[_qp] * _u[_qp] * _test[_i][_qp];
}

Real ReactionDiffusion::computeQpJacobian() {
  // d/d(u_j) of residual = ∇φⱼ · ∇ψᵢ + κ · φⱼ · ψᵢ
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp]
       + _kappa[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}
```

**Key patterns to note:**
- `getMaterialProperty<Real>("kappa")` in the constructor — stores a `const` reference
  that is populated at runtime before `computeQpResidual` is called
- `_kappa[_qp]` in the residual — indexed by the current quadrature point
- The Jacobian is `∂R_i/∂u_j` — differentiate `computeQpResidual()` with respect to
  the unknown DOF `u_j = _phi[_j][_qp]`

---

### Day 16: Converting to ADKernel

Read `docs/developer-guide.md` section 9 (Automatic Differentiation).

Hand-coding `computeQpJacobian()` is error-prone for nonlinear physics. The AD system
uses dual numbers to compute exact Jacobians automatically.

**Differences when using `ADKernel`:**
- Inherit from `ADKernel` instead of `Kernel`
- Return type of `computeQpResidual()` is `ADReal` instead of `Real`
- The `_u`, `_grad_u`, `_kappa` values become `ADReal` / `ADRealVectorValue` —
  they carry derivative information automatically
- You do **not** implement `computeQpJacobian()` — AD derives it for you

**Rewrite as `ADReactionDiffusion.h`:**

```cpp
#pragma once
#include "ADKernel.h"

class ADReactionDiffusion : public ADKernel
{
public:
  static InputParameters validParams();
  ADReactionDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _kappa;  // note: ADMaterialProperty, not MaterialProperty
};
```

**Source `ADReactionDiffusion.C`:**

```cpp
#include "ADReactionDiffusion.h"
registerMooseObject("MyPhysicsApp", ADReactionDiffusion);

// validParams identical to non-AD version...

ADReactionDiffusion::ADReactionDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _kappa(getADMaterialProperty<Real>("kappa"))  // note: getADMaterialProperty
{}

ADReal ADReactionDiffusion::computeQpResidual() {
  return _grad_u[_qp] * _grad_test[_i][_qp] + _kappa[_qp] * _u[_qp] * _test[_i][_qp];
}
```

The residual body is identical to the non-AD version. The AD machinery propagates
derivatives through the dual-number arithmetic, producing the exact Jacobian at no
additional coding cost.

Look at `tutorials/darcy_thermo_mech/step02_darcy_pressure/src/kernels/DarcyPressure.C`
for a complete, working ADKernel example:

```cpp
ADReal DarcyPressure::computeQpResidual() {
  return (_permeability / _viscosity) * _grad_test[_i][_qp] * _grad_u[_qp];
}
```

---

### Day 17: Writing Tests Using the Test Harness

Read `docs/developer-guide.md` section 11 (Testing).

**Test harness basics:** Tests live in subdirectories of `test/tests/`. Each directory
has a `tests` file (HIT format) and one or more input files. Run tests with:

```bash
cd MyPhysicsApp && ./run_tests -j4
```

**A minimal test spec (`test/tests/kernels/reaction_diffusion/tests`):**

```
[Tests]
  [steady_test]
    type = Exodiff
    input = reaction_diffusion.i
    exodiff = reaction_diffusion_out.e
    issues = '#1'
    design = 'ReactionDiffusion.md'
    requirement = 'The system shall solve the reaction-diffusion equation
                   with a constant coefficient and match the gold file.'
  []
  [jacobian_test]
    type = PetscJacobianTester
    input = reaction_diffusion.i
    ratio_tol = 1e-7
    difference_tol = 1e-6
    requirement = 'The Jacobian computed by ReactionDiffusion shall match
                   the finite-difference Jacobian to within tolerance.'
  []
[]
```

**Jacobian testing with `PetscJacobianTester`:**

The most important test for any kernel is whether the Jacobian is correct. MOOSE provides
a built-in test type that compares your analytical `computeQpJacobian()` to a
finite-difference Jacobian computed by PETSc:

```
-snes_test_jacobian      # PETSc flag to check the Jacobian
-snes_test_jacobian_view # print the comparison
```

`PetscJacobianTester` automates this. A passing Jacobian test confirms:
- Newton will converge quadratically (not just linearly)
- No typos in your `computeQpJacobian()` implementation

Look at the real test spec in the framework:
`test/tests/kernels/simple_diffusion/tests` — note the required fields `issues`, `design`,
`requirement`. These are enforced by MOOSE's SQA (Software Quality Assurance) system.

**Generating gold files:**

```bash
cd test/tests/kernels/reaction_diffusion
# Run once to create the output
../../MyPhysicsApp-dbg -i reaction_diffusion.i
# Copy output as the gold standard
mkdir gold && cp reaction_diffusion_out.e gold/
```

On future runs, `Exodiff` compares output to gold within a tolerance.

---

### Success Check: Phase 3

**Gate:** Your `ReactionDiffusion` (or `ADReactionDiffusion`) kernel must pass the
`PetscJacobianTester`. This means:

1. Build your app with `METHOD=dbg`
2. Run `./run_tests --re=jacobian_test`
3. See `OK`

For the non-AD version: the Jacobian test failing means your `computeQpJacobian()` has
a typo. Compare each term against the mathematical derivative `∂R_i/∂u_j`.

**Checklist after Phase 3:**
- [ ] Created an application with stork
- [ ] Implemented `computeQpResidual()` and `computeQpJacobian()` for a custom kernel
- [ ] Know the difference between `Kernel` and `ADKernel` and when to use each
- [ ] `PetscJacobianTester` passes for your kernel
- [ ] Know what `registerMooseObject` does and why it is required
- [ ] Can write a `tests` spec with `Exodiff` and `PetscJacobianTester` types

---

## Phase 4: Materials and BCs — Building a Complete App

**Objective:** Implement a temperature-dependent material, custom boundary conditions, and
postprocessors. Assemble a complete mini-app with all components tested.

**Primary resources:**
- `docs/developer-guide.md` sections 5-8 (Material, BCs, Postprocessor, UserObject)
- `modules/heat_transfer/include/` — real module for reference
- `tutorials/darcy_thermo_mech/step03_darcy_material/` — material tutorial

---

### Day 18: Writing a Material

Read `docs/developer-guide.md` section 5 (Writing a Material).

Materials compute properties at quadrature points. The key override is
`computeQpProperties()`.

**Example: temperature-dependent thermal conductivity**

The Wiedemann-Franz law models `k(T) = k₀ · (1 + α(T - T_ref))`:

**Header `include/materials/LinearConductivity.h`:**

```cpp
#pragma once
#include "Material.h"

class LinearConductivity : public Material
{
public:
  static InputParameters validParams();
  LinearConductivity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Input parameters
  const Real & _k0;
  const Real & _alpha;
  const Real & _T_ref;

  // Coupled variable (temperature)
  const VariableValue & _T;

  // Declared property — other objects retrieve this
  MaterialProperty<Real> & _conductivity;
};
```

**Source `src/materials/LinearConductivity.C`:**

```cpp
#include "LinearConductivity.h"
registerMooseObject("MyPhysicsApp", LinearConductivity);

InputParameters LinearConductivity::validParams() {
  InputParameters params = Material::validParams();
  params.addClassDescription("Temperature-dependent conductivity k = k0*(1 + alpha*(T - T_ref))");
  params.addParam<Real>("k0", 1.0, "Reference conductivity at T_ref");
  params.addParam<Real>("alpha", 0.0, "Temperature coefficient");
  params.addParam<Real>("T_ref", 300.0, "Reference temperature");
  params.addCoupledVar("temperature", "The temperature variable");
  return params;
}

LinearConductivity::LinearConductivity(const InputParameters & parameters)
  : Material(parameters),
    _k0(getParam<Real>("k0")),
    _alpha(getParam<Real>("alpha")),
    _T_ref(getParam<Real>("T_ref")),
    _T(coupledValue("temperature")),
    _conductivity(declareProperty<Real>("thermal_conductivity"))  // declare the output
{}

void LinearConductivity::computeQpProperties() {
  _conductivity[_qp] = _k0 * (1.0 + _alpha * (_T[_qp] - _T_ref));
}
```

Study `modules/heat_transfer/include/materials/HeatConductionMaterial.h` and
`tutorials/darcy_thermo_mech/step03_darcy_material/include/materials/PackedColumn.h` for
real examples.

**Exercise 4.1:** Implement `LinearConductivity`. Write a test where `T = sin(πx)` (use
`FunctionIC`) and verify that the conductivity field has the expected spatial variation
by checking it with `ElementAverageMaterialProperty` postprocessor.

---

### Day 19: Writing Boundary Conditions

Read `docs/developer-guide.md` section 6 (Writing Boundary Conditions).

**Three BC patterns:**

**1. Dirichlet (NodalBC) — prescribes the value:**

```cpp
class MyDirichletBC : public NodalBC {
protected:
  Real computeQpResidual() override { return _u[_qp] - _value; }
};
```

**2. Neumann (IntegratedBC) — prescribes the flux:**

```cpp
class MyNeumannBC : public IntegratedBC {
protected:
  Real computeQpResidual() override { return -_flux * _test[_i][_qp]; }
};
```

The sign convention: the kernel residual uses integration-by-parts, so the natural BC
term is `- ∮ (∇u·n) ψᵢ dS`. Prescribing `∇u·n = flux`, the BC residual is `-flux * ψᵢ`.

**3. Robin (IntegratedBC) — a combination:**

A convective heat transfer BC: `k ∂T/∂n = h(T - T∞)` where `h` is a heat transfer
coefficient and `T∞` is the ambient temperature. Rearranged: `k ∂T/∂n - hT = -hT∞`.
The residual contribution:

```cpp
Real MyRobinBC::computeQpResidual() {
  // Recall: residual form is R = 0
  // The kernel already contributes ∫ k ∇T·∇ψᵢ on interior.
  // This BC adds ∫ h(T - T∞) ψᵢ dS (from the integration-by-parts boundary term).
  return _h * (_u[_qp] - _T_inf) * _test[_i][_qp];
}
Real MyRobinBC::computeQpJacobian() {
  return _h * _phi[_j][_qp] * _test[_i][_qp];
}
```

Look at `framework/include/bcs/ADRobinBC.h` and `framework/src/bcs/ADRobinBC.C` for the
real AD implementation. Look at `modules/heat_transfer/include/bcs/ConvectiveHeatFluxBC.h`
for a production convective BC.

**Exercise 4.2:** Implement a Robin (convective) BC class. Test it against the analytic
solution for 1D heat conduction with convection at one end: prescribe `T = T_hot` on
`x = 0` and `k dT/dx = h(T - T_inf)` on `x = L`. The steady solution is:
```
  T(x) = T_inf + (T_hot - T_inf) * (1 - x/(L + k/h))   when k/h << L
```

---

### Day 20: Writing Postprocessors and UserObjects

Read `docs/developer-guide.md` sections 7-8.

**Postprocessors** are scalar quantities computed from the solution. They inherit from
one of several base classes depending on when they execute:

```
ElementIntegralPostprocessor  — integrates over all elements
NodalPostprocessor            — evaluates at nodes
GeneralPostprocessor          — arbitrary scalar computation
SideIntegralPostprocessor     — integrates over boundaries
```

Read `framework/include/postprocessors/ElementIntegralPostprocessor.h`. The key override
is `computeQpIntegral()` — same quadrature-point pattern as kernels.

**Example: total flux through a boundary**

```cpp
class TotalBoundaryFlux : public SideIntegralVariablePostprocessor {
protected:
  Real computeQpIntegral() override {
    return -_diffusivity[_qp] * _grad_u[_qp] * _normals[_qp];
    //      D · ∇u · n   at boundary quadrature point
  }
  const MaterialProperty<Real> & _diffusivity;
};
```

**UserObjects** are a more general mechanism — they can communicate between objects
during a solve. `ElementUserObject` loops over elements; `GeneralUserObject` runs once.

**Exercise 4.3:** Write an `ElementIntegralPostprocessor` subclass that computes the
total energy in the domain: `∫_Ω ρ c_p T dΩ`. Use a coupled variable for `T` and
material properties for `ρ` and `c_p`. Verify using the analytic result for a uniform
temperature field.

---

### Days 21-22: Testing Deep Dive

Read `docs/developer-guide.md` section 11 (Testing) in full.

**All test types in MOOSE:**

| Type | Purpose |
|---|---|
| `Exodiff` | Compare Exodus output field by field, within tolerance |
| `CSVDiff` | Compare CSV postprocessor output (used for scalar verification) |
| `PetscJacobianTester` | Verify Jacobian against finite-difference approximation |
| `RunException` | Expect a specific error message |
| `PythonUnitTest` | Run Python unittest for tooling |
| `Unittest` | Run C++ Google Test unit tests |

**Convergence testing (the gold standard):**

For a manufactured solution, refine the mesh and confirm the error decreases at the
expected rate. For first-order Lagrange elements and an L2 error norm:
- Linear equation: order 2 convergence (error halves when mesh doubles)
- The ratio `error(h) / error(h/2)` should approach 4.0

Write a Python script or `CSVDiff` test that checks the convergence rate:

```
[Postprocessors]
  [l2_error]
    type = ElementL2Error
    variable = u
    function = exact_solution  # ParsedFunction with the MMS solution
  []
[]
```

Run on meshes `nx = 4, 8, 16, 32`. Compute `log(err_coarse / err_fine) / log(2)` for
each refinement. This is your convergence rate.

**Exercise 4.4 (Key exercise):** Verify your `ReactionDiffusion` kernel achieves second-
order spatial convergence using the Method of Manufactured Solutions (MMS):

1. Choose `u_exact = sin(πx)·sin(πy)` on `[0,1]²`
2. Compute the forcing function `f = -∇²u_exact + κ·u_exact` analytically
3. Add `BodyForce` with `f` as a `ParsedFunction`
4. Apply `DirichletBC` using `FunctionDirichletBC` with the exact solution
5. Compute `ElementL2Error` at meshes `nx = 4, 8, 16, 32, 64`
6. Verify the convergence rate is ≈ 2.0

---

### Success Check: Phase 4

**Gate:** Complete mini-app with all of the following tested and passing:
- Custom `ReactionDiffusion` kernel (or `ADReactionDiffusion`) — `PetscJacobianTester` passes
- Custom `LinearConductivity` material — integration test passes
- Custom Robin BC — `Exodiff` test against analytic solution passes
- Custom energy postprocessor — `CSVDiff` test against analytic value passes
- MMS convergence test — rate ≥ 1.9 for all refinement levels

**Checklist after Phase 4:**
- [ ] Know the three BC types (Dirichlet/NodalBC, Neumann/IntegratedBC, Robin/IntegratedBC) and when to use each
- [ ] Can implement a Material that couples to a variable and declares a property
- [ ] Can write an `ElementIntegralPostprocessor`
- [ ] Have run a convergence study and verified second-order spatial accuracy
- [ ] Know all five test types and when to use each
- [ ] Understand what goes in the `issues`, `design`, `requirement` fields of a test spec

---

## Phase 5: Advanced Topics — Mastering the Framework

**Objective:** Understand MultiApps, MeshGenerators, debugging techniques, and real module
code. After this phase you will be ready to design your own module.

---

### Day 23: MultiApp and Transfers

Read `docs/architecture.md` section 9 (MultiApp and Transfer System).

MultiApps enable multiscale and operator-split coupling. A parent application contains
sub-applications, each solving their own problem. Data flows between them via Transfers.

**File to read:** `test/tests/multiapps/full_solve_multiapp/parent.i`

The parent runs a steady diffusion, and also runs a transient sub-app:
```
[MultiApps]
  [full_solve]
    type = FullSolveMultiApp
    execute_on = initial
    positions = '0 0 0'
    input_files = sub.i
  []
[]
```

`execute_on = initial` means the sub-app runs once before the parent solve. Other options:
`timestep_begin`, `timestep_end`, `linear`, `nonlinear`.

Read `test/tests/multiapps/full_solve_multiapp/sub.i` — it is a standalone transient
problem that happens to run as a sub-app.

**Transfer types:**

| Type | Direction | What it moves |
|---|---|---|
| `MultiAppCopyTransfer` | parent ↔ child | copies a variable field |
| `MultiAppInterpolationTransfer` | parent → child | interpolates between dissimilar meshes |
| `MultiAppPostprocessorTransfer` | parent ↔ child | moves a scalar postprocessor value |

**Exercise 5.1:** Create a two-app coupling problem:
- Parent: steady thermal problem on a coarse mesh
- Sub-app: finer mesh with a reaction equation driven by the parent's temperature
- Transfer: parent's temperature `T` → sub-app's `T_ext` variable at each coupling step

---

### Day 24: MeshGenerators

Read `docs/developer-guide.md` section 15 (Writing a MeshGenerator).

`MeshGenerator` objects build the mesh before the solve. They form a pipeline: the
output of one generator can be the input of another.

```
[Mesh]
  [coarse]
    type = GeneratedMeshGenerator
    dim = 2   nx = 4   ny = 4
  []
  [refined]
    type = RefineBlockGenerator
    input = coarse
    block = '0'
    refinement = 2
  []
[]
```

Read `framework/include/meshgenerators/GeneratedMeshGenerator.h`. Notice:
- Inherits from `MeshGenerator`
- Overrides `std::unique_ptr<MeshBase> generate()`
- Stores parameters as member variables

**Exercise 5.2:** Write a `MeshGenerator` that takes a `GeneratedMesh` and adds a
subdomain assignment: elements with `x > 0.5` get `subdomain_id = 1`, others get
`subdomain_id = 0`. This lets you apply different materials to different regions.

---

### Day 25: Debugging Techniques

Read `docs/developer-guide.md` section 16 (Debugging Techniques).

**Debug build:**
```bash
METHOD=dbg make -j8
./MyApp-dbg -i problem.i
```
Debug builds enable assertions (`mooseAssert`), bounds checking, and symbol information
for debuggers. Always develop with `METHOD=dbg`.

**PETSc options for diagnostics:**

```
# Print linear iteration counts each step
petsc_options = '-ksp_monitor'

# Check the Jacobian (slow — only for debugging)
petsc_options = '-snes_test_jacobian'

# Print matrix structure
petsc_options = '-mat_view'

# Print residual norm each Newton step
petsc_options = '-snes_monitor'
```

**MOOSE diagnostic flags:**

```bash
# Print the parsed input file (useful to verify substitutions)
./MyApp-dbg -i problem.i --show-input

# List all registered object types
./MyApp-dbg --show-objects

# Print the action execution order
./MyApp-dbg -i problem.i --show-actions

# Show mesh information
./MyApp-dbg -i problem.i --show-mesh-info
```

**GDB / LLDB:** Build with `METHOD=dbg`. Set a breakpoint in your kernel:
```bash
gdb ./MyApp-dbg
(gdb) break MyKernel::computeQpResidual
(gdb) run -i problem.i
```

**mooseError and mooseAssert:**
```cpp
// In your kernel:
mooseAssert(_kappa[_qp] > 0, "Negative conductivity is unphysical: " << _kappa[_qp]);
if (_u[_qp] < 0)
  mooseWarning("Negative concentration detected at ", _q_point[_qp]);
```

---

### Days 26-27: Study the heat_transfer Module

The `modules/heat_transfer/` module is the best-documented, most mature physics module
in MOOSE. Study it as a template for your own module.

**Directory structure:**
```
modules/heat_transfer/
  include/
    kernels/        HeatConduction.h, ADHeatConduction.h, HeatSource.h, ...
    bcs/            ConvectiveHeatFluxBC.h, ADConvectiveHeatFluxBC.h, FunctionRadiativeBC.h, ...
    materials/      HeatConductionMaterial.h, ...
    postprocessors/ ...
    actions/        ...
    physics/        ...
  src/              (mirrors include/ with .C files)
  test/             (test inputs and gold files)
  doc/              (MooseDocs documentation)
```

**Read these files in order:**

1. `modules/heat_transfer/include/kernels/HeatConduction.h` — notice it inherits from
   `Diffusion` and adds a material property `_diffusion_coefficient`

2. `modules/heat_transfer/include/kernels/ADHeatConduction.h` — inherits from `ADDiffusion`,
   uses `ADMaterialProperty<Real> & _thermal_conductivity`

3. `modules/heat_transfer/src/kernels/HeatConduction.C` — only 50 lines; `computeQpResidual`
   multiplies the parent's `_grad_u * _grad_test` by the material property

4. `modules/heat_transfer/include/materials/HeatConductionMaterial.h` — a template class
   `HeatConductionMaterialTempl<is_ad>` that handles both AD and non-AD cases with a
   `typedef` pair at the bottom

5. `modules/heat_transfer/include/bcs/ConvectiveHeatFluxBC.h` — the Robin-type convective BC
   `h(T - T_inf)` as an `IntegratedBC`

**Exercise 5.3:** For each file you read, answer:
- What PDE term does this object represent mathematically?
- What base class does it extend, and why that choice?
- What input parameters does it expose (read `validParams()`)?
- How does it couple to other objects (material properties, coupled variables)?

---

### Day 28: Study solid_mechanics Module

`modules/solid_mechanics/` is larger and more complex than `heat_transfer`. Study it
to understand how a mature, multi-physics module is structured.

Key files to read:

1. `modules/solid_mechanics/include/kernels/ADStressDivergenceTensors.h` — the primary
   kernel for the strong form `∇·σ + b = 0`. Notice it uses tensor-valued material properties.

2. `modules/solid_mechanics/include/kernels/Gravity.h` — a simple body force kernel for
   gravitational loading; a clean three-method example.

The solid mechanics module shows how to use higher-rank tensor material properties
(`MaterialProperty<RankTwoTensor>`, `MaterialProperty<RankFourTensor>`) and how modules
can build on each other (solid mechanics depends on `misc` and `solid_properties`).

---

### Success Check: Phase 5

**Gate:** Given any kernel or material in `modules/heat_transfer/` or
`modules/solid_mechanics/`:
- Identify which base class it extends and why
- Explain what PDE term it represents
- Identify all material properties it consumes and produces
- Identify all input parameters it exposes
- Find the corresponding test and understand what it verifies

**Checklist after Phase 5:**
- [ ] Understand MultiApp coupling with `FullSolveMultiApp` and `TransientMultiApp`
- [ ] Know the Transfer types and when to use each
- [ ] Can write a `MeshGenerator` with a `generate()` override
- [ ] Know five PETSc flags for diagnosis and when to use them
- [ ] Have studied at least 10 files in `modules/heat_transfer/`
- [ ] Can explain the structure of a module (include/, src/, test/, doc/)
- [ ] Understand how `ADHeatConduction` and `HeatConductionMaterial` work together

---

## Phase 6: Building a Physics Module — Your Own Module

**Objective:** Design and implement a minimal but complete physics module: nonlinear
Poisson equation with radiation boundary condition. Write documentation. Write 5+
tests including convergence verification.

---

### Days 29-30: Module Design

**Target physics:** The nonlinear Poisson equation with radiation BC:

```
  -∇·(D(u) ∇u) = f             in Ω
  u = g                          on Γ_D    (Dirichlet)
  D(u) ∂u/∂n + σ ε (u⁴ - u_s⁴) = 0   on Γ_R    (radiation, Stefan-Boltzmann)
```

This is a useful real-world problem (radiation-conduction coupling). It is nonlinear in
two ways: the diffusivity `D(u)` depends on `u`, and the BC has a `u⁴` term.

**Module directory structure:**

```
NonlinearPoisson/
  Makefile             (copy from stork/Makefile.module, update LIBRARY_NAME)
  include/
    base/NonlinearPoissonApp.h
    kernels/NonlinearPoissonKernel.h
    materials/NonlinearConductivity.h
    bcs/RadiationBC.h
    postprocessors/TotalRadiationLoss.h
  src/
    base/NonlinearPoissonApp.C
    kernels/NonlinearPoissonKernel.C
    materials/NonlinearConductivity.C
    bcs/RadiationBC.C
    postprocessors/TotalRadiationLoss.C
  test/
    tests/
      kernels/
        nonlinear_poisson/
          tests
          nonlinear_poisson.i
          gold/
      bcs/
        radiation_bc/
          tests
          radiation_bc.i
          gold/
      convergence/
        mms_convergence.i
        tests
        gold/
  doc/
    content/
      modules/nonlinear_poisson/index.md
      modules/nonlinear_poisson/kernels/NonlinearPoissonKernel.md
```

**Day 29 exercise:** Write the header files for all five objects before writing any `.C`
files. For each header:
1. Choose the correct base class
2. Declare all data members
3. Write the `validParams()` signature and list all parameters in a comment

**Day 30 exercise:** Design the weak form on paper:

For the kernel: multiply `-∇·(D(u)∇u) = f` by `ψᵢ` and integrate by parts:
```
  R_i = ∫_Ω D(u) ∇u · ∇ψᵢ dΩ - ∫_∂Ω_R σε(u⁴ - u_s⁴) ψᵢ dS - ∫_Ω f ψᵢ dΩ
```

For the Jacobian: differentiate `R_i` with respect to `u_j = φⱼ`:
```
  J_ij = ∫_Ω [D'(u)φⱼ ∇u + D(u)∇φⱼ] · ∇ψᵢ dΩ
```

On the radiation BC face:
```
  J_ij^BC = ∫_∂Ω_R σε · 4u³ · φⱼ · ψᵢ dS
```

---

### Days 31-33: Implementation

**Day 31: Kernel and Material**

Implement `NonlinearPoissonKernel` as an `ADKernel` to avoid hand-coding the complex
Jacobian. The kernel uses a material property `D(u)`:

```cpp
ADReal NonlinearPoissonKernel::computeQpResidual() {
  return _D[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}
```

The material `NonlinearConductivity` declares `D` and computes it from `u`:

```cpp
void NonlinearConductivity::computeQpProperties() {
  // e.g., D = D0 * (1 + beta * u)  (linearized Arrhenius)
  _D[_qp] = _D0 * (1.0 + _beta * _T[_qp]);
}
```

Since the kernel is an `ADKernel`, the `_D[_qp]` must be an `ADMaterialProperty<Real>`
declared with `declareADProperty<Real>` in the material. The AD type propagates dual
numbers through the material into the kernel's residual, yielding an exact Jacobian.

**Day 32: Radiation BC**

The radiation BC `σ ε (u⁴ - u_s⁴) = 0` on `Γ_R`. Use `ADIntegratedBC`:

```cpp
class RadiationBC : public ADIntegratedBC {
  ADReal computeQpResidual() override {
    return _sigma_eps * (std::pow(_u[_qp], 4) - std::pow(_u_surroundings, 4))
           * _test[_i][_qp];
  }
  const Real _sigma_eps;  // σ * ε product
  const Real _u_surroundings;  // surroundings temperature T_s
};
```

Verify the Jacobian with `PetscJacobianTester`. The `u⁴` term means the Jacobian has
`4u³` — easy to get wrong by hand, trivially correct with AD.

**Day 33: Postprocessor and Documentation**

Implement `TotalRadiationLoss` as a `SideIntegralPostprocessor` that integrates
`σ ε (u⁴ - u_s⁴)` over the radiation boundary.

Write documentation stubs for each object. MOOSE uses MooseDocs markdown:

```markdown
# NonlinearPoissonKernel

!syntax description /Kernels/NonlinearPoissonKernel

## Description

Implements the diffusion term for the nonlinear Poisson equation...

## Example Input File Syntax

!listing test/tests/kernels/nonlinear_poisson/nonlinear_poisson.i block=Kernels

## Input Parameters

!syntax parameters /Kernels/NonlinearPoissonKernel
```

---

### Days 34-35: Testing with Convergence Verification

**Five required tests:**

**Test 1 — Basic solve (`Exodiff`):**
Solve `−∇·(D∇u) = 0` with constant `D`, Dirichlet BCs. Verify the output matches gold.

**Test 2 — Jacobian test (`PetscJacobianTester`):**
Run `NonlinearPoissonKernel` with nonlinear `D(u)`. Confirm that the AD-generated
Jacobian matches the finite-difference Jacobian within tolerance `1e-7`.

**Test 3 — Radiation BC test (`Exodiff`):**
1D rod, insulated ends except the right boundary has radiation BC. Verify steady-state
temperature and total radiation loss via `TotalRadiationLoss` postprocessor.

**Test 4 — MMS Convergence test (custom Python check):**
Choose `u_exact(x,y) = sin(πx)·sin(πy)` with `D(u) = 1 + u`.
Compute the manufactured source `f = -∇·(D(u_exact)∇u_exact)`.
Solve on meshes `nx = 4, 8, 16, 32, 64`. Compute `ElementL2Error`.
Verify convergence rate ≥ 1.9 (second order). Write a `CSVDiff` test with the
expected L2 errors pre-computed.

**Test 5 — Error handling (`RunException`):**
Pass a negative emissivity `epsilon < 0` to `RadiationBC`. Verify the expected
`mooseError` message appears.

**Template `tests` file:**

```
[Tests]
  [basic_solve]
    type = Exodiff
    input = nonlinear_poisson.i
    exodiff = nonlinear_poisson_out.e
    issues = '#1'
    design = 'modules/nonlinear_poisson/kernels/NonlinearPoissonKernel.md'
    requirement = 'The NonlinearPoissonKernel shall solve the nonlinear
                   Poisson equation with spatially-varying conductivity.'
  []
  [jacobian]
    type = PetscJacobianTester
    input = nonlinear_poisson.i
    ratio_tol = 1e-7
    issues = '#1'
    design = 'modules/nonlinear_poisson/kernels/NonlinearPoissonKernel.md'
    requirement = 'The Jacobian produced by NonlinearPoissonKernel shall
                   match the finite-difference Jacobian to within 1e-7.'
  []
  [mms_convergence]
    type = CSVDiff
    input = mms_convergence.i
    csvdiff = mms_convergence_out.csv
    issues = '#2'
    design = 'modules/nonlinear_poisson/kernels/NonlinearPoissonKernel.md'
    requirement = 'The NonlinearPoissonKernel shall achieve second-order
                   spatial convergence on a series of uniformly-refined meshes.'
  []
[]
```

---

### Success Check: Phase 6

**Gate:** All five tests pass:
```bash
cd NonlinearPoisson && ./run_tests -j4
```
Output:
```
tests/kernels/nonlinear_poisson/basic_solve ... OK
tests/kernels/nonlinear_poisson/jacobian ... OK
tests/bcs/radiation_bc/radiation_bc ... OK
tests/convergence/mms_convergence ... OK
tests/error_handling/negative_emissivity ... OK
```

**Checklist after Phase 6:**
- [ ] Module builds from scratch with `make -j8`
- [ ] All five tests pass (including the `PetscJacobianTester` for the nonlinear kernel)
- [ ] Convergence rate verified ≥ 1.9 on MMS test
- [ ] Each object has a brief MooseDocs documentation stub
- [ ] Module directory structure mirrors the `modules/heat_transfer/` layout

---

## Phase 7: Production Skills — Going Pro

**Objective:** Understand parallel execution, profiling, code review practices, and
the CI system. Prepare a contribution.

---

### Day 36: Large-Scale Parallel Runs and Profiling

**Running in parallel:**

```bash
# MPI parallel run — replace N with your core count
mpiexec -n 4 ./MyApp-opt -i problem.i

# Check parallel performance
./MyApp-opt -i problem.i --show-mesh-info
```

For scalability testing:
1. Run the same problem on 1, 2, 4, 8 cores
2. Measure wall time with `--timing` or the perf_graph output block
3. Ideal strong scaling: 4x cores → 4x speedup

**Performance profiling:**

Add `[Outputs] perf_graph = true []` to see a detailed breakdown of where time is spent:

```
Performance Graph:
  solve:                   3.2 s
    computeResidual:        1.8 s
      MaterialWarehouse:     0.4 s
      KernelLoop:            1.4 s
    computeJacobian:        1.2 s
    linearSolve:            0.2 s
```

If `KernelLoop` dominates: your kernels are expensive. If `MaterialWarehouse` dominates:
your materials do heavy computation. If `linearSolve` dominates: consider a better
preconditioner.

**PETSc solver tuning:**

```
# BoomerAMG (algebraic multigrid) — fast for elliptic problems
petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
petsc_options_value  = 'hypre boomeramg 0.7'

# Direct solver — reliable for small problems
petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
petsc_options_value  = 'lu mumps'

# Block Jacobi + ILU — good for larger parallel runs
petsc_options_iname = '-pc_type -sub_pc_type'
petsc_options_value  = 'bjacobi ilu'
```

---

### Day 37: Reading a Real MOOSE Pull Request

Go to `https://github.com/idaholab/moose/pulls`. Find a recently merged PR that adds a
new object (kernel, BC, material). Read:
- The PR description
- The changed files (`.h`, `.C`, test inputs, `tests` specs, documentation)
- The reviewer comments

For each changed file ask:
1. What was the design choice made and why?
2. What test covers this new feature?
3. Does the documentation match the implementation?

This exercise calibrates your understanding against real production code contributions.

---

### Day 38: Contributing Guidelines and CI

Read the contributing guide at `framework/doc/content/framework/contributing.md` (or the
online version at `https://mooseframework.inl.gov/framework/contributing.html`).

**Key requirements for any MOOSE contribution:**
- Every new object must have a test with `issues`, `design`, and `requirement` fields
- Tests must pass on the CI system (Civet, at `https://civet.inl.gov`)
- Code must pass clang-format (enforced by pre-commit hook):
  ```bash
  scripts/install-format-hook.sh  # run once in the repo root
  ```
- Python code must pass Ruff linter: `ruff check python/`
- New objects need MooseDocs documentation pages

**The CI system (Civet):**
- Every PR triggers an automated build-and-test pipeline
- Tests run on Linux, macOS, and HPC environments
- Results are reported on the PR as pass/fail
- Your contribution is not ready for merge until CI passes

**SQA (Software Quality Assurance) requirements:**

MOOSE is used in nuclear safety calculations, so it follows strict SQA practices
(based on ASME NQA-1). Every test must have:
- `issues` — the GitHub issue number that motivated this work
- `design` — the path to the documentation page that specifies this feature
- `requirement` — a sentence stating what the system shall do

---

### Success Check: Phase 7

**Gate (mock contribution):** Prepare a contribution to your `NonlinearPoisson` module
that adds one new feature (e.g., a `VolumetricRadiationSource` kernel for participating
media: `∫_Ω κ(u_r⁴ - u⁴) ψᵢ dΩ`). Your contribution must include:
1. Implementation (`.h` and `.C` files)
2. Test with all four required fields (`issues`, `design`, `requirement`, test type)
3. Documentation page
4. Passing `PetscJacobianTester`
5. Passing `clang-format` check

**Checklist after Phase 7:**
- [ ] Understand how to run in parallel with MPI and interpret parallel performance
- [ ] Know the five most-used PETSc solver option combinations
- [ ] Have read a real MOOSE PR and understood every changed file
- [ ] Know the three SQA required fields for every test spec
- [ ] Can run `clang-format` on your code and fix all style issues

---

## Appendix: Quick Reference Card

### Input File Syntax (HIT Format)

```
# Most common blocks and their most-used parameters

[Mesh]
  type = GeneratedMesh    # or GeneratedMeshGenerator (newer syntax)
  dim = 2                 # 1, 2, or 3
  nx = 10                 # elements in x
  xmin = 0  xmax = 1.0   # domain bounds
[]

[Variables]
  [u]
    order = FIRST         # FIRST (default), SECOND, THIRD
    family = LAGRANGE     # LAGRANGE (default), MONOMIAL (DG), HERMITE
    initial_condition = 0 # or use [ICs] block
  []
[]

[Kernels]
  [diff]
    type = Diffusion      # or ADDiffusion, MatDiffusion, ...
    variable = u          # which variable this kernel acts on
  []
[]

[BCs]
  [left]
    type = DirichletBC    # or FunctionDirichletBC, NeumannBC, ADRobinBC
    variable = u
    boundary = left       # GeneratedMesh names: left right top bottom front back
    value = 0
  []
[]

[Materials]
  [mat]
    type = GenericConstantMaterial  # simplest option for constant properties
    prop_names  = 'D  rho'
    prop_values = '1  1000'
  []
[]

[Executioner]
  type = Steady           # or Transient, Eigenvalue
  solve_type = PJFNK      # PJFNK (matrix-free, good for large), NEWTON (full, better conv.)
  petsc_options_iname = '-pc_type'
  petsc_options_value  = 'hypre'
  nl_rel_tol = 1e-10      # Newton relative tolerance
  nl_max_its = 50         # max Newton iterations
  l_tol = 1e-5            # linear solver tolerance
  # For Transient:
  # num_steps = 100  dt = 0.01  scheme = BDF2
[]

[Postprocessors]
  [avg]
    type = ElementAverageValue
    variable = u
  []
[]

[Outputs]
  exodus = true           # Exodus II (ParaView)
  csv = true              # CSV for postprocessors
  perf_graph = true       # performance profiling table
[]
```

---

### Most-Used Class Names and Their Purpose

| Class | Header | Purpose |
|---|---|---|
| `MooseObject` | `base/MooseObject.h` | Root of all MOOSE objects; provides `InputParameters` |
| `Kernel` | `kernels/Kernel.h` | Volume integral contribution for one scalar variable |
| `ADKernel` | `kernels/ADKernel.h` | Same as Kernel but with automatic Jacobian via AD |
| `KernelBase` | `kernels/KernelBase.h` | Shared base for all kernel types; provides `_qp`, `_i`, `_JxW`, etc. |
| `Diffusion` | `kernels/Diffusion.h` | `∫ ∇u·∇ψᵢ dΩ` — the Laplacian operator |
| `Reaction` | `kernels/Reaction.h` | `∫ λ u ψᵢ dΩ` — consuming reaction term |
| `TimeDerivative` | `kernels/TimeDerivative.h` | `∫ (du/dt) ψᵢ dΩ` — transient term |
| `BodyForce` | `kernels/BodyForce.h` | `∫ f ψᵢ dΩ` — source term |
| `IntegratedBC` | `bcs/IntegratedBC.h` | Boundary integral for Neumann/Robin BCs |
| `NodalBC` | `bcs/NodalBC.h` | Dirichlet BC (value prescribed at nodes) |
| `ADIntegratedBC` | `bcs/ADIntegratedBC.h` | IntegratedBC with automatic Jacobian |
| `ADRobinBC` | `bcs/ADRobinBC.h` | Coefficient×u Robin BC |
| `Material` | `materials/Material.h` | Computes material properties at QPs |
| `ADMaterial` | `materials/ADMaterial.h` | Material with AD support for use with ADKernels |
| `GenericConstantMaterial` | `materials/GenericConstantMaterial.h` | Constant properties from input |
| `ElementIntegralPostprocessor` | `postprocessors/ElementIntegralPostprocessor.h` | `∫_Ω f dΩ` scalar |
| `ElementAverageValue` | `postprocessors/ElementAverageValue.h` | Volume-averaged variable |
| `ElementL2Error` | `postprocessors/ElementL2Error.h` | L2 error against analytic function |
| `GeneralUserObject` | `userobjects/GeneralUserObject.h` | Arbitrary scalar computation |
| `Steady` | `executioners/Steady.h` | Single-solve executioner |
| `Transient` | `executioners/Transient.h` | Time-stepping executioner |
| `Factory` | `base/Factory.h` | Creates objects from registered type strings |
| `MooseApp` | `base/MooseApp.h` | Application root: owns Factory, ActionWarehouse |
| `FEProblemBase` | `problems/FEProblemBase.h` | Owns NonlinearSystem, Assembly, all warehouses |
| `AddKernelAction` | `actions/AddKernelAction.h` | Action that processes `[Kernels]` block |
| `GeneratedMeshGenerator` | `meshgenerators/GeneratedMeshGenerator.h` | Structured mesh generation |
| `MeshGenerator` | `meshgenerators/MeshGenerator.h` | Base for all mesh generators |
| `FullSolveMultiApp` | `multiapps/FullSolveMultiApp.h` | Sub-app that runs to steady state |
| `TransientMultiApp` | `multiapps/TransientMultiApp.h` | Sub-app that steps with parent |

---

### Common PETSc Solver Option Combinations

```bash
# 1. BoomerAMG (algebraic multigrid) — best for elliptic/parabolic, large meshes
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value  = 'hypre boomeramg'

# 2. Direct solver (MUMPS) — reliable, memory-intensive, good up to ~1M DOFs
petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
petsc_options_value  = 'lu mumps'

# 3. ILU — general-purpose, works on any problem, poor parallelism
petsc_options_iname = '-pc_type -pc_factor_levels'
petsc_options_value  = 'ilu 4'

# 4. Jacobian-free Newton-Krylov (JFNK) — saves memory, does not form Jacobian matrix
solve_type = PJFNK
petsc_options_iname = '-pc_type'
petsc_options_value  = 'hypre'  # still need a preconditioner for the linear system

# 5. Fieldplit (for multiphysics) — applies different PCs to different variables
petsc_options_iname = '-pc_type -pc_fieldsplit_type -fieldsplit_0_pc_type -fieldsplit_1_pc_type'
petsc_options_value  = 'fieldsplit additive hypre hypre'
```

---

### Debugging Flags Reference

| Flag / Option | Type | What it does |
|---|---|---|
| `METHOD=dbg` | Make variable | Enables assertions, bounds checks, debug symbols |
| `METHOD=devel` | Make variable | Optimized but with extra runtime checks |
| `--show-input` | CLI flag | Print the fully-expanded input file (after substitutions) |
| `--show-objects` | CLI flag | List all registered object type names |
| `--show-actions` | CLI flag | Print the action execution sequence |
| `--show-mesh-info` | CLI flag | Print mesh statistics (elements, nodes, partitioning) |
| `-snes_monitor` | PETSc option | Print residual norm at each Newton step |
| `-snes_test_jacobian` | PETSc option | Compare analytic vs finite-difference Jacobian |
| `-ksp_monitor` | PETSc option | Print linear iteration count |
| `-ksp_view` | PETSc option | Print the KSP solver configuration |
| `-mat_view` | PETSc option | Print the sparse matrix structure |
| `-log_view` | PETSc option | Print PETSc performance profile |
| `perf_graph = true` | Output param | Print MOOSE's built-in timing table |
| `mooseAssert(cond, msg)` | C++ macro | Abort with message if condition is false (dbg only) |
| `mooseError(msg)` | C++ macro | Always abort with message |
| `mooseWarning(msg)` | C++ macro | Print warning without aborting |

---

### Key Repo Paths Reference

```
framework/
  include/kernels/       Kernel.h, ADKernel.h, Diffusion.h, Reaction.h, TimeDerivative.h, ...
  include/bcs/           IntegratedBC.h, NodalBC.h, NeumannBC.h, ADRobinBC.h, ...
  include/materials/     Material.h, ADMaterial.h, GenericConstantMaterial.h, ...
  include/postprocessors/ ElementIntegralPostprocessor.h, ElementL2Error.h, ...
  include/userobjects/   GeneralUserObject.h, ElementUserObject.h, ...
  include/executioners/  Steady.h, Transient.h, ...
  include/base/          MooseObject.h, Factory.h, MooseApp.h, ...
  include/actions/       Action.h, AddKernelAction.h, ...
  include/loops/         ComputeResidualThread.h, ComputeJacobianThread.h, ...
  include/meshgenerators/ GeneratedMeshGenerator.h, MeshGenerator.h, ...
  src/kernels/           Diffusion.C, Reaction.C, TimeDerivative.C, ...

modules/heat_transfer/   Best reference module: kernels, materials, bcs, postprocessors
modules/solid_mechanics/ Tensor mechanics, contact, dynamics

stork/                   App scaffolding template
tutorials/
  darcy_thermo_mech/     Best user-facing tutorial: 11 steps from diffusion to multiphysics
  tutorial01_app_development/  Best developer tutorial: 10 steps building a complete app

test/tests/
  kernels/simple_diffusion/  simple_diffusion.i — the canonical first problem
  kernels/2d_diffusion/      neumannbc.i, matdiffusion.i — BC and material examples
  kernels/ad_simple_diffusion/  ad_simple_diffusion.i — AD version
  kernels/ad_transient_diffusion/  Transient + AD
  multiapps/full_solve_multiapp/   parent.i + sub.i — MultiApp example
  postprocessors/element_average_value/ — Postprocessor example

docs/
  architecture.md        Class hierarchy, Factory, Action system, solve loop
  developer-guide.md     Step-by-step developer tutorial with code examples
  user-guide.md          Input file syntax reference
  quick-start.md         63 worked simulation examples
```
