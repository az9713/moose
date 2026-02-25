# MOOSE Framework Developer Guide

**Audience:** Developers with C/C++/Java experience who are new to scientific computing.

This guide walks you through everything needed to write, test, and debug MOOSE-based applications.
No prior finite element or PDE knowledge is assumed, but you will pick up the necessary concepts
as you go.

---

## Table of Contents

1. [Prerequisites and Mental Model](#1-prerequisites-and-mental-model)
2. [Setting Up Development Environment](#2-setting-up-development-environment)
3. [Creating a New Application](#3-creating-a-new-application)
4. [Writing Your First Kernel](#4-writing-your-first-kernel)
5. [Writing a Material](#5-writing-a-material)
6. [Writing Boundary Conditions](#6-writing-boundary-conditions)
7. [Writing a Postprocessor](#7-writing-a-postprocessor)
8. [Writing a UserObject](#8-writing-a-userobject)
9. [Automatic Differentiation (AD)](#9-automatic-differentiation-ad)
10. [The Build System in Depth](#10-the-build-system-in-depth)
11. [Testing](#11-testing)
12. [Code Style and Conventions](#12-code-style-and-conventions)
13. [The validParams() Chain](#13-the-validparams-chain)
14. [Writing a MultiApp](#14-writing-a-multiapp)
15. [Writing a MeshGenerator](#15-writing-a-meshgenerator)
16. [Debugging Techniques](#16-debugging-techniques)

---

## 1. Prerequisites and Mental Model

### C++17 Concepts You Must Know

MOOSE relies heavily on modern C++. If any of the following is unfamiliar, spend an hour reviewing
it before reading on.

#### Inheritance and Virtual Functions

MOOSE is a framework of base classes. You write new physics by inheriting from a MOOSE base class
and overriding pure virtual methods.

```cpp
// MOOSE declares a pure virtual in its base class:
class Kernel {
public:
  virtual Real computeQpResidual() = 0;  // You MUST override this
  virtual Real computeQpJacobian() { return 0; }  // Optional override
};

// You provide the physics:
class MyKernel : public Kernel {
public:
  Real computeQpResidual() override;   // Required - you implement this
};
```

The `= 0` suffix makes a method pure virtual: the base class declares the interface, your subclass
supplies the body. MOOSE calls your override during the solve; you never call it yourself directly.

The `override` keyword on your subclass is mandatory in MOOSE style. It lets the compiler catch
typos in method signatures - if the signature does not match the base class, compilation fails.

#### Templates

MOOSE uses templates extensively for type-generic code. The most common pattern you encounter is
templated material properties:

```cpp
// Declare a property that holds a Real (double)
MaterialProperty<Real> & _conductivity = declareProperty<Real>("thermal_conductivity");

// Declare a property that holds a 3x3 tensor
MaterialProperty<RealTensorValue> & _stress = declareProperty<RealTensorValue>("stress");
```

You do not need to write your own templates to use MOOSE, but you need to be able to read template
syntax to understand error messages and API documentation.

#### Smart Pointers

MeshGenerators return `std::unique_ptr<MeshBase>`, which is a pointer that automatically deletes
its target when it goes out of scope. The rule in MOOSE is simple: when a method signature says
`unique_ptr`, you must move-assign it, not copy it:

```cpp
std::unique_ptr<MeshBase> generate() override
{
  auto mesh = buildMeshBaseObject();
  // ... populate mesh ...
  return mesh;  // moves ownership to the caller
}
```

Reference-counted `std::shared_ptr` appears in MultiApp ownership. You rarely need to create these
yourself; MOOSE creates them internally.

#### const References

Almost everything MOOSE passes to your code arrives as a `const` reference. This is a performance
choice: no copies are made. The `&` means "reference, not copy"; `const` means "you cannot
modify it". When you see `const VariableValue & _u`, that means `_u` is an alias to the real data
sitting inside MOOSE's assembly system, and you may read it but not write to it.

---

### Brief PDE and Weak Form Primer

You are solving a Partial Differential Equation (PDE). For example, steady heat conduction in a
domain (call it Omega) with no source:

```
-div( D * grad(u) ) = 0     in Omega
```

Here `u` is temperature, `D` is thermal conductivity, and `div`/`grad` are the divergence and
gradient operators from calculus. You want MOOSE to find the function `u(x, y, z)` that satisfies
this equation everywhere inside the domain plus the boundary conditions you specify.

Computers cannot work with continuous functions directly. The Galerkin finite element method
converts the PDE into a system of algebraic equations that a computer can solve. Here is the idea
in plain English:

1. **Discretise the domain.** Cut the geometry into small elements (triangles, quads, tetrahedra,
   hexahedra). The nodes at element corners and edges are where we will track unknown values.

2. **Represent u approximately.** Express the unknown function as a sum of simple "hat" functions
   (called shape functions or basis functions, denoted `phi_j`). Each hat function is 1.0 at node
   j and 0.0 everywhere else, so the sum `u_h = sum_j u_j * phi_j` is just an interpolation of
   the nodal values `u_j`.

3. **Multiply by a test function and integrate.** Take the original PDE, multiply both sides by
   another "hat" function `psi_i` (the test function, also called `phi_i`), and integrate over the
   domain. The result is called the weak form. For the diffusion equation after integration by
   parts:

   ```
   integral( D * grad(u_h) . grad(psi_i) ) dOmega = 0    for each i
   ```

   This gives one equation per node. Stack them into a system `K * u = f` where `K` is the
   stiffness matrix. Solving that linear system gives you the nodal values of `u`.

4. **Use quadrature to evaluate the integrals.** On each element MOOSE picks a small set of
   "quadrature points" with associated weights. The integral over an element is approximated as:

   ```
   sum over quadrature points: integrand(x_qp) * weight(x_qp) * det(Jacobian)
   ```

   The determinant of the mapping Jacobian is stored in `_JxW[_qp]`. The product of this with the
   coordinate-system scaling factor `_coord[_qp]` gives the volumetric weight for each quadrature
   point.

**Key insight:** Your `computeQpResidual()` method is the integrand. MOOSE handles the summation
loop over quadrature points and the summation over elements. You only need to return the value at
one quadrature point.

---

### The Central Mental Model

> **MOOSE calls your `computeQpResidual()` at each quadrature point inside each element. You return
> a single number. MOOSE multiplies it by the test function `_test[_i][_qp]` and accumulates it
> into the residual vector.**

Everything else is bookkeeping around this loop. Once this model is clear, reading any MOOSE kernel
becomes straightforward.

For the Jacobian (the derivative of the residual with respect to the solution), MOOSE calls your
`computeQpJacobian()` at each quadrature point. You return `d(residual) / d(u)` evaluated at the
current quadrature point. MOOSE multiplies by `_phi[_j][_qp] * _test[_i][_qp]` and accumulates
into the stiffness matrix.

---

## 2. Setting Up Development Environment

### Prerequisites

Before building MOOSE you need:

- A working C++17 compiler (GCC >= 9 or Clang >= 9)
- MPI (OpenMPI or MPICH)
- Python >= 3.6
- CMake >= 3.11 (for some dependencies)

MOOSE ships its own copies of PETSc, libMesh, and other dependencies as git submodules. Do not
install system-wide PETSc expecting MOOSE to use it; the submodule versions are built with options
specifically chosen for MOOSE.

### Building MOOSE

The canonical build sequence is:

```bash
# Step 1: Obtain the source
cd $HOME/projects
git clone https://github.com/idaholab/moose.git
cd moose

# Step 2: Build PETSc (this takes 20-40 minutes)
cd scripts
./update_and_rebuild_petsc.sh

# Step 3: Build libMesh (this takes 10-20 minutes)
./update_and_rebuild_libmesh.sh

# Step 4: Set the environment variable so MOOSE knows where it lives
export MOOSE_DIR=$HOME/projects/moose

# Step 5: Build the MOOSE framework
cd $MOOSE_DIR/framework
make -j8 METHOD=opt
```

The `-j8` flag runs 8 compile jobs in parallel. Adjust to the number of CPU cores on your machine.

If anything fails in steps 2-3, do not skip to step 4. MOOSE's build system is designed to use
the exact library versions it built; mismatched versions cause cryptic link errors.

### Build Types (METHOD)

MOOSE supports four build flavors controlled by the `METHOD` variable:

| METHOD | Compiler flags | Use case |
|--------|---------------|----------|
| `opt`  | `-O2 -DNDEBUG` | Production runs and performance testing. Assertions disabled. Fast. |
| `dbg`  | `-O0 -g`       | Day-to-day development. Full assertions, all debug checks active. Slow. |
| `devel`| `-O2 -g`       | Performance investigation with some debug info. Assertions still on. |
| `oprof`| `-O2 -pg`      | Profiling with gprof. Rarely needed today; use a sampling profiler instead. |

Start new development with `METHOD=dbg`. When a run is too slow to reproduce your test problem
at reasonable speed, switch to `METHOD=opt`. Never run final results with `METHOD=dbg` because the
performance difference can be 10x to 100x.

The binary name reflects the method: `moose-opt`, `moose-dbg`, `moose-devel`.

### Verifying the Build

The framework ships with a test application called `test`. Build and run it:

```bash
cd $MOOSE_DIR/test
make -j8 METHOD=opt

# Run a specific test to verify a working build
./run_tests -j8 --re=simple_diffusion
```

A passing run prints something like:
```
test:kernels/simple_diffusion/test .................................... OK
----------------------------------------------
Ran 1 test in 1.2 seconds
1 passed, 0 skipped, 0 failed
```

If this test fails, your build is broken. Do not proceed until it passes.

### Editor Setup

#### compile_commands.json

Modern editors (VS Code, CLion, Neovim with clangd) use `compile_commands.json` to understand
which compiler flags apply to each file. This enables accurate jump-to-definition, auto-complete,
and error highlighting.

Generate it from your application directory:

```bash
cd /path/to/your/app
make compile_commands.json
```

This command must be the only make target when invoked; mixing it with other targets is an error
(the build system enforces this). The file is created in the application root. Point your editor's
language server to it.

Note: `compile_commands.json` is built with `METHOD=dbg` by default (the build system overrides
your setting when generating it). This ensures all debug-mode flags appear in the database.

#### clang-format

MOOSE uses clang-format to enforce consistent code style. The configuration file is at
`$MOOSE_DIR/.clang-format` (and stork copies it to new apps). To format a file:

```bash
clang-format -i src/kernels/MyKernel.C
clang-format -i include/kernels/MyKernel.h
```

Configure your editor to run clang-format on save. Most editors have a plugin for this. The
stork-generated apps also include a git pre-commit hook installer:

```bash
cd /path/to/your/app
./scripts/install-format-hook.sh
```

After installation, every `git commit` will automatically format staged files.

### Environment Variables

| Variable | Meaning | Example |
|----------|---------|---------|
| `MOOSE_DIR` | Root of the MOOSE source tree | `export MOOSE_DIR=$HOME/projects/moose` |
| `LIBMESH_DIR` | libMesh installation directory | Usually `$MOOSE_DIR/libmesh/installed`; rarely set manually |
| `MOOSE_JOBS` | Default parallel job count for `run_tests` | `export MOOSE_JOBS=8` |

You almost never need to set `LIBMESH_DIR` manually; the Makefile infers it from `MOOSE_DIR`. Set
`MOOSE_DIR` in your shell's startup file (`.bashrc` or `.zshrc`) so it is available in all
sessions.

---

## 3. Creating a New Application

### Using stork.sh

MOOSE provides a scaffolding script called `stork.sh` that generates a complete application
skeleton from a template:

```bash
# Navigate to the directory that will CONTAIN your app (not inside MOOSE)
cd $HOME/projects

# Run stork - the name MUST be in CamelCase
$MOOSE_DIR/scripts/stork.sh MyApp
```

Stork validates the name (no special characters, cannot start with a digit), copies the template
directory tree, renames all files and all occurrences of "Stork" inside those files to your chosen
name, initialises a git repository, and reports the location.

After creation:

```bash
cd $HOME/projects/my_app    # stork converts CamelCase to snake_case for the directory
make -j8 METHOD=opt
./run_tests
```

If both commands succeed, your empty application builds and the framework passes its sanity check
within your app.

### Files Created by stork.sh

```
my_app/
  Makefile                  - Top-level build entry point
  testroot                  - TestHarness configuration
  run_tests                 - Convenience script to invoke the TestHarness
  .clang-format             - Code style configuration (copied from MOOSE)
  .gitignore                - Standard MOOSE ignores
  include/
    base/
      MyAppApp.h            - Application class declaration
  src/
    base/
      MyAppApp.C            - Application class definition; registerAll() lives here
    main.C                  - Binary entry point (rarely modified)
  test/
    tests/                  - Put your test input files here
    include/base/
      MyAppTestApp.h        - Test app class declaration
    src/base/
      MyAppTestApp.C        - Test app class definition
  doc/                      - Documentation configuration
  scripts/                  - Helper scripts including install-format-hook.sh
  unit/                     - Optional unit test directory
```

#### MyAppApp.h

```cpp
#pragma once

#include "MooseApp.h"

class MyAppApp : public MooseApp
{
public:
  static InputParameters validParams();

  MyAppApp(const InputParameters & parameters);
  virtual ~MyAppApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};
```

This declares the application class. It inherits from `MooseApp`, the root of all MOOSE
applications. The three static methods handle registration.

#### MyAppApp.C

```cpp
#include "MyAppApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
MyAppApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

MyAppApp::MyAppApp(const InputParameters & parameters) : MooseApp(parameters)
{
  MyAppApp::registerAll(_factory, _action_factory, _syntax);
}

MyAppApp::~MyAppApp() {}

void
MyAppApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<MyAppApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"MyAppApp"});
  Registry::registerActionsTo(af, {"MyAppApp"});
}

void
MyAppApp::registerApps()
{
  registerApp(MyAppApp);
}
```

The `registerAll` method is where MOOSE learns about all the objects your application provides.
When you write a new Kernel, Material, BC, or any other object, you add `registerMooseObject` to
its `.C` file (described in section 12). The `Registry::registerObjectsTo` call harvests all those
registrations.

You should not normally modify `registerAll`. Just add `registerMooseObject("MyAppApp", MyClass)`
in each new class's `.C` file and the system picks it up automatically.

#### Makefile Structure

```makefile
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk    # Sets up compiler, flags, dependency scanning
include $(FRAMEWORK_DIR)/moose.mk    # Sets moose_LIB, libmesh_LIB, etc.

# ... module toggles ...

APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := my_app
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk   # Discovers sources, builds library and binary
```

The important line is the final `include $(FRAMEWORK_DIR)/app.mk`. This file (located at
`framework/app.mk`) contains the logic that:

1. Uses `find` to recursively discover all `.C` files under `src/`
2. Computes corresponding `.o` object file names
3. Links everything into `libmy_app-$(METHOD).la` and then `my_app-$(METHOD)`

You do not need to list individual source files. Drop a `.C` file into any subdirectory of `src/`
and make will find and compile it on the next build.

#### testroot

The `testroot` file is read by the TestHarness when you invoke `./run_tests`. It configures
globally applied test settings:

```
app_name = my_app
allow_warnings = false
allow_unused = false
allow_override = false
```

- `app_name` must match the executable name prefix (the `-opt` suffix is appended automatically).
- `allow_warnings = false` causes a test to fail if the executable prints any MOOSE warning. This
  is the recommended setting during development to catch problems early.
- `allow_unused = false` causes a test to fail if the input file contains parameters that the
  object does not recognise. This catches typos in input files.
- `allow_override = false` causes a test to fail if the same parameter is set twice in the input.

---

## 4. Writing Your First Kernel

A Kernel in MOOSE represents a volumetric term in the PDE's weak form. Each Kernel is associated
with exactly one variable and contributes to the residual and Jacobian for that variable's degrees
of freedom.

We will implement the reaction-diffusion equation:

```
-div( D * grad(u) ) + lambda * u = 0
```

where `D` is a diffusion coefficient and `lambda` is a reaction rate. Both are simple scalar
constants specified in the input file.

### Step 1: Create the Header File

Create `include/kernels/ReactionDiffusion.h`:

```cpp
//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Line 1: Include the base class.
// Kernel is declared in framework/include/kernels/Kernel.h.
// It provides access to _u, _grad_u, _test, _grad_test, _phi, _grad_phi,
// _q_point, _JxW, _qp, _i, and _j.
#include "Kernel.h"

// Line 2: Declare the class.
// The class name MUST match the filename (without the .h extension).
// This is enforced by convention and by registerMooseObject.
class ReactionDiffusion : public Kernel
{
public:
  // Line 3: Every MOOSE object must have this static method.
  // It returns an InputParameters object that describes what the input
  // file block for this object may contain.
  static InputParameters validParams();

  // Line 4: Standard constructor signature for all MOOSE objects.
  // You receive a const reference to the parsed InputParameters.
  ReactionDiffusion(const InputParameters & parameters);

protected:
  // Line 5: The only method you are required to override.
  // MOOSE calls this for every (element, quadrature point, test function index)
  // combination. You return a Real (double precision floating point number).
  // The 'override' keyword tells the compiler to verify the signature matches
  // the base class; compilation fails if it does not.
  virtual Real computeQpResidual() override;

  // Line 6: Optional but strongly recommended.
  // You return d(residual)/d(u) at the current quadrature point.
  // If you skip this, MOOSE falls back to finite-difference Jacobian
  // estimation, which is 10-100x slower.
  virtual Real computeQpJacobian() override;

  // Line 7: Parameters read from the input file.
  // These are set once in the constructor and never change during the solve.
  // Declare them const to prevent accidental modification.
  const Real _D;        // diffusion coefficient
  const Real _lambda;   // reaction rate
};
```

### Step 2: Create the Source File

Create `src/kernels/ReactionDiffusion.C`:

```cpp
//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ReactionDiffusion.h"

// This macro registers your class with the MOOSE object factory.
// First argument: the application name string (must match what you pass to
//                 Registry::registerObjectsTo in your App's registerAll()).
// Second argument: the class name.
// After this registration the input file type = ReactionDiffusion works.
registerMooseObject("MyAppApp", ReactionDiffusion);

// validParams() describes the input file syntax for this object.
// It is a static method, so it runs at startup when MOOSE scans available
// objects. It does NOT run once per simulation step.
InputParameters
ReactionDiffusion::validParams()
{
  // CRITICAL: Always call the parent class validParams() first.
  // This adds parameters that Kernel already handles (variable, block, etc.)
  // without which the object cannot function.
  InputParameters params = Kernel::validParams();

  // addClassDescription sets the help text shown in --dump output and docs.
  params.addClassDescription(
      "Implements -div(D*grad(u)) + lambda*u = 0. "
      "Diffusion coefficient D and reaction rate lambda are input parameters.");

  // addRequiredParam means the user MUST provide this in the input file.
  // If omitted, MOOSE errors before the solve begins.
  params.addRequiredParam<Real>("D", "Diffusion coefficient (positive real number)");

  // addParam provides a default value used when the parameter is absent.
  // The third argument is the default, the fourth is documentation.
  params.addParam<Real>("lambda", 0.0, "Reaction rate coefficient");

  return params;
}

// Constructor. Always pass parameters to the base class.
// Read your parameters using getParam<Type>("name").
ReactionDiffusion::ReactionDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _D(getParam<Real>("D")),
    _lambda(getParam<Real>("lambda"))
{
  // Validate inputs at construction time. Errors here abort before the solve.
  if (_D <= 0.0)
    mooseError("ReactionDiffusion: D must be positive, got ", _D);
}

// computeQpResidual() is the heart of every Kernel.
//
// Context when this method is called:
//   _qp   : current quadrature point index (0-based integer)
//   _i    : current test function index (0-based integer, loops over all
//           test functions associated with the current element)
//
// Available member variables (all indexed by [_qp] or [_i][_qp]):
//   _u[_qp]           : solution value at this quadrature point
//   _grad_u[_qp]      : gradient of solution at this quadrature point
//                       (a RealGradient, i.e., a 3-component vector)
//   _test[_i][_qp]    : value of test function i at this quadrature point
//   _grad_test[_i][_qp]: gradient of test function i at this quadrature point
//   _q_point[_qp]     : physical (x,y,z) coordinates of this quadrature point
//   _JxW[_qp]         : quadrature weight * mapping Jacobian determinant
//                       (the "volume" attributed to this quadrature point)
//   _coord[_qp]       : coordinate-system scaling factor (1 for Cartesian,
//                       r for axisymmetric, r^2*sin(theta) for spherical)
//
// IMPORTANT: You do NOT multiply by _test[_i][_qp] yourself here.
// MOOSE does that multiplication for you AFTER your method returns.
// You return the integrand WITHOUT the test function factor.
//
// Wait - but the weak form HAS a test function factor. Why not include it?
// Because the loop over test function index _i is managed by MOOSE, and
// it multiplies your return value by _test[_i][_qp] internally. This design
// lets you focus purely on the physics expression.
//
// For the diffusion term -div(D*grad(u)), after integration by parts the
// weak form is: D * grad(u) . grad(test_i)
// The '.' is the dot product of two RealGradient vectors.
//
// For the reaction term +lambda*u, the weak form is: lambda * u * test_i
// But since MOOSE multiplies by test_i externally, you return lambda * u.
Real
ReactionDiffusion::computeQpResidual()
{
  // Diffusion term: D * grad(u)[qp] . grad(test[i])[qp]
  // _grad_u[_qp] * _grad_test[_i][_qp] is the dot product of two RealGradient
  // objects. The * operator between RealGradient objects is overloaded to be
  // the dot product.
  //
  // Reaction term: lambda * u[qp]
  // MOOSE multiplies the returned value by _test[_i][_qp] and _JxW[_qp]
  // and _coord[_qp] after this method returns.
  return _D * _grad_u[_qp] * _grad_test[_i][_qp] + _lambda * _u[_qp] * _test[_i][_qp];
}

// computeQpJacobian() returns d(computeQpResidual()) / d(u)
//
// Additional context:
//   _j   : current shape function index (0-based integer)
//   _phi[_j][_qp]      : value of shape function j at this quadrature point
//   _grad_phi[_j][_qp] : gradient of shape function j at this quadrature point
//
// The Jacobian entry K[i][j] = integral( d(R_i)/d(u_j) )
// In the FEM context, u_h = sum_j u_j * phi_j, so d(u_h)/d(u_j) = phi_j
// and d(grad(u_h))/d(u_j) = grad(phi_j).
//
// Differentiate computeQpResidual() with respect to u_j:
//   d/du_j ( D * grad(u) . grad(test_i) + lambda * u * test_i )
// = D * grad(phi_j) . grad(test_i) + lambda * phi_j * test_i
//
// Again, MOOSE multiplies your return value by _test[_i][_qp] * _JxW * _coord
// after this method returns, so you return only the pre-test-function part.
// For the diffusion term that means returning D * grad_phi_j . grad_test_i.
// For the reaction term that means returning lambda * phi_j * test_i -- but
// since MOOSE multiplies by _test[_i][_qp] externally, you only return
// lambda * phi_j.
Real
ReactionDiffusion::computeQpJacobian()
{
  return _D * _grad_phi[_j][_qp] * _grad_test[_i][_qp] + _lambda * _phi[_j][_qp] * _test[_i][_qp];
}
```

### Explaining Every Protected Member Variable

These member variables are inherited from `Kernel` and `KernelBase`. They are available in
`computeQpResidual()` and `computeQpJacobian()`.

#### The Solution Variable

| Variable | Type | Description |
|----------|------|-------------|
| `_u` | `const VariableValue &` | Values of the solution variable at all quadrature points of the current element. Index with `_u[_qp]` to get a `Real`. |
| `_grad_u` | `const VariableGradient &` | Gradient of the solution variable at all quadrature points. Index with `_grad_u[_qp]` to get a `RealGradient` (3-component vector). |

`_u` and `_grad_u` hold the values from the current Newton iterate. During a transient solve,
these are the values at the current time step. If you need values from the previous time step, use
`_u_old[_qp]` and `_u_older[_qp]` (available via `Coupleable`).

#### Test Functions (Rows of the System)

| Variable | Type | Description |
|----------|------|-------------|
| `_test` | `const VariableTestValue &` | Test function values. `_test[_i][_qp]` is the value of test function `i` at quadrature point `qp`. |
| `_grad_test` | `const VariableTestGradient &` | Test function gradients. `_grad_test[_i][_qp]` is the gradient (a `RealGradient`) of test function `i` at quadrature point `qp`. |

The outer loop over `_i` (the test function index) corresponds to a row of the residual vector.
MOOSE runs this loop automatically. Inside `computeQpResidual()` the current row index is `_i`.

#### Shape Functions (Columns of the System / How u Is Represented)

| Variable | Type | Description |
|----------|------|-------------|
| `_phi` | `const VariablePhiValue &` | Shape function values. `_phi[_j][_qp]` is the value of shape function `j` at quadrature point `qp`. |
| `_grad_phi` | `const VariablePhiGradient &` | Shape function gradients. `_grad_phi[_j][_qp]` is the gradient of shape function `j` at quadrature point `qp`. |

The inner loop over `_j` (shape function index) corresponds to a column of the Jacobian matrix.
Inside `computeQpJacobian()` the current column index is `_j`. Note that for the same variable,
the test functions and shape functions are often identical (both are the same basis functions), but
they serve different conceptual roles.

#### Quadrature Point Geometry

| Variable | Type | Description |
|----------|------|-------------|
| `_qp` | `unsigned int` | Current quadrature point index. Always use `_qp` as the array index; it ranges from 0 to `_qrule->n_points() - 1`. |
| `_q_point` | `const MooseArray<Point> &` | Physical (x,y,z) coordinates of each quadrature point. `_q_point[_qp]` is a `Point` object with `.x()`, `.y()`, `.z()` accessors. Use this when the integrand explicitly depends on position (e.g., manufactured solutions). |
| `_JxW` | `const MooseArray<Real> &` | Quadrature weight multiplied by the Jacobian determinant of the mapping from reference element to physical element. This is the "volume" weight for each quadrature point. MOOSE multiplies your returned residual by `_JxW[_qp] * _coord[_qp]` after your method returns. |
| `_coord` | `const MooseArray<Real> &` | Coordinate-system scaling. For Cartesian problems this is always 1.0. For axisymmetric (RZ) problems it is the radial coordinate `r`. For spherical it is `r^2 * sin(theta)`. |
| `_i` | `unsigned int` | Current test function index (row). Available in `computeQpResidual()` and `computeQpJacobian()`. |
| `_j` | `unsigned int` | Current shape function index (column). Available only in `computeQpJacobian()`. |

#### The Implicit Quadrature Point Loop

You never write a loop over `_qp`. MOOSE's assembly system calls your method once per quadrature
point, incrementing `_qp` between calls. Here is the conceptual pseudo-code that MOOSE runs
(simplified):

```
for each element in mesh:
    for each quadrature point qp in element:
        _qp = qp
        for each test function index i:
            _i = i
            residual[dof(i)] += computeQpResidual() * _test[i][qp] * _JxW[qp] * _coord[qp]
            for each shape function index j:
                _j = j
                jacobian[dof(i)][dof(j)] += computeQpJacobian() * _test[i][qp] * _phi[j][qp]
                                            * _JxW[qp] * _coord[qp]
```

Wait - look again at `computeQpResidual()` for `ReactionDiffusion`. The diffusion term already
includes `_grad_test[_i][_qp]`, and the reaction term already includes `_test[_i][_qp]`. Then
the loop multiplies again by `_test[_i][_qp]`? That would double-count it.

The answer is that the convention in MOOSE is: **your `computeQpResidual()` returns the full
integrand, including the test function factor**. MOOSE then multiplies only by `_JxW[_qp] *
_coord[_qp]` - not by `_test`. Look at the `Diffusion` kernel in the framework for confirmation:

```cpp
Real
Diffusion::computeQpResidual()
{
  return _grad_u[_qp] * _grad_test[_i][_qp];  // includes _grad_test
}

Real
Diffusion::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];  // includes _grad_test
}
```

The `_test` and `_grad_test` factors ARE in your expression. MOOSE applies `_JxW * _coord` on top
of your result, but the test function part is YOUR responsibility. This is the correct interpretation.

### Step 3: Writing the Input File

Create `test/tests/reaction_diffusion/reaction_diffusion_test.i`:

```ini
# Solve -div(0.5 * grad(u)) + 2.0 * u = 0
# on a 2D unit square with Dirichlet BCs

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [u]
    # Default family = LAGRANGE, order = FIRST
    # These are first-order Lagrange (bilinear) basis functions.
  []
[]

[Kernels]
  [reaction_diffusion]
    type = ReactionDiffusion     # Must match the class name exactly
    variable = u                 # The variable this kernel acts on
    D = 0.5                      # Diffusion coefficient
    lambda = 2.0                 # Reaction rate
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'          # Preconditioned Jacobian-Free Newton-Krylov
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'hypre'
[]

[Outputs]
  exodus = true
[]
```

Run this with:

```bash
cd test
../my_app-opt -i tests/reaction_diffusion/reaction_diffusion_test.i
```

---

## 5. Writing a Material

### When to Use Material vs Hardcoding in a Kernel

Hardcode a coefficient directly in a Kernel when:
- It is a simple constant that will never vary spatially or in time
- No other object ever needs the same value

Use a Material when:
- The property is a function of the solution (e.g., temperature-dependent conductivity)
- Multiple Kernels or BCs need the same property
- You want users to be able to swap constitutive models without changing Kernel code
- The property depends on position, element ID, or external data

The Material/Kernel separation mirrors the separation of concerns in object-oriented programming:
Materials describe what the material IS (its physical properties), Kernels describe what the
governing equation says (how those properties enter the math).

### Declaring and Consuming Properties

A Material class declares properties with `declareProperty<T>()` and fills them in
`computeQpProperties()`. A Kernel (or other object) consumes them with `getMaterialProperty<T>()`.

The string name passed to both calls must match exactly. It is a runtime lookup, not a
compile-time check. A typo in the string causes a MOOSE error before the solve begins.

### Full Example: Temperature-Dependent Thermal Conductivity

We implement a material where `k(T) = k0 + k1 * T`.

#### Header: `include/materials/ThermalConductivity.h`

```cpp
#pragma once

#include "Material.h"

class ThermalConductivity : public Material
{
public:
  static InputParameters validParams();

  ThermalConductivity(const InputParameters & parameters);

protected:
  // This is the method you override to fill properties.
  // It is called for every quadrature point on every element.
  virtual void computeQpProperties() override;

  // Parameters
  const Real _k0;     // baseline conductivity at T = 0
  const Real _k1;     // linear temperature coefficient

  // The temperature variable - we read it but do not own it.
  // coupledValue() returns a const VariableValue & that MOOSE fills
  // with the current solution values before calling computeQpProperties().
  const VariableValue & _temperature;

  // The property we declare (write to).
  // declareProperty returns a writable MaterialProperty<Real> &.
  // We store it as a reference so we can assign to it in computeQpProperties.
  MaterialProperty<Real> & _conductivity;

  // If another object needs d(conductivity)/d(temperature) for the Jacobian,
  // declare a second property for it.
  MaterialProperty<Real> & _d_conductivity_dT;
};
```

#### Source: `src/materials/ThermalConductivity.C`

```cpp
#include "ThermalConductivity.h"

registerMooseObject("MyAppApp", ThermalConductivity);

InputParameters
ThermalConductivity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes k(T) = k0 + k1*T for use in heat conduction kernels.");
  params.addRequiredParam<Real>("k0", "Baseline conductivity (W/m/K)");
  params.addParam<Real>("k1", 0.0, "Linear temperature coefficient (W/m/K^2)");
  // addCoupledVar tells MOOSE that this material reads a solution variable.
  // The user writes  temperature = T  in the input file block.
  params.addCoupledVar("temperature", "The temperature variable");
  return params;
}

ThermalConductivity::ThermalConductivity(const InputParameters & parameters)
  : Material(parameters),
    _k0(getParam<Real>("k0")),
    _k1(getParam<Real>("k1")),
    // coupledValue returns a const reference to the array of values at
    // quadrature points for the coupled variable. The string "temperature"
    // is the PARAMETER name (declared in validParams with addCoupledVar),
    // not the variable name in the input file.
    _temperature(coupledValue("temperature")),
    // declareProperty creates a MaterialProperty and returns a writable ref.
    // The string "thermal_conductivity" is the property NAME that other
    // objects use in getMaterialProperty("thermal_conductivity").
    _conductivity(declareProperty<Real>("thermal_conductivity")),
    _d_conductivity_dT(declareProperty<Real>("d_thermal_conductivity_dT"))
{
}

void
ThermalConductivity::computeQpProperties()
{
  // _qp is the current quadrature point index, inherited from MaterialBase.
  // Assign to the property at the current quadrature point.
  _conductivity[_qp] = _k0 + _k1 * _temperature[_qp];
  _d_conductivity_dT[_qp] = _k1;
}
```

#### Consuming the Property in a Kernel

```cpp
// In HeatConductionKernel.h:
class HeatConductionKernel : public Kernel
{
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;

  // getMaterialProperty returns a const reference. You cannot write to it.
  const MaterialProperty<Real> & _conductivity;
  const MaterialProperty<Real> & _d_conductivity_dT;
};

// In HeatConductionKernel.C constructor:
HeatConductionKernel::HeatConductionKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _conductivity(getMaterialProperty<Real>("thermal_conductivity")),
    _d_conductivity_dT(getMaterialProperty<Real>("d_thermal_conductivity_dT"))
{
}

Real
HeatConductionKernel::computeQpResidual()
{
  return _conductivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}
```

#### Input File Block

```ini
[Materials]
  [k]
    type = ThermalConductivity
    k0 = 50.0
    k1 = 0.1
    temperature = T    # name of the variable in [Variables]
  []
[]
```

### Stateful Properties

Some physics require history. For example, plasticity needs the plastic strain from the previous
time step, and explicit damage models need the damage from the previous step.

Stateful properties have an "old" value (previous time step) and an "older" value (two steps ago).
To use them, declare the property normally with `declareProperty<T>()` and then also declare it
as stateful with `getMaterialPropertyOld<T>()` in the constructor:

```cpp
// In constructor:
_strain(declareProperty<Real>("plastic_strain")),
_strain_old(getMaterialPropertyOld<Real>("plastic_strain"))

// In computeQpProperties():
void
PlasticMaterial::computeQpProperties()
{
  // _strain_old[_qp] is the value from the PREVIOUS timestep.
  // Compute the new value based on the old one:
  _strain[_qp] = _strain_old[_qp] + computeIncrementalStrain();
}
```

You also need to implement `initQpStatefulProperties()` to provide initial values:

```cpp
void
PlasticMaterial::initQpStatefulProperties()
{
  _strain[_qp] = 0.0;  // zero initial plastic strain
}
```

`initQpStatefulProperties()` is called once before the first time step to initialise the "old"
value storage.

### AD Materials

`ADMaterial` is a type alias for `Material` (see `framework/include/materials/ADMaterial.h`). There
is no separate base class; the distinction is in how you declare and use properties:

```cpp
// Declare an AD property (value carries derivative information)
ADMaterialProperty<Real> & _conductivity = declareADProperty<Real>("thermal_conductivity");

// Consume an AD property in an ADKernel
const ADMaterialProperty<Real> & _conductivity = getADMaterialProperty<Real>("thermal_conductivity");
```

When you use AD materials with AD kernels, the derivatives are propagated automatically through
the material evaluation, giving you accurate Jacobian contributions from the material model
without any hand derivation. See section 9 for a full treatment of AD.

---

## 6. Writing Boundary Conditions

MOOSE provides two families of boundary conditions. The choice depends on whether the condition
is enforced pointwise at nodes or integrated over boundary faces.

### NodalBC: Pointwise Boundary Conditions

`NodalBC` is used for strong (Dirichlet) conditions. At each node on the specified boundary, MOOSE
replaces the residual equation for that DOF with your expression. The solve then drives your
expression to zero, which enforces the condition exactly.

The `DirichletBC` built into MOOSE is the canonical example: it returns `_u[_qp] - value`, so the
solver drives `u - value -> 0`, meaning `u -> value`.

Key members inherited from `NodalBC`:

| Member | Type | Description |
|--------|------|-------------|
| `_u` | `const VariableValue &` | Value of the variable at the current node. Index with `_u[_qp]` where `_qp` is always 0 for nodal BCs. |
| `_current_node` | `const Node * const &` | Pointer to the current libMesh node object. Access coordinates via `(*_current_node)(0)` (x), `(*_current_node)(1)` (y), `(*_current_node)(2)` (z). |
| `_qp` | `const unsigned int` | Always 0 for NodalBC. The "quadrature point" concept degenerates to the node itself. |

#### Example: Dirichlet BC with prescribed value

```cpp
// Header
class MyDirichletBC : public NodalBC
{
public:
  static InputParameters validParams();
  MyDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  const Real _value;
};

// Source
InputParameters
MyDirichletBC::validParams()
{
  InputParameters params = NodalBC::validParams();
  params.addClassDescription("Dirichlet BC that sets u = value on a boundary.");
  params.addParam<Real>("value", 0.0, "Value to enforce at the boundary.");
  return params;
}

MyDirichletBC::MyDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters), _value(getParam<Real>("value"))
{
}

Real
MyDirichletBC::computeQpResidual()
{
  // NodalBCs enforce residual == 0 at the node.
  // Return (u - prescribed_value). When the solver drives residual to zero,
  // it drives u to the prescribed value.
  return _u[_qp] - _value;
}
```

#### Input file block

```ini
[BCs]
  [left_wall]
    type = MyDirichletBC
    variable = u
    boundary = left     # boundary name from the mesh
    value = 0.0
  []
[]
```

The default `computeQpJacobian()` in `NodalBC` returns 1, which is correct for the identity
condition above (`d(u - value)/d(u) = 1`). Override it if your residual has a nonlinear
dependence on `_u`.

### IntegratedBC: Boundary-Integrated Conditions

`IntegratedBC` is used for natural (Neumann) conditions and flux boundary conditions. These terms
are assembled via integration over boundary faces, exactly like element assembly but restricted to
boundary sides.

Key members inherited from `IntegratedBC` (beyond the standard Kernel members):

| Member | Type | Description |
|--------|------|-------------|
| `_normals` | `const MooseArray<Point> &` | Outward unit normal vectors at boundary quadrature points. `_normals[_qp]` is a `Point` object. Use it for flux terms like `flux * _normals[_qp]`. |
| `_u`, `_grad_u` | same as Kernel | Solution values and gradients on the boundary face. |
| `_test`, `_grad_test` | same as Kernel | Test function values and gradients on the boundary face. |
| `_phi`, `_grad_phi` | same as Kernel | Shape function values and gradients on the boundary face. |

#### Example: Neumann (constant flux) BC

The Neumann condition applies a specified outward flux: `-D * grad(u) . n = g` where `g` is a
given flux value. After moving to the left side of the weak form, the boundary integral term is
`-g * test_i`.

```cpp
class NeumannBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  NeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  const Real _flux;   // prescribed outward flux magnitude
};

InputParameters
NeumannBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Applies a constant outward flux -D*grad(u).n = flux on a boundary.");
  params.addParam<Real>("flux", 0.0, "Outward flux magnitude.");
  return params;
}

NeumannBC::NeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _flux(getParam<Real>("flux"))
{
}

Real
NeumannBC::computeQpResidual()
{
  // The boundary integral contribution to the residual.
  // Sign convention: this subtracts the boundary term from the domain integral.
  // A positive flux value means material is leaving the domain.
  return -_flux * _test[_i][_qp];
}
```

#### Example: Convective (Robin) BC

A convective BC applies `h * (u - T_inf)` as a heat flux:

```cpp
Real
ConvectiveFluxBC::computeQpResidual()
{
  // h * (u - T_ambient) * test_i
  // When positive, this drains heat from the boundary (u > T_ambient means outflow).
  return _h * (_u[_qp] - _T_ambient) * _test[_i][_qp];
}

Real
ConvectiveFluxBC::computeQpJacobian()
{
  // d/du_j [ h * (u - T_ambient) * test_i ]
  // = h * phi_j * test_i
  return _h * _phi[_j][_qp] * _test[_i][_qp];
}
```

#### Input file block

```ini
[BCs]
  [right_flux]
    type = NeumannBC
    variable = u
    boundary = right
    flux = 10.0
  []
  [top_convective]
    type = ConvectiveFluxBC
    variable = u
    boundary = top
    h = 5.0            # convection coefficient
    T_ambient = 25.0   # far-field temperature
  []
[]
```

---

## 7. Writing a Postprocessor

Postprocessors are objects that compute a single scalar value from the simulation state. The value
is typically written to a CSV file and printed to the console. Examples include the total
integrated heat flux, the maximum temperature, the L2 norm of the solution, or the number of
nonlinear iterations.

### Four Families of Postprocessor

| Base Class | Loop Type | Override Method |
|------------|-----------|----------------|
| `ElementPostprocessor` | Over elements and their quadrature points | `execute()`, or inherit `ElementIntegralPostprocessor` and override `computeQpIntegral()` |
| `NodalPostprocessor` | Over nodes | `execute()` |
| `SidePostprocessor` | Over boundary faces | `execute()` |
| `GeneralPostprocessor` | No geometric loop; you compute whatever you want | `execute()` |

All four families require you to implement:
- `initialize()`: Reset accumulators to zero before each execution
- `execute()`: Accumulate data for the current geometric entity
- `finalize()`: Perform MPI reduction (see section 8 for details)
- `getValue() const`: Return the scalar result

The cleanest way to write an element-integral postprocessor is to inherit from
`ElementIntegralPostprocessor` (or `ElementIntegralVariablePostprocessor` if you already have a
coupled variable). These base classes provide `initialize`, `execute`, `threadJoin`, and `finalize`
implementations that handle the accumulator pattern correctly. You only override
`computeQpIntegral()`.

### Full Example: L2 Norm of Error

This postprocessor computes `||u - u_exact||_{L2} = sqrt( integral( (u - u_exact)^2 ) )`.
The built-in `ElementL2Error` in the framework does exactly this. We reproduce it here
with annotations.

#### Header: `include/postprocessors/SolutionL2Error.h`

```cpp
#pragma once

// ElementIntegralVariablePostprocessor already handles the loop over
// elements and quadrature points and accumulates the integral.
// It provides _u[_qp] and _q_point[_qp] and the FunctionInterface.
#include "ElementIntegralVariablePostprocessor.h"

// Forward declaration of Function to avoid including the full header
// in this header file. The full include goes in the .C file.
class Function;

class SolutionL2Error : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  SolutionL2Error(const InputParameters & parameters);

  // getValue() takes the sqrt of the accumulated integral.
  virtual Real getValue() const override;

protected:
  // computeQpIntegral() returns the integrand at one quadrature point.
  // The base class accumulates sum_qp( computeQpIntegral() * _JxW[_qp] * _coord[_qp] ).
  virtual Real computeQpIntegral() override;

  // Reference to the exact solution function object.
  const Function & _exact_solution;
};
```

#### Source: `src/postprocessors/SolutionL2Error.C`

```cpp
#include "SolutionL2Error.h"
#include "Function.h"

registerMooseObject("MyAppApp", SolutionL2Error);

InputParameters
SolutionL2Error::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addClassDescription("Computes the L2 norm of the error against an analytic solution.");
  // addRequiredParam<FunctionName> registers a parameter that holds the name of a [Functions]
  // block. MOOSE resolves it to a Function object reference automatically.
  params.addRequiredParam<FunctionName>("exact_solution",
                                       "The analytic solution to compare against.");
  return params;
}

SolutionL2Error::SolutionL2Error(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters),
    // getFunction() looks up the Function object by the name given in the input file.
    // It returns a const reference. You do not manage its lifetime.
    _exact_solution(getFunction("exact_solution"))
{
}

Real
SolutionL2Error::getValue() const
{
  // The base class accumulates integral( integrand ), so here it holds
  // integral( (u - u_exact)^2 ). Taking the sqrt gives the L2 norm.
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real
SolutionL2Error::computeQpIntegral()
{
  // _u[_qp]: solution at current quadrature point (from coupled variable)
  // _exact_solution.value(_t, _q_point[_qp]): analytic value at same point
  // _t: current time (available via TransientInterface)
  // _q_point[_qp]: physical coordinates of the quadrature point
  Real diff = _u[_qp] - _exact_solution.value(_t, _q_point[_qp]);
  return diff * diff;
}
```

#### Input File Block

```ini
[Functions]
  [exact]
    type = ParsedFunction
    expression = 'sin(pi*x) * sin(pi*y)'
  []
[]

[Postprocessors]
  [l2_error]
    type = SolutionL2Error
    variable = u
    exact_solution = exact
    execute_on = 'TIMESTEP_END'   # when to compute this postprocessor
  []
[]

[Outputs]
  csv = true      # writes postprocessor values to a CSV file
  [screen]
    type = Console
    postprocessor_screen_output = true
  []
[]
```

### CSV Output and Python Analysis

When `csv = true` is set in `[Outputs]`, MOOSE writes a file `<input_name>_out.csv` with one row
per output time and one column per postprocessor. The first column is always `time`.

Example CSV content:
```
time,l2_error
0,0
0.1,0.00234
0.2,0.00189
1.0,0.00012
```

Read it in Python:

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('my_simulation_out.csv')
plt.semilogy(df['time'], df['l2_error'])
plt.xlabel('Time')
plt.ylabel('L2 Error')
plt.title('Convergence History')
plt.savefig('convergence.png')
```

---

## 8. Writing a UserObject

### Difference from Postprocessor

A Postprocessor computes a single scalar and writes it to CSV. A UserObject is more general: it
can compute any data structure, provide spatial queries (`spatialValue()`), and communicate with
other objects during the solve.

Use a UserObject when:
- You need to store a complex data structure (not just one number)
- Other simulation objects query you during the solve for data
- You need fine control over execution timing
- You are implementing a lookup table, a spatial field from experimental data, or a custom
  inter-element communication pattern

The boundary between UserObject and Postprocessor blurs because `ElementPostprocessor` actually
inherits from `ElementUserObject`. A Postprocessor IS a UserObject with a `getValue()` method.

### Execute-On Timing

Every UserObject (and many other MOOSE objects) has an `execute_on` parameter that controls when
it runs. The available flags, declared in `framework/include/base/Moose.h`, are:

| Flag | When it fires |
|------|--------------|
| `EXEC_INITIAL` | Once at the beginning, before the first time step |
| `EXEC_TIMESTEP_BEGIN` | At the start of each time step, before the solve |
| `EXEC_TIMESTEP_END` | At the end of each time step, after convergence |
| `EXEC_LINEAR` | After each linear solve within each Newton iteration |
| `EXEC_NONLINEAR` | After each Newton iteration |
| `EXEC_FINAL` | Once after the last time step |
| `EXEC_ALWAYS` | Every time any of the above fires |
| `EXEC_CUSTOM` | Only when explicitly triggered by code |

Multiple flags can be combined: `execute_on = 'INITIAL TIMESTEP_END'`.

### Writing a Thread-Safe UserObject

When running with multiple threads (`--n-threads 4`), MOOSE spawns multiple copies of your
UserObject (one per thread). Each thread processes a subset of elements. After all threads finish,
MOOSE calls `threadJoin(const UserObject & other)` on the primary thread's copy, passing each
other copy as `other`. Your `threadJoin` must merge `other`'s data into `this`.

After all `threadJoin` calls complete, MOOSE calls `finalize()`. This is the correct place to
perform MPI communication to combine results across MPI ranks.

#### Full Example: Volume-Weighted Average of a Variable

```cpp
// Header
class VolumeWeightedAverage : public ElementUserObject
{
public:
  static InputParameters validParams();
  VolumeWeightedAverage(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  // Other objects can call this after execution to get the result.
  Real getValue() const { return _total_volume > 0.0 ? _sum / _total_volume : 0.0; }

protected:
  const VariableValue & _variable;

  // Accumulators
  Real _sum;           // sum of (variable_value * element_volume)
  Real _total_volume;  // total mesh volume
};

// Source
void
VolumeWeightedAverage::initialize()
{
  // Reset at the beginning of each execution pass.
  _sum = 0.0;
  _total_volume = 0.0;
}

void
VolumeWeightedAverage::execute()
{
  // Loop over quadrature points on the current element.
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    _sum += _variable[qp] * _JxW[qp] * _coord[qp];
    _total_volume += _JxW[qp] * _coord[qp];
  }
}

void
VolumeWeightedAverage::threadJoin(const UserObject & y)
{
  // Cast to this class to access its private accumulators.
  const VolumeWeightedAverage & other = static_cast<const VolumeWeightedAverage &>(y);
  _sum += other._sum;
  _total_volume += other._total_volume;
}

void
VolumeWeightedAverage::finalize()
{
  // _communicator is the MPI communicator, available from MooseObject.
  // sum() reduces across all MPI ranks using MPI_SUM.
  gatherSum(_sum);
  gatherSum(_total_volume);
}
```

Note that `gatherSum()` (provided by the parallel communication interface) handles both the
case of serial runs (no-op) and parallel MPI runs (calls `MPI_Allreduce` internally).

#### Input File Block

```ini
[UserObjects]
  [avg_temperature]
    type = VolumeWeightedAverage
    variable = T
    execute_on = 'TIMESTEP_END'
  []
[]
```

#### Consuming a UserObject in Another Object

```cpp
// In the constructor of a Kernel or BC that needs the average:
const VolumeWeightedAverage & _avg_T = getUserObject<VolumeWeightedAverage>("avg_temperature_uo");

// In computeQpResidual():
Real avg = _avg_T.getValue();
```

The input file parameter `avg_temperature_uo` maps to the name of the `[UserObjects]` block entry.

---

## 9. Automatic Differentiation (AD)

### The Problem AD Solves

Writing `computeQpJacobian()` by hand is tedious and error-prone. For complex material models
with many coupled variables, deriving the Jacobian analytically takes hours and mistakes produce
slow convergence or solver divergence. The errors are hard to find.

Automatic Differentiation (AD) solves this by computing derivatives automatically, exactly (up to
floating point precision), without symbolic manipulation.

### Dual Numbers

MOOSE's AD system is built on "dual numbers" (also called forward-mode AD or scalar-valued
derivatives). A dual number carries two parts:

```
ADReal x = { value: 3.14,
             derivatives: [dx/da, dx/db, dx/dc, ...] }
```

The derivative part is a fixed-size array with one entry per degree of freedom in the current
element. Every arithmetic operation on `ADReal` objects propagates both the value AND the
derivatives using the chain rule automatically. This is done at compile time through operator
overloading; there is no runtime interpretation.

Example: if `a` and `b` are `ADReal` variables with derivative arrays, then `a * b` produces a
new `ADReal` whose derivative array holds `b * da + a * db` by the product rule - automatically.

By the time you call `computeQpResidual()` and return an `ADReal`, its derivative array contains
the exact Jacobian entries for all DOFs in the element. MOOSE extracts those entries and populates
the Jacobian matrix without any further action from you.

### ADReal and ADKernel

The `ADKernel` base class (a typedef for `ADKernelTempl<Real>`) is a drop-in replacement for
`Kernel` with these differences:

1. `computeQpResidual()` returns `ADReal` instead of `Real`
2. You do not implement `computeQpJacobian()` - it is computed automatically
3. `_u`, `_grad_u` are `ADReal` (or arrays of `ADReal`) instead of `Real`
4. Material properties must be consumed with `getADMaterialProperty` to carry derivatives

#### Full Example: AD version of Diffusion

```cpp
// Header
#pragma once
#include "ADKernel.h"

class ADDiffusionKernel : public ADKernel
{
public:
  static InputParameters validParams();
  ADDiffusionKernel(const InputParameters & parameters);

protected:
  // Returns ADReal - automatically differentiable
  virtual ADReal computeQpResidual() override;
};

// Source
ADDiffusionKernel::ADDiffusionKernel(const InputParameters & parameters)
  : ADKernel(parameters)
{
}

ADReal
ADDiffusionKernel::computeQpResidual()
{
  // _grad_u[_qp] is now ADRealVectorValue (a 3-component vector of ADReal).
  // The dot product propagates all derivative information automatically.
  // No computeQpJacobian() needed.
  return _grad_u[_qp] * _grad_test[_i][_qp];
}
```

Compare this with the non-AD version in `framework/src/kernels/Diffusion.C`:

```cpp
Real Diffusion::computeQpResidual() { return _grad_u[_qp] * _grad_test[_i][_qp]; }
Real Diffusion::computeQpJacobian() { return _grad_phi[_j][_qp] * _grad_test[_i][_qp]; }
```

The AD version eliminates `computeQpJacobian()` entirely. For simple kernels this saves little
effort. For kernels with coupled variables, material properties, and multiple terms, the savings
are dramatic.

### Why AD Eliminates Hand-Coded Jacobians

Consider a kernel with a nonlinear material property `k(u)` (conductivity depending on the
solution variable). The residual is `k(u) * grad(u) . grad(test)`. The Jacobian requires:

```
d/du_j [ k(u) * grad(u) . grad(test) ]
= dk/du * phi_j * grad(u) . grad(test) + k(u) * grad(phi_j) . grad(test)
```

You must derive `dk/du` analytically and implement it. If the material model changes, you update
two places (the material property AND the Jacobian). With AD, the derivative propagates through
`k(u)` automatically as long as `k` is also computed using `ADReal` arithmetic.

### PJFNK Fallback

MOOSE supports Jacobian-Free Newton-Krylov (JFNK) in the executioner:

```ini
[Executioner]
  type = Steady
  solve_type = PJFNK     # Uses preconditioner but approximates Jacobian-vector products
[]
```

With PJFNK, MOOSE approximates the Jacobian-vector product `J * v` using finite differences of
the residual: `(F(u + eps*v) - F(u)) / eps`. This avoids computing the Jacobian explicitly.
The solve still converges but:
- Each Newton iteration requires more linear iterations (no exact Jacobian information)
- The finite-difference approximation introduces noise that can stall convergence

PJFNK is a useful fallback when you cannot compute the Jacobian (e.g., when calling external
legacy code) but AD is almost always the better solution for new kernels.

### Derivative Array Size

The AD system pre-allocates a fixed-size derivative array per `ADReal`. The size is set at
MOOSE compile time and defaults to 50. This covers most single-physics problems (50 DOFs per
element is typical for second-order 3D hexahedral elements with 2-3 variables).

If you have many coupled variables or high-order elements, you may need to increase this:

```bash
# During libMesh configuration:
./configure --with-derivative-size=100
```

Setting the derivative size too large wastes memory and slows down AD slightly. Setting it too
small causes a runtime error when MOOSE tries to add more derivatives than the array can hold.
Check the number of DOFs per element in your problem and set accordingly.

---

## 10. The Build System in Depth

### How METHOD Changes Compiler Flags

`build.mk` queries libmesh-config with the current `METHOD` to get compiler flags:

```makefile
libmesh_CXXFLAGS := $(shell METHOD=$(METHOD) $(libmesh_config) --cxxflags)
```

libMesh was itself built with the same METHOD options, so the flags match. Typical flag sets:

| METHOD | Key flags |
|--------|----------|
| opt    | `-O2 -DNDEBUG` - disables `assert()` and MOOSE's internal assertions |
| dbg    | `-O0 -g -DDEBUG` - enables full assertion checking, no optimisation |
| devel  | `-O2 -g -DDEBUG` - optimised but with assertions and debug symbols |

When METHOD=opt, Eigen assertions are also disabled. This means code that accesses a vector out of
bounds will corrupt memory silently instead of crashing with a useful message. Always develop
with METHOD=dbg.

### How app.mk Discovers Sources Automatically

`app.mk` uses shell `find` to collect all `.C` files under `src/`:

```makefile
SRC_DIRS := $(APPLICATION_DIR)/src
srcfiles  := $(shell find $(SRC_DIRS) -regex "[^\#~]*\.C" $(find_excludes))
```

The regex `[^\#~]*\.C` matches any path ending in `.C` that does not contain `#` or `~`
(editor backup files). The `.C` extension is uppercase to distinguish MOOSE C++ source files
from C files (`.c`).

Object files are placed in a build directory:

```
build/
  objects/
    opt/
      src/
        kernels/
          ReactionDiffusion.opt.o
```

The METHOD appears in the object file suffix, so you can build multiple methods without
cleaning between them.

### Adding External Libraries

To link an external library, add to your application's `Makefile` before the `include app.mk` line:

```makefile
# Include path for headers
ADDITIONAL_INCLUDES := -I/path/to/external/include

# Library to link
ADDITIONAL_LIBS := -L/path/to/external/lib -lexternal

# Pass these to app.mk
APP_INCLUDES += $(ADDITIONAL_INCLUDES)
LIBS += $(ADDITIONAL_LIBS)
```

For header-only libraries, only `ADDITIONAL_INCLUDES` is needed.

If the library requires special compiler flags (e.g., OpenMP or CUDA), add them via:

```makefile
libmesh_CXXFLAGS += -fopenmp
```

### conf_vars.mk for Optional Features

When you run `./configure` in the MOOSE root (optional but useful for controlling which optional
packages are found), it writes a file `conf_vars.mk` that is automatically included by
`build.mk`:

```makefile
-include $(MOOSE_DIR)/conf_vars.mk
```

This file sets variables like `HAVE_SLEPC`, `HAVE_SUPERLU`, or `PETSC_HAVE_HDF5` to yes/no based
on what was found during configure. Application Makefiles can test these:

```makefile
ifeq ($(HAVE_SLEPC),yes)
  # enable eigenvalue solver features
endif
```

Without running configure, `conf_vars.mk` does not exist and the `-include` silently succeeds with
all optional features disabled.

---

## 11. Testing

MOOSE has a rigorous test harness. Every new feature should have at least one test before it is
considered finished. Tests are the safety net that lets you refactor with confidence.

### Test Spec Format

Tests live in `test/tests/` and are organised into subdirectories by feature. Each directory
has a file named `tests` (no extension) in HIT format (the same format as input files):

```ini
[Tests]
  [my_kernel_test]
    type = Exodiff                            # test class from TestHarness
    input = my_kernel.i                       # input file relative to this directory
    exodiff = my_kernel_out.e                 # Exodus output file to compare
    issues = '#123'                           # GitHub issue numbers
    design = 'kernels/MyKernel.md'           # documentation file path
    requirement = 'The kernel shall correctly compute ... '  # SRS requirement text
  []
[]
```

The `issues`, `design`, and `requirement` fields are used by MOOSE's software quality assurance
system. They are required for framework contributions but optional for application development.

#### Test Types

| Type | What it does |
|------|-------------|
| `Exodiff` | Compares ExodusII output against a gold file using `exodiff` |
| `CSVDiff` | Compares CSV output against a gold file |
| `RunApp` | Runs the app and checks the exit code is 0 |
| `RunException` | Runs the app and verifies a specific error message appears |
| `PetscJacobianTester` | Runs the app with PETSc's Jacobian checker and verifies tolerance |
| `AnalyzeJacobian` | Uses MOOSE's built-in finite-difference Jacobian checker |

### Gold Files and Regeneration

Gold files are the accepted-correct outputs stored in a `gold/` subdirectory. The first time you
create a test, you run the simulation manually and copy the output to `gold/`:

```bash
# Run the simulation
./my_app-opt -i test/tests/my_feature/my_test.i

# The output file is in the current directory (or as specified in [Outputs])
# Copy it to the gold directory
cp my_test_out.e test/tests/my_feature/gold/my_test_out.e
```

After a change that intentionally changes the output (e.g., a bug fix that changes a result),
regenerate the gold file:

```bash
./run_tests --re=my_test --update-all
```

This re-runs the test, takes the new output as correct, and overwrites the gold file. Review the
diff before committing.

### RunException for Error Testing

`RunException` verifies that your code fails with the correct error message:

```ini
[test_bad_input]
  type = RunException
  input = bad_input.i
  expect_err = 'D must be positive'   # regex pattern matched against stderr
[]
```

The `expect_err` value is a Python regular expression. The test passes if the executable exits
with a nonzero exit code AND the expected pattern appears in its error output. This ensures your
validation messages remain correct as code evolves.

### PetscJacobianTester for Jacobian Verification

Before committing a new Kernel, verify the Jacobian is correct. MOOSE provides a built-in
mechanism via PETSc's Jacobian check:

```ini
# In the input file, add to [Executioner]:
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options = '-snes_test_jacobian'     # Enables the check
  petsc_options_iname = '-snes_test_jacobian_view'
  petsc_options_value = ''
[]
```

Or use the `PetscJacobianTester` test type:

```ini
[jacobian_test]
  type = PetscJacobianTester
  input = my_kernel.i
  ratio_tol = 1e-7         # tolerance for ratio of AD Jacobian to FD Jacobian
  difference_tol = 1e0     # tolerance for absolute difference
[]
```

PETSc computes the Jacobian twice: once using your analytical/AD `computeQpJacobian()` and once
using finite differences of the residual. It then compares the two. A ratio near 1.0 and a
difference near 0.0 mean your Jacobian is correct.

### Running Specific Tests

```bash
# Run all tests in the current app
./run_tests -j8

# Run tests matching a regular expression (re is a Python regex)
./run_tests -j8 --re=reaction_diffusion

# Run tests matching a pattern in a specific subdirectory
./run_tests -j8 --re=kernels/

# Run tests tagged with a specific group
./run_tests -j8 --group=my_feature_group

# Run with a specific METHOD
./run_tests -j8 METHOD=dbg

# Dry run: shows which tests would run without running them
./run_tests --dry-run --re=simple_diffusion

# Update all gold files for matching tests
./run_tests --update-all --re=reaction_diffusion
```

The `--re` flag is a Python `re.search` pattern applied to the full test name, which is of the
form `path/to/tests:test_name`. Using `--re` is case-sensitive by default.

To run a single specific test case:

```bash
./run_tests --re='kernels/reaction_diffusion:my_kernel_test'
```

---

## 12. Code Style and Conventions

These conventions are enforced by `.clang-format`, code review, and in some cases by the build
system itself.

### .clang-format Settings

```yaml
BasedOnStyle: LLVM
TabWidth: 2
ColumnLimit: 100
UseTab: Never
BreakBeforeBraces: Allman        # opening brace on its own line
AlwaysBreakAfterDefinitionReturnType: TopLevel  # return type on its own line
PointerAlignment: Middle         # spaces on both sides: Real * ptr
SortIncludes: false              # preserve include order as written
IndentCaseLabels: true
```

Key points for a developer coming from Java or Python:
- **2-space indent**, not 4.
- **100-column line limit**, not 80.
- **Allman braces**: the opening `{` goes on its own line for functions, classes, and
  control structures. Inside a function, it is on its own line for control structures.
- **Return type on its own line** for top-level function definitions:
  ```cpp
  Real                                      // return type alone on this line
  MyKernel::computeQpResidual()             // function signature on next line
  {
    return _u[_qp] * _test[_i][_qp];
  }
  ```

### File Extensions

| Extension | Meaning |
|-----------|---------|
| `.h` | Header files (class declarations, inline methods) |
| `.C` | C++ source files (class definitions) |
| `.c` | Plain C source files (rare in MOOSE) |

The uppercase `.C` is intentional. Unix/macOS filesystems are case-sensitive; Windows is not.
When working on Windows within WSL or MSYS2, be careful: creating `myfile.c` when the build
system expects `MyFile.C` will cause missing-symbol errors that are hard to diagnose.

### Naming Conventions

| Element | Convention | Example |
|---------|-----------|---------|
| Class name | PascalCase | `ReactionDiffusion`, `ThermalConductivity` |
| Filename | Same as class, exact case | `ReactionDiffusion.h`, `ReactionDiffusion.C` |
| Member variable | `_camelCase` prefix with underscore | `_conductivity`, `_total_volume` |
| Local variable | `camelCase` | `localValue`, `errorNorm` |
| Method name | `camelCase` | `computeQpResidual`, `getTemperature` |
| Parameter name (in `validParams`) | `snake_case` | `"diffusion_coefficient"` |
| Constant | `ALL_CAPS_SNAKE_CASE` | `MAX_ITERATIONS` |

The underscore prefix for member variables is strictly enforced in code review. It prevents naming
collisions between member variables and local variables or parameters, and makes it immediately
visible in a method body which names refer to member state.

### registerMooseObject Placement

The `registerMooseObject` macro must appear exactly once per class, at the top of the `.C` file,
immediately after the includes and before the first function definition:

```cpp
#include "ReactionDiffusion.h"

registerMooseObject("MyAppApp", ReactionDiffusion);  // here, and only here

InputParameters
ReactionDiffusion::validParams()
{ ... }
```

Never put `registerMooseObject` in a header file. Including the header from multiple `.C` files
would result in multiple registrations and a runtime error.

### LGPL License Header

Every `.h` and `.C` file in a MOOSE application must begin with the LGPL license header:

```cpp
//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
```

For your own application's proprietary files, replace this with your organisation's license
header. Check with your legal team what is required.

### Class Name Must Match Filename

The class name declared in a header MUST be identical (including case) to the header filename
without the `.h` extension. This is because `registerMooseObject` uses the class name as the
type string in the factory, and the factory is case-sensitive. A mismatch compiles fine but
causes the input file parser to fail with "unknown object type".

---

## 13. The validParams() Chain

### Always Call Parent First

Every `validParams()` implementation must call the parent class version first:

```cpp
InputParameters
MyKernel::validParams()
{
  // This line is not optional. It adds all base-class parameters.
  InputParameters params = Kernel::validParams();

  // Now add your own parameters
  params.addRequiredParam<Real>("my_param", "Documentation string");

  return params;
}
```

If you forget the parent call, your object will be missing parameters like `variable`, `block`,
`use_displaced_mesh`, `save_in`, and dozens of others. The solve may fail with confusing errors.

### Parameter API

```cpp
// Required parameter - simulation errors before solve if missing
params.addRequiredParam<Real>("D", "Diffusion coefficient");

// Optional parameter with default
params.addParam<Real>("lambda", 0.0, "Reaction rate (defaults to 0)");

// Parameter with no default - the user may optionally provide it
// Use hasParam("name") to check before getParam in the constructor
params.addParam<std::string>("output_file", "Optional output filename");

// Class description - shown in --dump and documentation
params.addClassDescription("Implements reaction-diffusion.");

// Coupled variable - declares that this object reads a solution variable
params.addCoupledVar("temperature", "The temperature variable to couple to");

// Multiple coupled variables - user can list several variable names
params.addCoupledVar("other_vars", "Additional variables (list)");

// Material property name parameter
// Allows the user to specify which material property to read
// by name in the input file, rather than hardcoding the string
params.addParam<MaterialPropertyName>("conductivity",
    "thermal_conductivity",
    "Name of the thermal conductivity material property");
```

The `MaterialPropertyName` parameter type is preferable to a raw `std::string` because MOOSE
generates better documentation and enables property name checking.

### Dump Syntax

To see all available input syntax for your application, including every parameter of every object:

```bash
./my_app-opt --dump
```

This prints the complete HIT grammar for your application, showing every object type and every
parameter with its documentation, default value, and type. Use it to:
- Check your `addClassDescription` text appears correctly
- Verify parameter defaults are right
- See what inherited parameters are available

Filter to a specific object:

```bash
./my_app-opt --dump | grep -A 50 'type = ReactionDiffusion'
```

---

## 14. Writing a MultiApp

### When to Split Into Sub-Apps

Use MultiApps when:
- You have physics on different domains (fuel pin + coolant channel)
- You have physics at different scales that decouple except through boundary conditions
- You want to use different time steps or mesh refinements for different physics
- You want to re-use an existing standalone application as a sub-component

Do not use MultiApps for tightly coupled physics that must converge simultaneously on the same
mesh. Use multiple Kernels in one application for that.

### TransientMultiApp vs FullSolveMultiApp

| Class | Use case |
|-------|---------|
| `TransientMultiApp` | Sub-app is a transient simulation that advances time alongside the parent |
| `FullSolveMultiApp` | Sub-app completes its full solve (all time steps) every time it is executed |

`TransientMultiApp` is appropriate for co-simulation where parent and sub-app march through
time together, exchanging data at each step. `FullSolveMultiApp` is appropriate when the sub-app
solves a steady-state problem driven by boundary conditions that change each parent time step.

#### Input File Structure

```ini
# parent_app.i

[MultiApps]
  [neutronics]
    type = TransientMultiApp
    input_files = neutronics_sub.i    # path to sub-app input file
    positions = '0 0 0'               # (x,y,z) position of sub-app origin
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [send_temperature]
    type = MultiAppInterpolationTransfer
    direction = to_multiapp          # parent -> sub-app
    multi_app = neutronics
    source_variable = T             # variable in parent
    variable = temperature          # variable in sub-app
  []

  [receive_power]
    type = MultiAppInterpolationTransfer
    direction = from_multiapp        # sub-app -> parent
    multi_app = neutronics
    source_variable = power_density # variable in sub-app
    variable = q                    # variable in parent
  []
[]
```

### Writing a Transfer

Transfers are the objects that move data between parent and sub-apps. The most common built-in
transfer types:

| Transfer Type | What it does |
|---------------|-------------|
| `MultiAppInterpolationTransfer` | Spatially interpolates field variable values between different meshes |
| `MultiAppCopyTransfer` | Copies variable DOF values when both apps use the same mesh |
| `MultiAppPostprocessorTransfer` | Transfers a scalar postprocessor value |
| `MultiAppNearestNodeTransfer` | Maps values using nearest-node pairing |

To write a custom Transfer, inherit from `MultiAppTransfer` and implement `execute()`. Inside
`execute()`, access the parent problem via `_fe_problem` and sub-app problems via
`getFromMultiApp()` or `getToMultiApp()`.

### Picard Iteration Setup

Picard (fixed-point) iteration couples parent and sub-app by repeated alternating solves until
convergence. Configure it in the executioner:

```ini
[Executioner]
  type = Transient
  fixed_point_max_its = 10           # maximum Picard iterations per time step
  fixed_point_rel_tol = 1e-6         # convergence tolerance (relative change)
  fixed_point_abs_tol = 1e-8         # convergence tolerance (absolute change)
[]
```

When `fixed_point_max_its > 1`, MOOSE automatically runs the Picard loop: solve parent, transfer
to sub-app, solve sub-app, transfer back, check convergence, repeat.

---

## 15. Writing a MeshGenerator

MeshGenerators create or transform the mesh before the simulation begins. They run in a directed
acyclic graph order based on their declared dependencies. The final generator in the chain
produces the mesh that the simulation uses.

### Inheriting MeshGenerator

```cpp
// Header
#pragma once
#include "MeshGenerator.h"

class MyMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();
  MyMeshGenerator(const InputParameters & parameters);

  // Must override. Returns the generated or modified mesh.
  // The return type is std::unique_ptr<MeshBase>.
  virtual std::unique_ptr<MeshBase> generate() override;

protected:
  // Parameters
  const Real _length;
  const unsigned int _nx;

  // Optional: reference to an input mesh from another generator
  // This is how you declare a dependency on a previous generator.
  std::unique_ptr<MeshBase> & _input_mesh;
};
```

```cpp
// Source
#include "MyMeshGenerator.h"

registerMooseObject("MyAppApp", MyMeshGenerator);

InputParameters
MyMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addClassDescription("Generates a 1D uniform mesh of given length and resolution.");
  params.addRequiredParam<Real>("length", "Domain length");
  params.addRequiredParam<unsigned int>("nx", "Number of elements");
  // Optional input from another generator - makes this a dependent generator
  params.addParam<MeshGeneratorName>("input", "Name of an input mesh generator to modify");
  return params;
}

MyMeshGenerator::MyMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _length(getParam<Real>("length")),
    _nx(getParam<unsigned int>("nx")),
    // getMesh declares that this generator depends on the "input" parameter's generator.
    // MUST be captured by reference. The unique_ptr will be populated before generate() is called.
    _input_mesh(getMesh("input", /* allow_invalid = */ true))
{
}

std::unique_ptr<MeshBase>
MyMeshGenerator::generate()
{
  // If an input mesh was provided, start from it; otherwise create a new mesh.
  std::unique_ptr<MeshBase> mesh;
  if (_input_mesh)
    mesh = std::move(_input_mesh);
  else
    mesh = buildMeshBaseObject();

  // Build the mesh using libMesh API
  // ... mesh manipulation code ...

  return mesh;
}
```

### Chaining with `input = other_generator`

A generator that takes an `input` parameter will receive the output mesh of `other_generator`
before its own `generate()` is called. This creates a pipeline:

```ini
[Mesh]
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []

  [refined_mesh]
    type = RefineBlockGenerator
    input = base_mesh              # depends on base_mesh
    block = '0'
    refinement = 2
  []

  [final]
    type = MyMeshGenerator
    input = refined_mesh           # depends on refined_mesh
    length = 10.0
    nx = 20
  []
[]
```

The execution order is automatically determined from the dependency graph.

### Mesh Metadata

Generators can attach metadata to the mesh that persists through the simulation and restarts.
This is useful for communicating geometry information (element counts, length scales, boundary IDs)
to other objects:

```cpp
// In the constructor of the generator (not in generate()):
Real & _element_count = declareMeshProperty<Real>("element_count", 0.0);

// In generate():
_element_count = mesh->n_elem();
setMeshProperty<Real>("element_count", static_cast<Real>(mesh->n_elem()));
```

Other objects read mesh metadata via `getMeshProperty<Real>("element_count")` from the
`MeshMetaDataInterface`.

---

## 16. Debugging Techniques

### METHOD=dbg for Assertions

Always compile with `METHOD=dbg` when developing. This enables:
- All C++ `assert()` macros
- MOOSE's internal `mooseAssert()` calls (index bounds checks, state validation, etc.)
- libMesh's bounds checking on vectors and arrays
- Eigen's bounds checking on matrix operations

A crash in `METHOD=dbg` produces a message like:

```
Assertion failed: (_qp < _u.size()), function computeQpResidual,
file /path/to/ReactionDiffusion.C, line 42.
```

In `METHOD=opt` the same error silently corrupts memory and produces wrong answers or crashes
much later with an unrelated error message.

### --show-input for Parsed Input Dump

After MOOSE parses your input file it has a complete in-memory representation of all objects,
their types, and their parameters. Print this with:

```bash
./my_app-opt -i test.i --show-input
```

The output is a HIT-formatted dump of the parsed input. Use it to:
- Verify that default parameters have been applied correctly
- Check that your `addParam` defaults appear as expected
- Understand what MOOSE internally adds to blocks (coordinate systems, etc.)

### PETSc Options for Solver Diagnostics

PETSc has a rich set of runtime options for diagnosing solver problems. Pass them in the input
file or on the command line with the `-pc_type` style syntax.

#### Monitor Residual Convergence

```ini
[Executioner]
  petsc_options = '-snes_monitor -ksp_monitor'
[]
```

`-snes_monitor` prints the nonlinear (SNES) residual norm at each Newton iteration:
```
0 SNES Function norm 2.34e+00
1 SNES Function norm 5.67e-03
2 SNES Function norm 1.23e-08
```

`-ksp_monitor` prints the linear (KSP) residual at each conjugate gradient / GMRES iteration.
If the linear solver stagnates, the residual will not decrease.

#### View Solver Configuration

```ini
[Executioner]
  petsc_options = '-snes_view'
[]
```

This prints the complete PETSc solver configuration including preconditioner type, solver
hierarchy, and option settings. Use it to verify that the solver is configured as you intend.

#### Check Jacobian

```ini
[Executioner]
  petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
[]
```

This computes both the analytical Jacobian and a finite-difference approximation, prints the
comparison, and highlights large discrepancies. It runs for one Newton step and exits, so you
do not need a converged solve.

### TopResidualDebugOutput

MOOSE can print the variables and elements with the largest residual contributions:

```ini
[Outputs]
  [top_residual]
    type = TopResidualDebugOutput
    num_variables = 5     # show top 5 variables by residual contribution
    residual_type = 'NONLINEAR'
  []
[]
```

This is invaluable when the Newton solver fails to converge. The output identifies which physical
variable and which part of the domain is causing problems. You can then focus mesh refinement or
physics debugging on that region.

### GDB with MOOSE

MOOSE can be debugged with GDB like any other C++ program. The most useful workflow:

```bash
# Build with debug info
make -j8 METHOD=dbg

# Run under GDB
gdb --args ./my_app-dbg -i test.i

# In GDB:
(gdb) run                    # start execution
# When it crashes:
(gdb) bt                     # print stack trace (backtrace)
(gdb) frame 3                # move to frame 3 in the stack
(gdb) print _qp              # print member variable _qp
(gdb) print _u[_qp]          # print the solution value
```

To set a breakpoint in your kernel:

```bash
(gdb) break ReactionDiffusion::computeQpResidual
(gdb) run
# GDB stops when the breakpoint is hit
(gdb) print _qp
(gdb) print _grad_u[_qp]
(gdb) continue               # continue until next hit
```

#### Debugging Parallel Runs

When debugging an MPI run that crashes only on a specific rank:

```bash
# Launch with 4 processes, attach GDB to rank 0
mpiexec -n 4 xterm -e gdb ./my_app-dbg -i test.i
```

This opens 4 xterm windows, each running GDB for one MPI rank. Alternatively, add a sleep at
the start of `main.C` and attach GDB to the specific process ID:

```bash
# In another terminal:
ps aux | grep my_app-dbg
gdb -p <pid_of_rank_to_debug>
```

### Common Error Messages and What They Mean

| Error message | Likely cause |
|---------------|-------------|
| `Unable to find 'MyKernel' registering in...` | `registerMooseObject` not present in `.C` file, or wrong app name string |
| `Unknown parameter: 'my_typo'` | Input file parameter does not match any `addParam` declaration |
| `Material property 'conductivity' requested by 'MyKernel', but not defined` | No Material declares that property, or the property name strings do not match |
| `nonlinear iterations exceeded max` | Newton solver not converging. Check Jacobian, try smaller time step, add `-snes_monitor` |
| `Inf or NaN found in residual` | Division by zero or physically invalid state. Use `METHOD=dbg` and check field values |
| `Cannot convert type...` | Template type mismatch in `getMaterialProperty`, `getParam`, or similar |

---

## Summary: The Object Hierarchy at a Glance

```
MooseApp
  └─ FEProblemBase
       ├─ Kernels (volumetric PDEs)
       │    └─ Kernel → computeQpResidual() [and Jacobian]
       │    └─ ADKernel → computeQpResidual() returning ADReal
       ├─ BoundaryConditions
       │    └─ NodalBC → computeQpResidual() (at nodes)
       │    └─ IntegratedBC → computeQpResidual() (on boundary faces)
       ├─ Materials
       │    └─ Material → computeQpProperties() fills MaterialProperty arrays
       ├─ AuxKernels (derived quantities, not unknowns)
       ├─ Postprocessors → scalar values, written to CSV
       ├─ UserObjects → general-purpose data objects
       ├─ MeshGenerators → build the mesh before solve
       └─ MultiApps/Transfers → couple multiple simulations
```

Each object type has its own loop and its own method that you override. The pattern is always the
same: inherit, declare parameters, construct, override the compute method, register.

---

## Quick Reference Card

### Kernel Template

```cpp
// MyKernel.h
#pragma once
#include "Kernel.h"
class MyKernel : public Kernel {
public:
  static InputParameters validParams();
  MyKernel(const InputParameters & p);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  const Real _param;
};

// MyKernel.C
#include "MyKernel.h"
registerMooseObject("MyAppApp", MyKernel);
InputParameters MyKernel::validParams() {
  auto p = Kernel::validParams();
  p.addClassDescription("...");
  p.addRequiredParam<Real>("param", "...");
  return p;
}
MyKernel::MyKernel(const InputParameters & p)
  : Kernel(p), _param(getParam<Real>("param")) {}
Real MyKernel::computeQpResidual() { return /* ... */; }
Real MyKernel::computeQpJacobian() { return /* ... */; }
```

### Key Member Variables

```
_u[_qp]              solution at QP
_grad_u[_qp]         gradient of solution at QP
_test[_i][_qp]       test function i at QP
_grad_test[_i][_qp]  gradient of test function i at QP
_phi[_j][_qp]        shape function j at QP
_grad_phi[_j][_qp]   gradient of shape function j at QP
_q_point[_qp]        (x,y,z) coordinates of QP
_JxW[_qp]            quadrature weight * Jacobian determinant
_coord[_qp]          coordinate system scale factor
_qp                  current quadrature point index
_i                   current test function index (residual row)
_j                   current shape function index (Jacobian column)
```

### Build Commands

```bash
cd framework && make -j8 METHOD=opt      # build framework
cd test && make -j8 && ./run_tests -j8   # verify build
cd my_app && make -j8 METHOD=opt         # build your app
make compile_commands.json               # for editor support
```

### Test Commands

```bash
./run_tests -j8                          # run all tests
./run_tests --re=pattern                 # filter by name
./run_tests --dry-run --re=pattern       # show which would run
./run_tests --update-all --re=pattern    # regenerate gold files
METHOD=dbg ./run_tests --re=my_test      # run in debug mode
```
