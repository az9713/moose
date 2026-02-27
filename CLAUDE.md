# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is MOOSE?

MOOSE (Multiphysics Object-Oriented Simulation Environment) is a finite-element, multiphysics framework developed by Idaho National Laboratory. It builds on PETSc (nonlinear solvers) and libMesh (finite elements) to provide a high-level C++ API for writing physics simulation codes.

## Build System

MOOSE uses GNU Make with autotools configuration. Key environment variables:
- `MOOSE_DIR` — root of the repository (auto-detected)
- `LIBMESH_DIR` — libMesh installation (defaults to `$MOOSE_DIR/libmesh/installed`)
- `METHOD` — build type: `opt` (default), `dbg` (debug with assertions), `devel`
- `MOOSE_JOBS` — parallel make jobs (default: 8)

### Building

```bash
# Configure (optional, needed for libtorch/mfem/neml2/kokkos)
./configure [--with-derivative-size=64] [--with-libtorch] [--with-mfem] [--with-kokkos]

# Build the framework
cd framework && make -j8

# Build the test application
cd test && make -j8

# Build a specific module
cd modules/heat_transfer && make -j8

# Build all modules (combined app)
cd modules/combined && make -j8

# Generate compile_commands.json (uses METHOD=dbg)
make compile_commands.json
```

### Build Files

- `framework/build.mk` — libMesh integration, METHOD selection, compiler detection
- `framework/moose.mk` — core build logic for libmoose (PCRE, HIT, optional deps)
- `framework/app.mk` — template included by all application Makefiles
- `modules/modules.mk` — module dependency resolution and enablement flags

## Running Tests

Tests use the Python TestHarness. Every testable directory has a symlink `run_tests -> scripts/run_tests` and a `testroot` file.

```bash
# Run framework tests
cd test && ./run_tests -j8

# Run a specific test by name
cd test && ./run_tests -j8 --re=simple_diffusion

# Run tests for a single module
cd modules/solid_mechanics && ./run_tests -j8

# Run C++ unit tests (Google Test)
cd unit && make -j8 && ./run_tests

# Run Python tests
cd python && ./run_tests
```

Test specs are defined in `tests` files (HIT format) within test directories. Example:
```
[Tests]
  [test_name]
    type = 'Exodiff'
    input = 'input_file.i'
    exodiff = 'input_file_out.e'
    issues = '#1234'
    design = 'path/to/doc.md'
    requirement = 'Description of what the test verifies.'
  []
[]
```

Common test types: `Exodiff`, `CSVDiff`, `RunException`, `PythonUnitTest`.

## Architecture

### Hierarchy

```
Framework (core kernel)
  └── Modules (optional physics: heat_transfer, solid_mechanics, navier_stokes, ...)
       └── Applications (user codes built on framework + selected modules)
```

### Framework (`framework/`)

The core library (`libmoose`). Source in `framework/src/`, headers in `framework/include/`, both organized into ~80 subdirectories by subsystem. Key subsystems:

- `base/` — `MooseApp`, `MooseObject`, `Factory`, `AppFactory`, core infrastructure
- `kernels/` — weak form residual contributions (FEM)
- `fvkernels/` — finite volume kernels
- `bcs/`, `fvbcs/`, `linearfvbcs/` — boundary conditions
- `materials/`, `functormaterials/` — material property system
- `mesh/`, `meshgenerators/` — mesh handling and generation
- `executioners/` — solver strategies (Steady, Transient, Eigenvalue)
- `systems/` — NonlinearSystem, LinearSystem, AuxiliarySystem wrappers
- `multiapps/`, `transfers/` — multi-app coupling
- `outputs/` — Exodus, VTK, CSV, Console output
- `actions/` — automatic input-driven object creation
- `parser/` — HIT input file parsing
- `postprocessors/`, `reporters/` — scalar and general output quantities

### Modules (`modules/`)

25 optional physics modules. Each mirrors the framework directory structure (`include/`, `src/`, `test/`, `doc/`). Modules declare dependencies in `modules/modules.mk` — e.g., `FSI` requires `NAVIER_STOKES + SOLID_MECHANICS`, `THERMAL_HYDRAULICS` requires 7 other modules.

### MOOSE Object Pattern

Every MOOSE object (kernel, BC, material, etc.) follows this pattern:

**Header** (`.h`):
```cpp
#pragma once
#include "Kernel.h"

class MyKernel : public Kernel
{
public:
  static InputParameters validParams();
  MyKernel(const InputParameters & parameters);
protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};
```

**Source** (`.C`):
```cpp
#include "MyKernel.h"
registerMooseObject("MyApp", MyKernel);

InputParameters
MyKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("...");
  // add parameters here
  return params;
}

MyKernel::MyKernel(const InputParameters & parameters) : Kernel(parameters) {}
// implement computeQpResidual, computeQpJacobian
```

Key conventions: `registerMooseObject` macro registers with the Factory. `validParams()` is static and chains from parent. Return type goes on its own line above the function name.

### Input File Format (HIT)

MOOSE uses Hierarchical Input Text format (`.i` files):
```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]
[Variables]
  [u]
  []
[]
[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]
```

### Python Tooling (`python/`)

- `TestHarness/` — test execution framework (parallel scheduling, result diffing)
- `MooseDocs/` — documentation generation system
- `mooseutils/` — CSV/image diffing, postprocessor readers, git utilities
- `pyhit/` — HIT format parser (reads/writes `.i` files)
- `peacock/` — GUI application
- `chigger/` — visualization tools
- `moosesqa/` — Software Quality Assurance tools

### Creating a New Application

Use the stork generator:
```bash
cd ~/projects
$MOOSE_DIR/scripts/stork.sh MyApp
```

This creates a complete application skeleton with Makefile, src/include dirs, test harness setup, and documentation scaffolding.

## Documentation

Comprehensive documentation is available in `docs/`:

| Document | Description | Audience |
|----------|-------------|----------|
| [Architecture](docs/architecture.md) | System design with ASCII diagrams — class hierarchy, Factory, Action system, solve loop, MultiApp, parallel architecture | Developers |
| [Developer Guide](docs/developer-guide.md) | Step-by-step tutorial — writing kernels, materials, BCs, postprocessors, AD, testing, debugging | New C++ developers |
| [User Guide](docs/user-guide.md) | Complete input file reference — every block type, parameters, running simulations, postprocessing | Simulation users |
| [Quick Start](docs/quick-start.md) | 68 progressive worked examples — from 1D diffusion to chemical reactions and geochemistry, with complete input files | Everyone |
| [Zero to Hero](docs/zero-to-hero.md) | 8-phase study plan (~8 weeks) — prerequisites through production skills, using this repo | Self-learners |
| [Modules Reference](docs/modules-reference.md) | All 29 physics modules — capabilities, key classes, dependencies, coupling patterns | All users |
| [Docker Guide](docs/docker-guide.md) | Running MOOSE on Windows with Docker — installation, volume mounts, MPI, troubleshooting | Windows users |

## Code Style

- **C++**: `.clang-format` based on LLVM style — 2-space indent, 100-column limit, Allman braces, no tabs, return type on own line for definitions
- **Python**: 4-space indent. Black formatter and Ruff linter configured in `pyproject.toml`
- **General**: UTF-8, LF line endings, trim trailing whitespace, final newline (`.editorconfig`)
- **License header**: All source files start with the MOOSE LGPL 2.1 license comment block (`//* This file is part of the MOOSE framework ...`)
- **Source file extension**: C++ implementation files use `.C` (uppercase), headers use `.h`
