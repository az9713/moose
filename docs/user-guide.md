# MOOSE User Guide: Running Simulations with Input Files

This guide is for engineers and scientists who want to run simulations using a MOOSE-based
application. It assumes you have a compiled MOOSE application binary and some familiarity with
finite element concepts, but does not require C++ expertise to follow. After reading this guide
you should be able to write, run, and interpret results from MOOSE input files independently.

---

## Table of Contents

1. [What MOOSE Does for You](#1-what-moose-does-for-you)
2. [The HIT Input File Format](#2-the-hit-input-file-format)
3. [The Mesh Block](#3-the-mesh-block)
4. [The Variables Block](#4-the-variables-block)
5. [The Kernels Block](#5-the-kernels-block)
6. [The BCs Block](#6-the-bcs-block)
7. [The ICs Block](#7-the-ics-block)
8. [The Materials Block](#8-the-materials-block)
9. [The Functions Block](#9-the-functions-block)
10. [The Executioner Block](#10-the-executioner-block)
11. [The Postprocessors Block](#11-the-postprocessors-block)
12. [AuxVariables and AuxKernels](#12-auxvariables-and-auxkernels)
13. [The Outputs Block](#13-the-outputs-block)
14. [MultiApps and Transfers](#14-multiapps-and-transfers)
15. [Running Simulations](#15-running-simulations)
16. [Postprocessing Output](#16-postprocessing-output)
17. [Restarting Simulations](#17-restarting-simulations)

---

## 1. What MOOSE Does for You

### The Problem MOOSE Solves

MOOSE (Multiphysics Object-Oriented Simulation Environment) solves systems of partial differential
equations (PDEs) on unstructured meshes using the finite element method. At the heart of every
MOOSE simulation is a nonlinear system of equations assembled from three conceptual pieces:

```
R(u) = 0
```

where the residual `R(u)` is built by adding contributions from:

- **Kernels**: terms from the PDE integrated over the domain volume (e.g., diffusion, time
  derivative, reaction)
- **Boundary Conditions (BCs)**: constraints or fluxes imposed on domain boundaries
- **Initial Conditions (ICs)**: the starting state of each variable at time zero

As a concrete example, the transient diffusion equation

```
du/dt - nabla dot (nabla u) = f
```

maps onto MOOSE blocks as follows:

| PDE term           | MOOSE block      | Object type        |
|--------------------|------------------|--------------------|
| `du/dt`            | `[Kernels]`      | `TimeDerivative`   |
| `-nabla dot nabla u` | `[Kernels]`    | `Diffusion`        |
| `f` (source)       | `[Kernels]`      | `BodyForce`        |
| `u = g` on boundary | `[BCs]`         | `DirichletBC`      |
| `du/dn = h` on boundary | `[BCs]`    | `NeumannBC`        |
| `u(x, 0) = u_0`   | `[ICs]`          | `FunctionIC`       |

In the weak form, each kernel object implements one `computeQpResidual()` function that returns a
single quadrature-point contribution. MOOSE loops over all elements, all quadrature points, and all
kernels to assemble the global residual vector and Jacobian matrix. You do not write those loops.

### The Simulation Workflow

```
input.i  --->  [MOOSE application binary]  --->  output.e  (Exodus mesh + fields)
                                           --->  output.csv (postprocessors)
                                           --->  console output (convergence history)
```

The complete workflow is:

1. Write an HIT-format input file (`.i` extension by convention)
2. Run `./myapp-opt -i input.i`
3. MOOSE reads and validates the input file, then:
   - Builds the mesh
   - Creates and initializes all objects (kernels, BCs, materials, etc.)
   - Applies initial conditions
   - Calls the executioner's solve loop (steady or transient)
   - During/after the solve, calls postprocessors and writes outputs
4. Open `output.e` in ParaView or load `output.csv` in Python

### What You Control vs What MOOSE Handles

| You specify (input file)                  | MOOSE handles automatically                |
|-------------------------------------------|--------------------------------------------|
| Mesh geometry and topology                | Quadrature rule selection                  |
| Which variables to solve for              | Shape function evaluation                  |
| Which kernels/BCs/ICs apply               | Residual and Jacobian assembly             |
| Material property values or expressions  | Parallel domain decomposition              |
| Solver type and tolerances                | Newton iterations and linear solves        |
| Output format and frequency               | Time stepping logic (for Transient)        |
| Postprocessors to compute                 | Checkpoint/restart bookkeeping             |

---

## 2. The HIT Input File Format

MOOSE input files use the **Hierarchical Input Text** (HIT) format. It is a simple block-parameter
language that is human readable and easy to write by hand.

### Basic Syntax

**Blocks** are delimited by `[BlockName]` (open) and `[]` (close):

```
[BlockName]
  parameter_name = value
[]
```

**Nested sub-blocks** follow the same pattern inside an outer block:

```
[OuterBlock]
  [sub_one]
    type = SomeType
    variable = u
  []
  [sub_two]
    type = AnotherType
    variable = u
    value = 3.14
  []
[]
```

The sub-block name (`sub_one`, `sub_two`) is the object's name within MOOSE. It appears in log
output and can be referenced by other objects. The `type` parameter selects which registered C++
class to instantiate.

**Older syntax** (still valid in existing files): sub-blocks used `[./name]` / `[../]` delimiters.
Both styles are accepted:

```
[Kernels]
  [./diff]          # older style
    type = Diffusion
    variable = u
  [../]

  [new_style_diff]  # current style
    type = Diffusion
    variable = u
  []
[]
```

### Parameter Types

| HIT type          | Example                                   | Notes                                      |
|-------------------|-------------------------------------------|--------------------------------------------|
| `Real`            | `value = 3.14`                            | Floating-point scalar                      |
| `Integer`         | `nx = 10`                                 | Integer scalar                             |
| `string`          | `boundary = left`                         | Unquoted or single-quoted                  |
| `boolean`         | `exodus = true`                           | `true`/`false` (case-insensitive)          |
| `vector of Real`  | `prop_values = '1.0 2.5 0.0'`            | Space-separated, enclosed in single quotes |
| `vector of string`| `prop_names = 'k rho cp'`                | Space-separated, enclosed in single quotes |
| `MooseEnum`       | `solve_type = PJFNK`                      | One value from a predefined list           |
| `FunctionName`    | `function = my_func`                      | Name of an object in `[Functions]`         |
| `PostprocessorName`| `postprocessor = avg_temp`               | Name of an object in `[Postprocessors]`    |

Vectors are always written inside single quotes with space separation:

```
[Materials]
  [props]
    type = GenericConstantMaterial
    prop_names  = 'thermal_conductivity density specific_heat'
    prop_values = '18.0 8000.0 0.466'
  []
[]
```

### Comments

The `#` character starts a line comment. Everything from `#` to the end of the line is ignored:

```
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100   # 100 elements in x gives adequate resolution
    ny = 10
    xmax = 0.304   # length of chamber (m)
  []
[]
```

### GlobalParams

The `[GlobalParams]` block sets default parameter values that apply to all objects in the file.
This is useful when many objects share the same parameter, such as `variable`:

```
[GlobalParams]
  variable = u
[]

[Kernels]
  [diff]
    type = Diffusion
    # variable = u  <-- inherited from GlobalParams, no need to repeat
  []
  [td]
    type = TimeDerivative
    # variable = u  <-- also inherited
  []
[]
```

If an object specifies the parameter explicitly, that value overrides the global default. From the
test file `test/tests/globalparams/global_param/global_param_test.i`:

```
[GlobalParams]
  variable = u
  dim = 2
[]
```

### Variable Substitution with `${var}`

You can define variables at the top of an input file and substitute them later using `${varname}`:

```
nx = 50
dt = 0.01

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${nx}
  []
[]

[Executioner]
  type = Transient
  dt = ${dt}
[]
```

This avoids duplicating magic numbers and makes parametric studies easier. You can also define a
variable with an empty default and set it from the command line (see below).

A real example from `test/tests/kernels/2d_diffusion/bodyforce.i`:

```
AD = ''

[Kernels]
  [diff]
    type = ${AD}Diffusion   # becomes 'Diffusion' since AD is empty
    variable = u
  []
[]
```

Setting `AD = AD` on the command line would select `ADDiffusion` instead.

### Command-Line Overrides

Any parameter in the input file can be overridden from the command line by specifying the full
path `BlockName/SubBlockName/parameter=value`:

```bash
# Override a mesh parameter
./myapp-opt -i input.i Mesh/gmg/nx=200

# Override a BC value
./myapp-opt -i input.i BCs/left/value=350

# Override executioner settings
./myapp-opt -i input.i Executioner/num_steps=100 Executioner/dt=0.05

# Override a top-level variable
./myapp-opt -i input.i nx=200
```

This is extremely useful for parameter sweeps in scripts without editing the input file.

### Multi-File Input

You can split an input across multiple files. MOOSE merges them in order, with later files
overriding earlier ones. This is useful for having a base configuration and a separate override:

```bash
./myapp-opt -i base.i override.i
```

Parameters in `override.i` take precedence over `base.i` for any conflicting keys.

---

## 3. The Mesh Block

The `[Mesh]` block defines the computational domain. MOOSE uses a **mesh generator pipeline**:
one or more sub-blocks, each of type `MeshGenerator`, that are chained together. The last generator
in the chain (or the one declared final) produces the mesh used by the simulation.

### GeneratedMeshGenerator

For simple rectangular/box domains, use `GeneratedMeshGenerator`. This is the most common starting
point for test problems:

```
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2          # 1, 2, or 3
    nx = 50          # elements in x
    ny = 50          # elements in y
    nz = 1           # elements in z (3D only)
    xmin = 0.0       # domain bounds (all default to 0 or 1)
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
  []
[]
```

For a 3D box:

```
[Mesh]
  [box]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 0.5
  []
[]
```

`GeneratedMeshGenerator` automatically creates named sidesets:

| `dim = 1` | `dim = 2`               | `dim = 3`                                      |
|-----------|-------------------------|------------------------------------------------|
| `left`    | `left`, `right`         | `left`, `right`                                |
| `right`   | `bottom`, `top`         | `bottom`, `top`                                |
|           |                         | `front`, `back`                                |

These names can be used directly in `[BCs]` without any additional setup.

### FileMeshGenerator

To read an existing mesh file (Exodus `.e`/`.exd`, Gmsh `.msh`, or other formats supported by
libMesh):

```
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = my_geometry.e
  []
[]
```

Exodus files are the most common format in the MOOSE ecosystem, as they support multiple time
steps, element blocks, and sidesets. You can export Exodus meshes from Cubit/Coreform Cubit, or
use MOOSE's own `--mesh-only` flag (see below) to produce an Exodus mesh from a generator pipeline.

For Gmsh `.msh` files:

```
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = mesh.msh
  []
[]
```

libMesh will automatically detect the format from the file extension.

### Mesh Generator Pipeline

You can chain generators to build complex meshes step by step. Each generator takes an `input`
parameter naming the generator whose output it should receive:

```
[Mesh]
  # Step 1: generate a base rectangle
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
    xmax = 2.0
  []

  # Step 2: rename a subdomain for block-restricted physics
  [fuel_region]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = fuel
    bottom_left = '0.0 0.0 0.0'
    top_right   = '1.0 1.0 0.0'
  []

  # Step 3: add a named sideset between the two subdomains
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = fuel_region
    primary_block = fuel
    paired_block = 0
    new_boundary = fuel_interface
  []
[]
```

The pipeline is evaluated in dependency order (not necessarily text order). MOOSE determines the
final mesh by identifying the generator that no other generator depends on.

### Subdomains (Blocks) and Sidesets (Boundaries)

MOOSE uses two types of geometric labels:

- **Subdomains** (also called "blocks"): volumetric regions identified by integer IDs or names.
  Used to restrict kernels and materials to a portion of the mesh.
- **Sidesets** (also called "boundaries"): surfaces (in 3D) or edges (in 2D) identified by integer
  IDs or names. Used to specify boundary conditions.

You reference these by name or ID in `[BCs]`, `[Kernels]`, and `[Materials]`:

```
[BCs]
  [wall_temp]
    type = DirichletBC
    variable = T
    boundary = 'left right top'   # multiple sidesets as a space-separated vector
    value = 300
  []
[]

[Materials]
  [fuel_props]
    type = GenericConstantMaterial
    block = fuel                   # restrict to the 'fuel' subdomain
    prop_names  = 'k rho cp'
    prop_values = '3.0 10960.0 247.0'
  []
[]
```

### Checking the Mesh with --mesh-only

To generate and inspect the mesh without running a solve, use the `--mesh-only` flag:

```bash
./myapp-opt -i input.i --mesh-only
```

This writes the final mesh to `input_mesh.e` (or a name derived from `file_base` if set). Open
this in ParaView to verify that sidesets and subdomains are named correctly before investing time
in a full solve.

---

## 4. The Variables Block

The `[Variables]` block declares the primary unknowns that MOOSE solves for. These are the fields
that appear in the nonlinear system `R(u) = 0`.

### Field Variables (Continuous Galerkin)

The default variable is a continuous Lagrange finite element field. The minimum declaration is:

```
[Variables]
  [u]
  []
[]
```

This creates a first-order LAGRANGE scalar variable named `u`. You can specify the polynomial
family and order explicitly:

```
[Variables]
  [temperature]
    family = LAGRANGE
    order  = FIRST       # linear elements
  []
  [pressure]
    family = LAGRANGE
    order  = SECOND      # quadratic elements, more accurate but more DOFs
  []
[]
```

**Available families and orders** (common ones):

| Family      | Orders available | Notes                                          |
|-------------|-----------------|------------------------------------------------|
| `LAGRANGE`  | `FIRST`, `SECOND` | Most common; continuous across element edges  |
| `HERMITE`   | `THIRD`         | C1-continuous; useful for beam/plate problems  |
| `MONOMIAL`  | `CONSTANT`, `FIRST`, `SECOND` | Discontinuous; used for DG and FV methods |
| `L2_LAGRANGE` | `FIRST`, `SECOND` | Discontinuous Lagrange (L2 projection)      |

For most heat conduction, diffusion, and structural mechanics problems, `LAGRANGE FIRST` is
appropriate. Use `LAGRANGE SECOND` when you need higher accuracy or when solving problems with
curved geometries.

### Finite Volume Variables

For finite volume discretizations, set `fv = true`. The variable stores one value per element
centroid (equivalent to `MONOMIAL CONSTANT`):

```
[Variables]
  [pressure]
    type = MooseVariableFVReal
    # or equivalently:
    fv = true
  []
[]
```

Finite volume variables are used with FV kernels (in `[FVKernels]`) and FV boundary conditions
(in `[FVBCs]`), which are separate blocks from the standard `[Kernels]` and `[BCs]`.

### Scaling

Variable scaling improves solver conditioning when variables have very different magnitudes (e.g.,
temperature in Kelvin alongside pressure in Pascals):

```
[Variables]
  [temperature]
    scaling = 1e-3    # divide internal residual by this factor
  []
  [pressure]
    scaling = 1e-5
  []
[]
```

Scaling can also be set automatically with `automatic_scaling = true` in `[Problem]`.

### Inline Initial Conditions

A simple initial condition can be set directly inside the variable block:

```
[Variables]
  [temperature]
    initial_condition = 300   # uniform initial value of 300 K
  []
[]
```

This is equivalent to a `ConstantIC` in the `[ICs]` block. For spatially varying initial
conditions, use the `[ICs]` block instead (see Section 7).

---

## 5. The Kernels Block

Each sub-block in `[Kernels]` represents **one term in the weak form** of the PDE. MOOSE
assembles the complete residual by summing the contributions of all active kernels.

### What a Kernel Represents

In the finite element method, a PDE is multiplied by a test function `psi_i` and integrated over
the domain. Each term in that integral becomes a kernel. For example, the steady Poisson equation

```
-nabla dot (D nabla u) = f
```

in weak form (after integration by parts) is

```
integral[ D (nabla u) dot (nabla psi_i) ] = integral[ f psi_i ]
```

The left side is the `MatDiffusion` kernel and the right side (moved to the residual) is the
`BodyForce` kernel with a negative sign. You combine them in `[Kernels]`:

```
[Kernels]
  [diffusion]
    type = MatDiffusion
    variable = u
    diffusivity = D    # name of a material property
  []
  [source]
    type = BodyForce
    variable = u
    value = -1.0       # constant forcing (positive body force adds to residual as -f)
  []
[]
```

### Built-In Kernels Reference

The following kernels are registered in `framework/src/kernels/` and are available in any MOOSE
application.

#### Diffusion Kernels

**`Diffusion`** — The Laplacian operator with unit diffusivity.

- Weak form: `(nabla psi_i, nabla u_h)`
- Represents: `-nabla^2 u = 0` (or more precisely, contributes `nabla^2 u` term to residual)
- Source: `framework/src/kernels/Diffusion.C`

```
[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]
```

**`ADDiffusion`** — Automatic differentiation version of `Diffusion`. Computes the Jacobian
automatically. Preferred in new code for correctness and convenience.

```
[Kernels]
  [diff]
    type = ADDiffusion
    variable = u
  []
[]
```

**`MatDiffusion`** — Diffusion with a spatially/temporally varying diffusivity taken from a
material property.

- Weak form: `(nabla psi_i, D nabla u_h)` where `D` is a material property
- Represents: `-nabla dot (D nabla u)`
- Source: `framework/src/kernels/MatDiffusion.C`

```
[Kernels]
  [diff]
    type = MatDiffusion
    variable = u
    diffusivity = thermal_conductivity   # must match a prop_name in [Materials]
  []
[]
```

Also available as `ADMatDiffusion` for automatic differentiation.

#### Time Derivative Kernels

**`TimeDerivative`** — First-order time derivative (for transient problems).

- Weak form: `(psi_i, du_h/dt)`
- Represents: `du/dt` term in a transient PDE
- Source: `framework/src/kernels/TimeDerivative.C`

```
[Kernels]
  [td]
    type = TimeDerivative
    variable = u
  []
[]
```

**`ADTimeDerivative`** — AD version of `TimeDerivative`. Use in combination with AD kernels.

```
[Kernels]
  [td]
    type = ADTimeDerivative
    variable = u
  []
[]
```

Both versions support an optional `lumping = true` parameter that applies mass lumping for
stabilization in explicit time integration.

#### Reaction and Source Kernels

**`Reaction`** — A linear consuming reaction term.

- Weak form: `(psi_i, lambda * u_h)`
- Represents: `lambda * u` added to the residual (a first-order decay/sink)
- Source: `framework/src/kernels/Reaction.C`

```
[Kernels]
  [rxn]
    type = Reaction
    variable = u
    rate = 1.0    # lambda; defaults to 1.0
  []
[]
```

`rate` is controllable (can be varied during a simulation via a `Control` object).

**`BodyForce`** — A volumetric source/sink term that can be a constant, a function, or a
postprocessor value.

- Weak form: `(psi_i, -f)` where `f = value * postprocessor * function(t, x)`
- Source: `framework/src/kernels/BodyForce.C`

```
[Kernels]
  [src]
    type = BodyForce
    variable = u
    value = 10.0              # constant multiplier (default 1.0)
    function = heat_source    # optional: time/space-varying function
  []
[]
```

```
[Kernels]
  [ramp_source]
    type = BodyForce
    variable = u
    postprocessor = power_pp  # value provided by a postprocessor each step
  []
[]
```

**`CoupledForce`** — A source term proportional to the value of another variable.

- Weak form: `(psi_i, -sigma * v)` where `v` is a coupled variable
- Source: `framework/src/kernels/CoupledForce.C`

```
[Kernels]
  [heat_source_from_neutronics]
    type = CoupledForce
    variable = temperature
    v    = neutron_flux    # the variable providing the source
    coef = 2.5e8           # scaling coefficient sigma
  []
[]
```

#### Advection Kernels

**`ConservativeAdvection`** — Conservative divergence form of advection.

- Weak form: `(-nabla psi_i, vec_v * u)`
- Represents: `nabla dot (vec_v u)` in conservative form
- Source: `framework/src/kernels/ConservativeAdvection.C`

The velocity can be supplied as a coupled variable, a vector variable, or a material property:

```
[Kernels]
  [advection]
    type = ConservativeAdvection
    variable = c
    velocity_variable = vel     # coupled VECTOR variable for velocity
    upwinding_type = none       # or 'full' for full upwinding
  []
[]
```

### Combining Kernels: A Complete Example

The transient advection-diffusion-reaction equation

```
dc/dt + nabla dot (v c) - nabla dot (D nabla c) + lambda c = f
```

is assembled by listing all contributing kernels:

```
[Kernels]
  [time_deriv]
    type = TimeDerivative
    variable = c
  []
  [advect]
    type = ConservativeAdvection
    variable = c
    velocity_variable = velocity
    upwinding_type = full
  []
  [diffuse]
    type = MatDiffusion
    variable = c
    diffusivity = diffusivity
  []
  [decay]
    type = Reaction
    variable = c
    rate = 0.01
  []
  [source]
    type = BodyForce
    variable = c
    function = source_func
  []
[]
```

MOOSE sums all kernel residuals into a single system. Each object is independent and can be
activated or deactivated by adding/removing it from the block.

### Block-Restricting Kernels

A kernel can be restricted to a specific subdomain with the `block` parameter:

```
[Kernels]
  [diff_fuel]
    type = MatDiffusion
    variable = temperature
    diffusivity = k_fuel
    block = fuel_block
  []
  [diff_coolant]
    type = MatDiffusion
    variable = temperature
    diffusivity = k_coolant
    block = coolant_block
  []
[]
```

---

## 6. The BCs Block

The `[BCs]` block specifies boundary conditions. Each sub-block applies one BC to one variable on
one or more sidesets.

### Dirichlet (Essential) Boundary Conditions

Dirichlet BCs fix the value of the variable on a boundary: `u = g`.

**`DirichletBC`** — Constant value Dirichlet BC.

- Enforces: `u = value` on the specified boundary
- Source: `framework/src/bcs/DirichletBC.C`

```
[BCs]
  [inlet]
    type = DirichletBC
    variable = pressure
    boundary = left
    value = 4000          # in Pa
  []
  [outlet]
    type = DirichletBC
    variable = pressure
    boundary = right
    value = 0
  []
[]
```

**`FunctionDirichletBC`** — Time and space-varying Dirichlet BC using a Function object.

- Enforces: `u = g(t, x)` on the boundary
- Source: `framework/src/bcs/FunctionDirichletBC.C`

```
[BCs]
  [heated_wall]
    type = FunctionDirichletBC
    variable = temperature
    boundary = left
    function = ramp_temperature   # name of a [Functions] sub-block
  []
[]
```

**`ADDirichletBC`** — AD version of `DirichletBC`. Use when all kernels use AD.

```
[BCs]
  [wall]
    type = ADDirichletBC
    variable = u
    boundary = 'left right'
    value = 0
  []
[]
```

You can apply the same BC to multiple sidesets by listing them in a single `boundary` vector:

```
[BCs]
  [zero_walls]
    type = DirichletBC
    variable = u
    boundary = 'left right top bottom'
    value = 0
  []
[]
```

### Neumann (Natural) Boundary Conditions

Neumann BCs specify the normal derivative (or normal flux) on a boundary:
`D * du/dn = h`. In the weak form, this appears as a boundary integral term.

**`NeumannBC`** — Constant flux BC.

- Enforces: `du/dn = value` (the constant prescribed gradient dotted with the outward normal)
- Source: `framework/src/bcs/NeumannBC.C`

```
[BCs]
  [heat_flux_in]
    type = NeumannBC
    variable = temperature
    boundary = left
    value = 1000    # W/m^2 (for a heat conduction problem)
  []
[]
```

A zero Neumann BC (no flux, insulating wall) is the **natural** BC in Galerkin FEM — you achieve
it simply by not listing any BC on that boundary.

**`FunctionNeumannBC`** — Time and space-varying flux BC.

- Enforces: `du/dn = h(t, x)` on the boundary
- Source: `framework/src/bcs/FunctionNeumannBC.C`

```
[BCs]
  [varying_flux]
    type = FunctionNeumannBC
    variable = temperature
    boundary = top
    function = flux_profile   # a Function object
  []
[]
```

### Periodic Boundary Conditions

Periodic BCs are declared inside a special `[Periodic]` sub-block of `[BCs]`. They link
two opposing boundaries so that the solution and its gradient are continuous across the
periodic interface. This is common for crystal plasticity, pattern formation, and unit cell
simulations.

From `test/tests/bcs/periodic/periodic_bc_test.i`:

```
[BCs]
  [Periodic]
    [x_direction]
      variable = u
      primary   = 3        # left boundary (sideset ID)
      secondary = 1        # right boundary (sideset ID)
      translation = '40 0 0'   # vector from primary to secondary
    []
    [y_direction]
      variable = u
      primary   = 0
      secondary = 2
      translation = '0 40 0'
    []
  []
[]
```

For `GeneratedMeshGenerator` meshes you can use sideset names (`left`, `right`, etc.) instead
of integer IDs.

### How Sidesets Map to BCs

Sidesets are named groups of element faces (in 3D) or edges (in 2D). Every BC object has a
`boundary` parameter that accepts one or more sideset names or IDs:

```
[BCs]
  [no_slip_walls]
    type = DirichletBC
    variable = velocity
    boundary = 'top bottom'   # applies to both top and bottom walls
    value = 0
  []
[]
```

Use `--mesh-only` and open the result in ParaView (filter: "Extract Block") to confirm that
your sidesets are correctly defined before writing BCs.

---

## 7. The ICs Block

The `[ICs]` block specifies how field variables are initialized at the start of a simulation
(time = 0). ICs are required for transient problems and optional (but sometimes useful) for steady
problems where a good initial guess aids convergence.

### ConstantIC

Sets a uniform initial value across the entire domain or a subdomain.

- Source: `framework/src/ics/ConstantIC.C`

```
[ICs]
  [temp_ic]
    type = ConstantIC
    variable = temperature
    value = 300     # all nodes start at 300 K
  []
[]
```

This is equivalent to `initial_condition = 300` directly in the `[Variables]` block.

### FunctionIC

Sets the initial value from a Function object, enabling spatially and/or temporally varying
initial conditions.

- Source: `framework/src/ics/FunctionIC.C`

```
[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = initial_profile   # defined in [Functions]
  []
[]

[Functions]
  [initial_profile]
    type = ParsedFunction
    expression = 'exp(-((x-0.5)^2 + (y-0.5)^2) / 0.01)'
  []
[]
```

From `test/tests/ics/function_ic/parsed_function.i`:

```
[ICs]
  [u_ic]
    type = FunctionIC
    variable = 'u'
    function = parsed_function
  []
[]
```

### BoundingBoxIC

Sets one value inside a rectangular bounding box and another value outside. Useful for
setting up phase-field simulations, diffuse interface problems, or patch tests.

- Source: `framework/src/ics/BoundingBoxIC.C`

From `test/tests/ics/bounding_box_ic/bounding_box_ic_test.i`:

```
[ICs]
  [u_ic]
    type = BoundingBoxIC
    variable = u
    x1 = 0.1   # lower-left corner of the box
    y1 = 0.1
    x2 = 0.6   # upper-right corner
    y2 = 0.6
    inside  = 2.3   # value inside the box
    outside = 4.6   # value outside the box
  []
[]
```

For 3D, also specify `z1` and `z2`.

### Inline IC Syntax

ICs can also be declared inline within the `[Variables]` block:

```
[Variables]
  [u]
    [InitialCondition]
      type = BoundingBoxIC
      x1 = 0.1
      y1 = 0.1
      x2 = 0.6
      y2 = 0.6
      inside  = 2.3
      outside = 4.6
    []
  []
[]
```

The inline syntax is older (from the `[./name][../]` era) but still valid. The separate `[ICs]`
block is generally clearer for complex setups.

---

## 8. The Materials Block

The `[Materials]` block defines material properties — spatially varying scalar, vector, or tensor
fields that kernels, BCs, and auxiliary kernels can query by name. A material property is computed
at every quadrature point on every element during each residual evaluation.

### GenericConstantMaterial

Creates one or more scalar material properties with constant values. The most common material
object for simple problems.

- Source: `framework/src/materials/GenericConstantMaterial.C`

```
[Materials]
  [steel]
    type = GenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '18.0 0.466 8000.0'    # W/(m K), J/(kg K), kg/m^3
  []
[]
```

The AD version `ADGenericConstantMaterial` is used when kernels use automatic differentiation:

```
[Materials]
  [steel]
    type = ADGenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '18.0 0.466 8000.0'
  []
[]
```

Both versions take the same parameters. The names in `prop_names` must match exactly the names
requested by kernels (e.g., `diffusivity = thermal_conductivity`).

### GenericFunctionMaterial

Creates material properties whose values are evaluated from Function objects. This allows
time-varying or spatially-varying material properties without writing C++ code.

- Source: `framework/src/materials/GenericFunctionMaterial.C`

```
[Materials]
  [conductivity]
    type = GenericFunctionMaterial
    prop_names  = 'k'
    prop_values = 'k_func'   # name(s) of [Functions] objects
  []
[]

[Functions]
  [k_func]
    type = ParsedFunction
    expression = '10.0 + 0.01 * t'   # linearly increasing conductivity
  []
[]
```

### ParsedMaterial

Computes a material property from a math expression that can reference:
- The current coordinates (`x`, `y`, `z`, `t`)
- Coupled variable values
- Other material properties (declared first)

- Source: `framework/src/materials/ParsedMaterial.C`

```
[Materials]
  [free_energy]
    type = ParsedMaterial
    expression = '(eta - 0.5)^2'
    coupled_variables = 'eta'
    outputs = exodus    # optional: output this property to the Exodus file
  []
[]
```

The `coupled_variables` list tells MOOSE which FE variables to expose in the expression.
Available variables in the expression are: `x`, `y`, `z`, `t`, any names listed in
`coupled_variables`, and any declared constants.

### DerivativeParsedMaterial

Like `ParsedMaterial` but also automatically computes analytic derivatives of the property
with respect to each coupled variable. Required for use with kernels that need material
property derivatives (e.g., phase-field Allen-Cahn/Cahn-Hilliard kernels).

- Source: `framework/src/materials/DerivativeParsedMaterial.C`

```
[Materials]
  [f_bulk]
    type = DerivativeParsedMaterial
    expression = 'W * eta^2 * (1-eta)^2'
    coupled_variables = 'eta'
    constant_names  = 'W'
    constant_expressions = '1.0'
  []
[]
```

MOOSE will compute `df/d_eta`, `d2f/d_eta2`, etc. automatically and make them available as
material properties named `df_bulk/deta`, `d2f_bulk/deta2`, etc.

### Block Restriction

Any material object can be restricted to specific subdomains using the `block` parameter:

```
[Materials]
  [fuel_props]
    type = GenericConstantMaterial
    block = 'fuel_block'     # only computed on elements in this subdomain
    prop_names  = 'k rho cp'
    prop_values = '3.5 10500 250'
  []
  [cladding_props]
    type = GenericConstantMaterial
    block = 'cladding_block'
    prop_names  = 'k rho cp'
    prop_values = '16.0 6550 330'
  []
[]
```

This is essential in multimaterial simulations where different regions have different physical
properties. The kernels using those properties should be restricted to the same blocks.

---

## 9. The Functions Block

The `[Functions]` block defines mathematical functions of space and time that can be referenced
by BCs, ICs, materials, kernels, and postprocessors. Functions are evaluated on demand rather
than stored at quadrature points.

### ParsedFunction

Evaluates an arbitrary mathematical expression. The expression can use coordinates (`x`, `y`,
`z`, `t`), standard math functions (`sin`, `cos`, `exp`, `log`, `sqrt`, etc.), and user-defined
variables and values.

- Registered as: `ParsedFunction` (alias for `MooseParsedFunction`)
- Source: `framework/src/functions/MooseParsedFunction.C`

```
[Functions]
  [ramp]
    type = ParsedFunction
    expression = 't'           # linearly increasing with time
  []

  [gaussian_source]
    type = ParsedFunction
    expression = 'exp(-((x-0.5)^2 + (y-0.5)^2) / sigma^2)'
    symbol_names  = 'sigma'
    symbol_values = '0.1'
  []

  [sinusoidal_bc]
    type = ParsedFunction
    expression = 'amplitude * sin(2 * pi * f * t)'
    symbol_names  = 'amplitude f'
    symbol_values = '100.0 0.5'
  []
[]
```

You can also reference other postprocessors or scalar variables in the expression via `vars`
and `vals` parameters (for coupling postprocessor values into function expressions).

### PiecewiseLinear

Performs linear interpolation between a set of (x, y) data pairs. The independent variable
is usually time (the default), but can be `x`, `y`, or `z`.

- Source: `framework/src/functions/PiecewiseLinear.C`

```
[Functions]
  [power_history]
    type = PiecewiseLinear
    x = '0  100  200  300  500'     # time values (seconds)
    y = '0  1e6  1e6  5e5  5e5'    # power values (W/m^3)
  []
[]
```

The `axis` parameter selects the independent variable:

```
[Functions]
  [wall_temperature_profile]
    type = PiecewiseLinear
    axis = x
    x = '0.0  0.1  0.5  1.0'
    y = '300  350  400  300'
  []
[]
```

Set `extrap = true` to extrapolate outside the data range (default is to clamp at the
endpoint values).

### ConstantFunction

Returns the same value everywhere and at all times. Mostly used as a placeholder or for
testing, since you could also use a constant directly in the parameter.

- Source: `framework/src/functions/ConstantFunction.C`

```
[Functions]
  [zero_flux]
    type = ConstantFunction
    value = 0.0
  []
[]
```

### Using Functions in BCs, ICs, and Materials

Functions are referenced by name wherever a `FunctionName` parameter is accepted:

```
[BCs]
  [time_varying_left]
    type = FunctionDirichletBC
    variable = u
    boundary = left
    function = ramp           # name of a [Functions] sub-block
  []
[]

[ICs]
  [u_ic]
    type = FunctionIC
    variable = u
    function = gaussian_source
  []
[]

[Materials]
  [k_mat]
    type = GenericFunctionMaterial
    prop_names  = 'k'
    prop_values = 'power_history'
  []
[]
```

---

## 10. The Executioner Block

The `[Executioner]` block controls the solver: how the nonlinear system is solved (steady or
transient), which linear and nonlinear solver algorithms are used, and what convergence tolerances
apply.

### Steady Executioner

Use `type = Steady` for problems with no time derivative. MOOSE performs one nonlinear solve
and writes final output.

```
[Executioner]
  type = Steady
  solve_type = NEWTON      # nonlinear solver algorithm
  nl_rel_tol = 1e-8        # nonlinear relative residual tolerance
  nl_abs_tol = 1e-10       # nonlinear absolute residual tolerance
  l_tol      = 1e-5        # linear solver relative tolerance

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
```

From `tutorials/darcy_thermo_mech/step01_diffusion/problems/step1.i`:

```
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
```

### Transient Executioner

Use `type = Transient` for time-dependent problems. You must have at least one `TimeDerivative`
(or equivalent) kernel in `[Kernels]`.

```
[Executioner]
  type = Transient
  num_steps = 100         # total number of time steps
  dt        = 0.01        # time step size (if using ConstantDT, the default)
  end_time  = 10.0        # stop after this time (alternative to num_steps)

  solve_type = PJFNK      # Preconditioned Jacobian-Free Newton-Krylov
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_tol      = 1e-5
  l_max_its  = 300        # maximum linear iterations per Newton step

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
```

You can control when the transient ends with either `num_steps` or `end_time` (or both; whichever
is reached first stops the simulation).

### Solve Types

| `solve_type` | Description                                                  | When to use              |
|--------------|--------------------------------------------------------------|--------------------------|
| `NEWTON`     | Full Newton with exact Jacobian assembled from hand-coded `computeQpJacobian` or AD | AD-based problems; ensures quadratic convergence |
| `PJFNK`      | Preconditioned Jacobian-Free Newton-Krylov: Jacobian-vector products via finite differences of residual | When exact Jacobian is unavailable or expensive; most common |
| `JFNK`       | Jacobian-Free Newton-Krylov (no preconditioner)              | Rarely preferred; PJFNK is better |
| `LINEAR`     | Single linear solve (for linear problems)                    | Linear elliptic problems  |

For new code using AD (automatic differentiation) objects, `NEWTON` is preferred because AD
provides exact Jacobians at low developer cost.

### PETSc Options

MOOSE uses PETSc as its underlying linear algebra library. Preconditioner and solver options
are passed through `petsc_options_iname` / `petsc_options_value` pairs:

```
[Executioner]
  type = Steady
  solve_type = PJFNK

  # Hypre BoomerAMG algebraic multigrid (excellent for elliptic problems)
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
```

```
[Executioner]
  type = Steady
  solve_type = PJFNK

  # ILU with fill level 2 (for advection-diffusion)
  petsc_options_iname = '-pc_type -pc_factor_levels'
  petsc_options_value = 'ilu 2'
[]
```

```
[Executioner]
  type = Steady
  solve_type = NEWTON

  # LU direct solve (for small problems or debugging)
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
```

### Common Solver Configurations

| Problem type                    | `solve_type` | Preconditioner           |
|---------------------------------|-------------|--------------------------|
| Steady diffusion / Laplacian    | `PJFNK`     | `hypre boomeramg`        |
| Steady diffusion with AD        | `NEWTON`    | `hypre boomeramg`        |
| Transient heat conduction       | `PJFNK`     | `hypre boomeramg`        |
| Transient heat conduction (AD)  | `NEWTON`    | `hypre boomeramg`        |
| Advection dominated             | `PJFNK`     | `ilu` with level 1-4     |
| Structural mechanics            | `PJFNK`     | `hypre boomeramg`        |
| Small debugging problem         | `NEWTON`    | `lu` (direct)            |

### Convergence Criteria

The nonlinear solver iterates until the residual satisfies:

```
||R(u^n)|| < nl_abs_tol                             (absolute tolerance)
||R(u^n)|| / ||R(u^0)|| < nl_rel_tol               (relative tolerance, reduction from initial)
```

Both criteria are checked; convergence is declared when either is met.

| Parameter    | Default  | Meaning                                           |
|--------------|----------|---------------------------------------------------|
| `nl_rel_tol` | `1e-8`   | Nonlinear relative residual tolerance             |
| `nl_abs_tol` | `1e-50`  | Nonlinear absolute residual tolerance             |
| `nl_max_its` | `50`     | Maximum nonlinear iterations per time step        |
| `l_tol`      | `1e-5`   | Linear (Krylov) solver relative tolerance         |
| `l_max_its`  | `10000`  | Maximum linear iterations                         |
| `l_abs_tol`  | `1e-50`  | Linear absolute tolerance                         |

A typical starting point for production simulations is `nl_rel_tol = 1e-8` and
`nl_abs_tol = 1e-10`. Tighten if your postprocessors show sensitivity to solver tolerance.

### Time Steppers

For transient simulations you can specify an adaptive time stepper inside the executioner.
The default is a constant time step (`ConstantDT`).

**`IterationAdaptiveDT`** — Adjusts `dt` based on the number of nonlinear iterations in the
previous step. Grows the time step when convergence is fast and cuts it when convergence is slow:

```
[Executioner]
  type = Transient
  end_time = 100.0
  nl_abs_tol = 1e-14

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 6       # target this many NL iterations per step
    growth_factor = 1.5          # multiply dt by this if fewer than optimal
    cutback_factor = 0.5         # multiply dt by this if more than optimal
    iteration_window = 2         # tolerance around optimal_iterations
  []
[]
```

**`FunctionDT`** — Time step size determined by a user-defined Function:

```
[Executioner]
  type = Transient
  end_time = 1000

  [TimeStepper]
    type = FunctionDT
    function = dt_func
  []
[]

[Functions]
  [dt_func]
    type = PiecewiseLinear
    x = '0   10  100  1000'
    y = '0.1  1   10   100'
  []
[]
```

### Time Integrators

The time integrator is selected with `scheme` in `[Executioner]` or via the `[TimeIntegrator]`
sub-block. The default is `implicit-euler` (backward Euler), which is first-order accurate and
unconditionally stable:

```
[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 50
  scheme = crank-nicolson    # second-order accurate, unconditionally stable
[]
```

Available schemes:

| `scheme`              | Order | Stability         | Notes                          |
|-----------------------|-------|-------------------|--------------------------------|
| `implicit-euler`      | 1     | Unconditional     | Default; robust, diffusive     |
| `crank-nicolson`      | 2     | Unconditional     | Less diffusion; may oscillate  |
| `bdf2`                | 2     | Unconditional     | Good default for smooth fields |
| `explicit-euler`      | 1     | Conditional       | Requires small `dt`; CFL condition |
| `newmark-beta`        | 2     | Tunable           | Structural dynamics            |

---

## 11. The Postprocessors Block

Postprocessors compute scalar quantities from the simulation state at specified times.
Results are printed to the console and (when `csv = true` in `[Outputs]`) written to a CSV file.

Each sub-block in `[Postprocessors]` produces one scalar value per output step:

```
[Postprocessors]
  [avg_temp]
    type = ElementAverageValue
    variable = temperature
  []
  [max_temp]
    type = NodalMaxValue
    variable = temperature
  []
[]
```

### Common Postprocessors Reference

| Type                                  | Description                                           | Required parameters              |
|---------------------------------------|-------------------------------------------------------|----------------------------------|
| `ElementAverageValue`                 | Volume-weighted average of a variable over the domain (or a block) | `variable` |
| `ElementIntegralVariablePostprocessor`| Volume integral of a variable                         | `variable`                       |
| `NodalMaxValue`                       | Maximum nodal value of a variable                     | `variable`                       |
| `NodalExtremeValue`                   | Maximum or minimum nodal value (with `value_type`)    | `variable`, `value_type`         |
| `SideIntegralVariablePostprocessor`   | Integral of a variable over a sideset                 | `variable`, `boundary`           |
| `SideAverageValue`                    | Average of a variable over a boundary                 | `variable`, `boundary`           |
| `SideDiffusiveFluxIntegral`           | Integral of diffusive flux through a boundary         | `variable`, `boundary`, `diffusivity` |
| `PointValue`                          | Variable value at a specific point                    | `variable`, `point`              |
| `NumNonlinearIterations`              | Number of NL iterations in the last time step         | (none)                           |
| `NumLinearIterations`                 | Number of linear iterations in the last time step     | (none)                           |
| `TimestepSize`                        | Current time step size `dt`                           | (none)                           |
| `ElementAverageMaterialProperty`      | Average of a material property over elements          | `mat_prop`                       |

Examples:

```
[Postprocessors]
  [avg_temperature]
    type = ElementAverageValue
    variable = temperature
    block = 'fuel_block'     # restrict to a subdomain
  []

  [max_temperature]
    type = NodalMaxValue
    variable = temperature
  []

  [outlet_heat_flux]
    type = SideIntegralVariablePostprocessor
    variable = temperature
    boundary = right
  []

  [center_temperature]
    type = PointValue
    variable = temperature
    point = '0.5 0.5 0.0'
  []

  [nl_iters]
    type = NumNonlinearIterations
    execute_on = 'timestep_end'
  []
[]
```

### Execute-On Timing

The `execute_on` parameter controls when a postprocessor is evaluated. Common values:

| `execute_on` value   | When it runs                                               |
|----------------------|------------------------------------------------------------|
| `initial`            | Before the first time step (at t=0, after ICs applied)    |
| `timestep_begin`     | At the start of each time step, before the solve          |
| `linear`             | After each linear iteration (expensive; use sparingly)    |
| `nonlinear`          | After each nonlinear iteration                            |
| `timestep_end`       | After a successful time step solve                        |
| `final`              | After the last time step                                  |
| `ALWAYS`             | Every time any output is written                          |

The default is usually `timestep_end` which is appropriate for most uses. You can specify
multiple values in a vector:

```
[Postprocessors]
  [solution_residual]
    type = NumNonlinearIterations
    execute_on = 'initial timestep_end'
  []
[]
```

---

## 12. AuxVariables and AuxKernels

Auxiliary variables are **not solved for** in the nonlinear system. Instead, they are computed
explicitly from other variables, material properties, or functions at defined times. They are
used to:

- Visualize intermediate quantities (e.g., stress components, flux vectors)
- Feed data into other objects (e.g., a velocity field computed from pressure gradient)
- Store material property values for output

### AuxVariables

Declared identically to primary variables, but in the `[AuxVariables]` block:

```
[AuxVariables]
  [heat_flux_x]
    order  = CONSTANT
    family = MONOMIAL      # element-centered, appropriate for flux quantities
  []
  [velocity]
    order  = CONSTANT
    family = MONOMIAL_VEC  # vector field, one vector per element
  []
[]
```

`CONSTANT MONOMIAL` is the standard choice for quantities that are naturally element-centered,
such as fluxes and stresses. Nodal aux variables (`LAGRANGE FIRST`) can also be used.

### AuxKernels

Each sub-block in `[AuxKernels]` computes one auxiliary variable using an explicit formula.
AuxKernels run in a separate loop after the nonlinear solve.

**`MaterialRealAux`** — Maps a real-valued material property to an auxiliary variable for
visualization or further use.

- Source: `framework/src/auxkernels/MaterialRealAux.C`

```
[AuxKernels]
  [conductivity_field]
    type = MaterialRealAux
    variable = conductivity_out    # must be declared in [AuxVariables]
    property = thermal_conductivity
    execute_on = 'timestep_end'
  []
[]
```

**`VariableGradientComponent`** — Extracts one component of a variable's gradient.

- Source: `framework/src/auxkernels/VariableGradientComponent.C`

```
[AuxVariables]
  [grad_T_x]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [grad_T_x_kernel]
    type = VariableGradientComponent
    variable = grad_T_x
    gradient_variable = temperature
    component = x            # 'x', 'y', or 'z'
    execute_on = 'timestep_end'
  []
[]
```

**`FunctionAux`** — Sets an aux variable to the value of a Function at each quadrature point.
Useful for comparing against an analytic solution or for setting a prescribed velocity field:

```
[AuxKernels]
  [exact_solution]
    type = FunctionAux
    variable = u_exact
    function = analytic_sol
  []
[]
```

**`ParsedAux`** — Sets an aux variable from a parsed expression (like `ParsedMaterial` but
for aux variables):

```
[AuxKernels]
  [heat_generation]
    type = ParsedAux
    variable = q_dot
    coupled_variables = 'neutron_flux temperature'
    expression = 'sigma_f * neutron_flux * 200e6 / (1 + alpha * (temperature - 300))'
    constant_names  = 'sigma_f alpha'
    constant_expressions = '0.05 0.002'
  []
[]
```

### Execute-On for AuxKernels

Like postprocessors, `execute_on` controls when the aux kernel runs:

```
[AuxKernels]
  [velocity_from_pressure]
    type = DarcyVelocity
    variable = velocity
    pressure = pressure
    execute_on = 'timestep_end'    # compute once after the nonlinear solve
  []
[]
```

From `tutorials/darcy_thermo_mech/step04_velocity_aux/problems/step4.i`:

```
[AuxVariables]
  [velocity]
    order = CONSTANT
    family = MONOMIAL_VEC
  []
[]

[AuxKernels]
  [velocity]
    type = DarcyVelocity
    variable = velocity
    execute_on = timestep_end
    pressure = pressure
  []
[]
```

---

## 13. The Outputs Block

The `[Outputs]` block controls what data MOOSE writes to disk and how often.

### Output Formats

**`exodus`** — Exodus II format (`.e` files). This is the standard output for field data. It
supports multiple time steps in a single file, element blocks, sidesets, and both nodal and
element-centered data. Open in ParaView or VisIt.

```
[Outputs]
  exodus = true
[]
```

**`csv`** — Comma-separated values for postprocessor and scalar data. Each row is one time step;
each column is one postprocessor or scalar variable.

```
[Outputs]
  csv = true
[]
```

**`checkpoint`** — Writes MOOSE checkpoint files that can be used to restart the simulation.
Checkpoints include the full solution state: all field variables, material states, and simulation
metadata.

```
[Outputs]
  checkpoint = true
[]
```

By default, MOOSE keeps the two most recent checkpoint sets. Specify `num_checkpoint_files` to
change this:

```
[Outputs]
  [cp]
    type = Checkpoint
    num_files = 5       # keep last 5 checkpoints
  []
[]
```

**`console`** — Controls what is printed to the terminal. Console output is always on by default.
You can configure it explicitly:

```
[Outputs]
  [screen]
    type = Console
    output_linear = true     # print residual at each linear iteration
    output_nonlinear = true  # print residual at each NL iteration
  []
[]
```

**`VTK`** — VTK XML format (`.pvtu` files), an alternative to Exodus for use with ParaView:

```
[Outputs]
  [vtk_out]
    type = VTKOutput
  []
[]
```

### Output Intervals and file_base

Control how often outputs are written with `interval` (for transient simulations):

```
[Outputs]
  [exodus_out]
    type = Exodus
    interval = 10       # write every 10th time step
  []
  csv = true            # write every step (default interval = 1)
[]
```

Set `file_base` to control the base filename for all outputs:

```
[Outputs]
  file_base = 'my_simulation'    # outputs: my_simulation.e, my_simulation.csv
  exodus = true
  csv    = true
[]
```

The default `file_base` is derived from the input file name (the part before `.i`).

### Full Output Block Example

```
[Outputs]
  file_base = 'heat_conduction_2d'
  execute_on = 'initial timestep_end'    # write initial state and after each step
  exodus    = true
  csv       = true
  checkpoint = true
  [screen]
    type = Console
    output_linear = false
    output_nonlinear = true
  []
[]
```

---

## 14. MultiApps and Transfers

MOOSE supports **multiphysics coupling** through the MultiApp system. A "parent" application runs
a primary simulation and spawns one or more "sub-applications", each running their own input file.
Data is exchanged between parent and sub-apps through Transfer objects.

### MultiApp Types

**`TransientMultiApp`** — Runs sub-apps that advance in time alongside the parent. Each sub-app
executes at the same time step (or sub-cycling).

- Source: `framework/src/multiapps/TransientMultiApp.C`

```
[MultiApps]
  [neutronics]
    type = TransientMultiApp
    app_type = OpenMCCoupling   # the C++ application class name
    input_files = neutronics.i
    positions = '0 0 0'         # where to place the sub-app mesh (relative to parent)
    execute_on = 'timestep_end'
  []
[]
```

**`FullSolveMultiApp`** — Runs sub-apps to full convergence (all time steps) at a single point
in the parent's solve, then stops. Useful for initialization or one-way coupling.

- Source: `framework/src/multiapps/FullSolveMultiApp.C`

```
[MultiApps]
  [steady_state_init]
    type = FullSolveMultiApp
    app_type = MyApp
    input_files = init.i
    execute_on = initial        # run once at the start of the parent simulation
    positions = '0 0 0'
  []
[]
```

Multiple sub-apps at different locations can be declared by providing multiple positions:

```
[MultiApps]
  [pins]
    type = TransientMultiApp
    input_files = 'pin_a.i pin_b.i pin_c.i'
    positions = '0 0 0   0.02 0 0   0.04 0 0'   # 3 sub-apps at different x positions
    execute_on = timestep_end
  []
[]
```

### Transfer Types

Transfers move field data between the parent and sub-apps:

**`MultiAppGeneralFieldNearestLocationTransfer`** — The most general transfer. Maps a variable
from one app to a variable in another using nearest-neighbor interpolation:

```
[Transfers]
  [temperature_to_neutronics]
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = neutronics        # direction: parent -> sub
    source_variable = temperature
    variable = T_fuel
  []

  [power_from_neutronics]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = neutronics      # direction: sub -> parent
    source_variable = power_density
    variable = heat_source
  []
[]
```

**`MultiAppPostprocessorTransfer`** — Transfers a single scalar postprocessor value between apps:

```
[Transfers]
  [max_temp_to_sub]
    type = MultiAppPostprocessorTransfer
    to_multi_app = thermo
    from_postprocessor = max_temperature
    to_postprocessor = inlet_temperature
  []
[]
```

**`MultiAppCopyTransfer`** — Copies a variable between apps that share the same mesh:

```
[Transfers]
  [copy_solution]
    type = MultiAppCopyTransfer
    from_multi_app = fine_mesh_app
    source_variable = u_fine
    variable = u
  []
[]
```

From `test/tests/transfers/from_full_solve/parent.i`:

```
[Transfers]
  [from_full]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = full_solve
    source_variable = u
    variable = from_full
  []
[]
```

### Picard (Fixed-Point) Iteration

When physics are tightly coupled (each app's solution affects the other), you need
Picard (fixed-point) iterations to converge the coupling. Enable this in the executioner:

```
[Executioner]
  type = Transient
  num_steps = 20
  dt = 0.1

  fixed_point_max_its = 30    # maximum Picard iterations per time step
  nl_abs_tol = 1e-14
  fixed_point_abs_tol = 1e-10  # convergence criterion for the Picard loop
  fixed_point_rel_tol = 1e-8
[]
```

From `test/tests/multiapps/picard/picard_adaptive_parent.i`:

```
[Executioner]
  type = Transient
  num_steps = 20
  dt = 0.1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  fixed_point_max_its = 30
  nl_abs_tol = 1e-14
[]
```

---

## 15. Running Simulations

### Basic Command

```bash
./myapp-opt -i input.i
```

Replace `myapp-opt` with your compiled binary name. The `-opt` suffix indicates a fully optimized
build. Development builds are `-dbg` (debug, very slow but informative) and `-oprof` (profiling).

### Parallel Execution with MPI

MOOSE is designed for parallel execution using MPI. Distribute the problem across `N` processors:

```bash
mpirun -n 4 ./myapp-opt -i input.i
```

```bash
mpirun -n 16 ./myapp-opt -i input.i
```

MOOSE automatically partitions the mesh and distributes work. Most large simulations scale well
to hundreds or thousands of cores. The output files are still merged into a single Exodus file
(or individual per-processor files for very large cases with `Nemesis` output).

For cluster job schedulers, replace `mpirun` with the appropriate launcher:

```bash
# SLURM
srun -n 64 ./myapp-opt -i input.i

# PBS
mpiexec -n 64 ./myapp-opt -i input.i
```

### Useful Command-Line Flags

| Flag                        | Effect                                                        |
|-----------------------------|---------------------------------------------------------------|
| `--mesh-only`               | Generate and write the mesh to an Exodus file, then exit. No solve is performed. |
| `--dump`                    | Print all available input parameters (with defaults) for every registered type, then exit. |
| `--dump [type_name]`        | Print parameters for a specific type, e.g. `--dump DirichletBC` |
| `--show-input`              | Print the fully resolved input file (with all substitutions applied) to the console, then exit. |
| `--help` / `-h`             | Print help and all command-line options, then exit.           |
| `--recover`                 | Restart from the latest checkpoint (see Section 17).          |
| `--recover path/to/cp`      | Restart from a specific checkpoint directory.                 |
| `-o output_base`            | Override `file_base` for output filenames.                    |
| `--n-threads=N`             | Use N OpenMP threads per MPI rank (thread-parallel assembly). |
| `--no-gdb-backtrace`        | Disable gdb backtrace on errors (useful when gdb is unavailable). |

### Reading Console Output

A typical transient run produces output like this:

```
Time Step 1, time = 0.1, dt = 0.1

 0 Nonlinear |R| = 3.524687e+02
      0 Linear |R| = 3.524687e+02
      1 Linear |R| = 2.104523e+01
      2 Linear |R| = 1.842317e+00
      3 Linear |R| = 1.032845e-01
      4 Linear |R| = 8.124532e-03
 1 Nonlinear |R| = 8.124532e-03
      0 Linear |R| = 8.124532e-03
      1 Linear |R| = 2.341872e-04
      2 Linear |R| = 8.423541e-06
 2 Nonlinear |R| = 8.423541e-06
Converged at iteration 2. |R| = 8.423541e-06 < 1e-05 = nl_rel_tol * |R_0|

Postprocessor Values:
+----------------+----------------+
| time           | avg_temp       |
+----------------+----------------+
| 0.100000       | 312.345678     |
+----------------+----------------+
```

Key items to monitor:

- **Nonlinear `|R|`** — Each row is one Newton iteration. Should decrease by roughly an order of
  magnitude per iteration for well-conditioned problems. If it stalls or grows, the problem may be
  poorly conditioned or your time step may be too large.
- **Linear `|R|`** — Each indented row is one Krylov (GMRES) iteration inside a Newton step.
  Should decrease steadily. If it requires many iterations, consider a stronger preconditioner.
- **Convergence message** — Tells you which tolerance criterion was satisfied.
- **Postprocessor table** — Scalar values at the end of each step. Useful for quick sanity checks.

If a time step fails to converge, MOOSE will cut the time step (when using `IterationAdaptiveDT`
or similar) and retry. A completely failed solve prints an error and the run may stop.

---

## 16. Postprocessing Output

### Opening Exodus Files in ParaView

The primary field output is the Exodus `.e` file. To open it:

1. Launch ParaView
2. File > Open > select `simulation_out.e`
3. Click **Apply** in the Properties panel to load the data
4. Use the **Pipeline Browser** to add filters (e.g., Clip, Warp by Scalar, Contour)
5. Select the variable to display from the dropdown menu in the toolbar
6. Use the time controls to step through time steps

Useful ParaView operations for MOOSE output:

- **Warp by Scalar** — Deform the mesh surface by a field value (good for visualizing
  temperature distributions)
- **Contour** — Draw isosurface/isoline at specified values
- **Clip** — Cut the domain to see interior values in 3D meshes
- **Plot Over Line** — Extract a 1D profile along an arbitrary line through the domain
- **Plot Data Over Time** — Plot a postprocessor value over time (for ParaView's CSV reader)

To visualize element-centered (CONSTANT MONOMIAL) data in ParaView, MOOSE automatically
projects it to nodes in the Exodus output so ParaView can render it smoothly.

### Reading CSV Files with Python

When `csv = true` in `[Outputs]`, MOOSE writes a CSV file with columns for time and each
postprocessor. The file is named `{file_base}.csv` by default.

Example CSV content:

```
time,avg_temperature,max_temperature,heat_flux_out
0,300.000000,300.000000,0.000000
0.1,301.234567,315.678901,1234.567
0.2,302.456789,328.901234,2345.678
...
```

Reading and plotting with pandas and matplotlib:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV
df = pd.read_csv('heat_conduction_out.csv')

# Basic plot
fig, axes = plt.subplots(2, 1, figsize=(8, 6))

axes[0].plot(df['time'], df['avg_temperature'], label='Average Temperature')
axes[0].plot(df['time'], df['max_temperature'], label='Max Temperature')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('Temperature (K)')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(df['time'], df['heat_flux_out'], color='red')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Heat Flux Out (W/m^2)')
axes[1].grid(True)

plt.tight_layout()
plt.savefig('results.png', dpi=150)
plt.show()
```

### Using mooseutils

The MOOSE repository includes a Python package `mooseutils` (in `python/mooseutils/`) that
provides convenience readers for MOOSE output formats.

Reading an Exodus file:

```python
import sys
sys.path.insert(0, '/path/to/moose/python')
import mooseutils

# Read exodus data
data = mooseutils.VectorPostprocessorReader('simulation_out.e')
```

The `mooseutils` package also includes `MooseDataFrame`, which reads MOOSE CSV output directly
into a pandas DataFrame with added metadata handling.

---

## 17. Restarting Simulations

MOOSE supports checkpointing and restarting long simulations. If a job is interrupted (power loss,
wall time limit, etc.), you can resume from the last checkpoint without losing progress.

### Setting Up Checkpoints

Add checkpoint output to your `[Outputs]` block:

```
[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  checkpoint = true    # write checkpoint files every time step
[]
```

For fine-grained control:

```
[Outputs]
  [cp]
    type = Checkpoint
    execute_on = 'timestep_end'
    num_files = 2       # keep the last 2 checkpoints (default)
  []
[]
```

From `test/tests/outputs/checkpoint/checkpoint_parent.i`:

```
[Outputs]
  execute_on = 'timestep_end'
  checkpoint = true
[]
```

MOOSE writes checkpoint files into a directory named `{file_base}_cp/`. Each checkpoint consists
of several files (mesh data, solution vectors, and a metadata file).

### Recovering from a Checkpoint

To restart a simulation from the latest checkpoint, run with the `--recover` flag:

```bash
./myapp-opt -i input.i --recover
```

MOOSE automatically finds the most recent valid checkpoint in `{file_base}_cp/` and restarts
from that state. The simulation continues from where it stopped.

To recover from a specific checkpoint:

```bash
./myapp-opt -i input.i --recover path/to/checkpoint_cp/0050
```

The number `0050` refers to the time step number in the checkpoint filename.

### Checkpoint File Structure

After a run, checkpoint files appear as:

```
simulation_cp/
  0010_mesh.cpr          # mesh topology at step 10
  0010.xdr               # solution data at step 10
  0020_mesh.cpr
  0020.xdr
  LATEST                 # symlink or file pointing to the latest checkpoint
```

MOOSE keeps only `num_files` (default 2) checkpoint sets at a time, automatically deleting older
ones to save disk space. If you need to keep more checkpoints (e.g., for analysis), increase
`num_files`.

### Restart vs Recover

| Action      | Flag          | Use case                                               |
|-------------|---------------|--------------------------------------------------------|
| **Recover** | `--recover`   | Continue an interrupted run from where it stopped. The input file and all parameters stay the same. |
| **Restart** | (via `[Mesh]` + `use_for_exodus_restart`) | Start a new simulation (possibly with different physics) using the solution from a previous run as the initial condition. |

For a restart from an Exodus file:

```
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = previous_run_out.e
    use_for_exodus_restart = true
  []
[]

[Variables]
  [temperature]
    initial_from_file_var = temperature   # read IC from the Exodus variable
    initial_from_file_timestep = LATEST   # use the last time step
  []
[]
```

---

## Appendix A: Complete Example — Transient Heat Conduction

The following is a complete, working input file for transient heat conduction on a 2D rectangle,
based on patterns verified in the MOOSE test suite and tutorial files.

```
# Transient 2D heat conduction on a rectangular domain
# PDE: rho*cp * dT/dt = nabla dot (k nabla T) + Q
#
# Domain: 0.304 m x 0.0257 m
# Left BC: T = 350 K (hot inlet)
# Right BC: T = 300 K (cold outlet)
# Top/Bottom: no flux (natural BC)
# Initial condition: T = 300 K (uniform)

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 10
    xmax = 0.304
    ymax = 0.0257
  []
[]

[Variables]
  [temperature]
    initial_condition = 300   # uniform IC at 300 K
  []
[]

[Kernels]
  [heat_conduction_time_deriv]
    type = TimeDerivative
    variable = temperature
  []
  [heat_conduction]
    type = MatDiffusion
    variable = temperature
    diffusivity = k_over_rho_cp   # combined property: k/(rho*cp)
  []
[]

[BCs]
  [inlet]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 350
  []
  [outlet]
    type = DirichletBC
    variable = temperature
    boundary = right
    value = 300
  []
[]

[Materials]
  [steel]
    type = GenericConstantMaterial
    prop_names  = 'k_over_rho_cp'
    prop_values = '4.838e-6'   # k/(rho*cp) = 18/(8000*0.466) m^2/s
  []
[]

[Postprocessors]
  [avg_temp]
    type = ElementAverageValue
    variable = temperature
  []
  [max_temp]
    type = NodalMaxValue
    variable = temperature
  []
[]

[Executioner]
  type = Transient
  num_steps = 50
  dt = 10.0
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  file_base = 'heat_conduction_2d'
  exodus = true
  csv    = true
[]
```

---

## Appendix B: Complete Example — Steady Diffusion with Material Properties

A steady-state diffusion problem with spatially varying diffusivity, using the AD (automatic
differentiation) kernel chain:

```
# Steady diffusion: -nabla dot (D(x) nabla u) = 0
# D = 1.0 + 0.5*x (varies linearly from 1.0 to 1.5 across the domain)
# Left BC: u = 1
# Right BC: u = 0

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 50
  []
[]

[Variables]
  [u]
    family = LAGRANGE
    order  = FIRST
  []
[]

[Kernels]
  [diff]
    type = ADMatDiffusion
    variable = u
    diffusivity = diffusivity
  []
[]

[BCs]
  [left]
    type = ADDirichletBC
    variable = u
    boundary = left
    value = 1
  []
  [right]
    type = ADDirichletBC
    variable = u
    boundary = right
    value = 0
  []
[]

[Materials]
  [diff_mat]
    type = ADGenericFunctionMaterial
    prop_names  = 'diffusivity'
    prop_values = 'diffusivity_func'
  []
[]

[Functions]
  [diffusivity_func]
    type = ParsedFunction
    expression = '1.0 + 0.5 * x'
  []
[]

[Postprocessors]
  [avg_u]
    type = ElementAverageValue
    variable = u
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

---

## Appendix C: Troubleshooting Common Issues

### The solver is not converging

**Symptoms**: Nonlinear `|R|` does not decrease, or oscillates, or converges to a large residual.

**Steps to diagnose**:
1. Run with a coarser mesh to check if the problem is well-posed.
2. Check your BC values and ensure you have physically meaningful initial conditions.
3. For transient problems, try reducing `dt` by a factor of 10.
4. Switch to `solve_type = NEWTON` with a direct solver (`-pc_type lu`) to check if the Jacobian
   is correct. If it converges with LU but not with AMG, the issue is in the preconditioner choice.
5. Run `--show-input` to confirm the input is being parsed as expected.
6. Add `output_nonlinear = true` to the Console output to see the residual at every NL iteration.

### Boundary condition not being applied

**Symptoms**: Field values at a boundary are not what you specified.

**Steps**:
1. Run with `--mesh-only` and inspect the mesh in ParaView to confirm the sideset names/IDs
   are correct.
2. Check for typos in the `boundary = ` parameter.
3. Verify that the sideset name matches exactly (case-sensitive) what is in the mesh.

### "Unknown parameter" error

**Example**: `Error: Unknown parameter 'diffusivity' in block 'Kernels/diff'`

**Fix**: Use `./myapp-opt --dump MatDiffusion` to see all valid parameters for that type.

### Material property not found

**Example**: `Error: Material property 'k' not defined on block 0`

**Fix**: Ensure the material providing `k` either has no block restriction (applies everywhere)
or has `block = 0` (or the appropriate block name). A common mistake is restricting a material
to one block while a kernel on another block requests the same property name.

### Output files are empty or missing

**Fix**: Check that `execute_on` in `[Outputs]` includes `timestep_end`. The default is usually
correct, but if you have overridden it with only `final`, you will only get one output at the end.

### Running out of memory in parallel

**Steps**:
1. Use `type = DistributedMesh` (enabled by default for large meshes) rather than replicated mesh.
2. Reduce the mesh size or use adaptive mesh refinement.
3. Use `--n-threads=N` to take advantage of shared memory within each node.

---

## Appendix D: Input File Syntax Quick Reference

```
# Comments start with #

# Top-level variable substitution
variable_name = value

[GlobalParams]
  variable = u    # applies to all objects unless overridden
[]

[Mesh]
  [generator_name]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []
[]

[Variables]
  [var_name]
    family = LAGRANGE   # optional; default is LAGRANGE
    order  = FIRST      # optional; default is FIRST
    initial_condition = 0.0   # optional inline IC
  []
[]

[Kernels]
  [kernel_name]
    type = SomeKernelType
    variable = var_name
    # ... type-specific parameters
  []
[]

[BCs]
  [bc_name]
    type = DirichletBC
    variable = var_name
    boundary = left
    value = 0
  []
[]

[ICs]
  [ic_name]
    type = ConstantIC
    variable = var_name
    value = 300
  []
[]

[Materials]
  [mat_name]
    type = GenericConstantMaterial
    block = subdomain_name   # optional block restriction
    prop_names  = 'k'
    prop_values = '1.0'
  []
[]

[Functions]
  [func_name]
    type = ParsedFunction
    expression = 'sin(x) * cos(t)'
  []
[]

[Postprocessors]
  [pp_name]
    type = ElementAverageValue
    variable = var_name
    execute_on = 'initial timestep_end'
  []
[]

[AuxVariables]
  [aux_name]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [auxk_name]
    type = MaterialRealAux
    variable = aux_name
    property = k
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = Steady               # or Transient
  solve_type = PJFNK          # or NEWTON, LINEAR
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  # Transient-only parameters:
  num_steps = 100
  dt = 0.01
  end_time = 10.0
  scheme = implicit-euler

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  file_base = 'my_simulation'
  execute_on = 'initial timestep_end'
  exodus    = true
  csv       = true
  checkpoint = true
[]
```
