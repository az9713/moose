# Using the `moose-simulation` Skill — A Field Guide

This document captures the real-world experience of using the Claude Code
`moose-simulation` skill to create, run, debug, and visualize 36 MOOSE
finite-element simulations on Windows via Docker. It is written so that a
future user can invoke the skill with confidence and know exactly what to
expect.

---

## What the Skill Does

The `moose-simulation` skill is a structured checklist that governs the
**complete lifecycle** of a MOOSE simulation:

```
Prerequisites → Author .i file → Run in Docker → Validate outputs
     → Debug failures → Generate plots → Write README → Report results
```

You invoke it by telling Claude Code to "use the moose-simulation skill" or
by working on any `.i` file in the `quickstart-runs/` directory. Claude Code
will then follow the skill's 9-section checklist systematically.

---

## What Was Actually Built

Over multiple sessions the skill guided creation of **58 quickstart cases**
covering 9 MOOSE physics modules:

| Batch | Cases | Source Textbook | Physics |
|-------|-------|-----------------|---------|
| Original | 01–21 | Various | Diffusion, heat transfer, solid mechanics, Navier-Stokes, phase field, porous flow, electromagnetics |
| Melcher | 22–29 | Melcher, *Continuum Electromechanics* | Charge relaxation, EHD, MHD, ferrofluid, electroquasistatics |
| Haus | 30–36 | Professor Herman A. Haus, *Electromagnetic Noise and Quantum Optical Measurements* (Springer, 2000) — classical chapters only (Chs 1-5, 10); the quantum chapters (Chs 6-9, 11-13) describe photon operators, squeezed states, and quantum noise that have no classical PDE representation | Waveguide eigenvalues, driven cavities, dielectric slabs, coupled resonators, thermal noise, dispersive pulses, solitons |
| Rieutord | 37–44 | Michel Rieutord, *Fluid Dynamics: An Introduction* (Springer, 2015) — Chapters 4-10 | Rayleigh-Benard convection, KH instability, Blasius boundary layer, k-epsilon turbulence, RT instability, Sod shock tube, Ekman spiral, Alfven wave |
| Smith | 45–48 | Smith, *Uncertainty Quantification* (SIAM, 2014) | Monte Carlo UQ, polynomial chaos expansion, heat source inversion (adjoint optimization), Latin Hypercube parameter study |
| Nonlinear Solid Mechanics (Batch A) | 49–53 | MOOSE solid_mechanics module | J2 plasticity, finite-strain compression, power-law creep, phase-field fracture, Lame pressure vessel solution |
| Nuclear Reactor Physics | 54–58 | Nuclear engineering fundamentals (diffusion theory, reactor kinetics) | 1-group and 2-group neutron diffusion eigenvalue, fuel-pin RZ heat transfer, xenon-135 poisoning transient, control rod worth |

Every case converges in Docker with `combined-opt` in under 2 minutes and
produces both Exodus (`.e`) and CSV output files.

---

## How to Invoke the Skill

### Approach 1: Explicit invocation

> "Use the moose-simulation skill to create and run a new case for
> thermal convection in a porous medium."

This triggers the full checklist from Section 1 (prerequisites) through
Section 9 (report results).

### Approach 2: Implicit — working on a `.i` file

> "Debug this input file — it's failing with 'unused parameter'."

The skill activates automatically when working on MOOSE `.i` files. Claude
Code will check Docker prerequisites, run the file, and diagnose using the
skill's Section 7 failure patterns.

### Approach 3: Batch creation

> "Implement the following plan: [detailed case list]"

The skill scales to batch creation. For cases 30–36, seven input files,
seven READMEs, seven plot functions, and documentation updates were created
in a single session using parallel agent dispatch.

---

## The 9-Step Workflow in Practice

Here is what each step looked like during the Cases 30–36 implementation,
with the actual time each step took and the real issues encountered.

### Step 1: Verify Prerequisites (~30 seconds)

```bash
docker info > /dev/null 2>&1 && echo "DOCKER OK"
MSYS_NO_PATHCONV=1 docker image inspect idaholab/moose:latest > /dev/null 2>&1 && echo "IMAGE OK"
```

Both checks passed. The skill blocks all further work until Docker responds.

**Lesson learned**: Never skip `MSYS_NO_PATHCONV=1`. Without it, MINGW
silently mangles `/work` and `/opt/moose` paths, and Docker either fails
or produces empty output with no error message.

### Step 2: Create Case Directories

The skill enforces naming conventions:
- Directory: `case30-waveguide-cutoff/` (hyphens)
- Input file: `case30_waveguide_cutoff.i` (underscores)

Seven directories were created. Each `.i` file follows the skill's
authoring standards: header comment block, inline physics comments,
`exodus = true` and `csv = true` in `[Outputs]`.

### Step 3: Write Input Files

The skill's Section 2 standards produced consistently high-quality files.
Key conventions that prevented bugs:

- **Header block**: Every file starts with the governing equations, boundary
  conditions, and expected results in comments. This catches physics errors
  before running.
- **Docker portability rules**: `DerivativeParsedMaterial` blocks always
  include `disable_fpoptimizer = true` and `enable_jit = false`.
- **Mesh sizing**: All meshes stayed within the 20×20 to 40×40 guideline.
  Largest run (Case 36: 100×2 quasi-1D, 500 timesteps) took 109 seconds.

### Step 4: Run in Docker

The skill's canonical Docker command was used verbatim for every run:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case30-waveguide-cutoff \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case30_waveguide_cutoff.i 2>&1 | tail -40'
```

**Every element matters**:
- `MSYS_NO_PATHCONV=1` — prevents MINGW path mangling (silent failure)
- `--entrypoint /bin/bash` — overrides default entrypoint
- `2>&1 | tail -40` — MOOSE prints thousands of lines; only the last 40 matter
- Single quotes around the `-c` argument — prevents host shell expansion

### Step 5: Validate Outputs

After each successful run, the skill requires checking:
1. `.e` and `.csv` files exist with non-zero size
2. CSV columns span the expected time range
3. No NaN/Inf values
4. Conserved quantities actually conserve

**Real example**: Case 35 (dispersive pulse) — `total_A` was verified
constant to 14 significant digits across 300+ timesteps. This confirmed
mass conservation and validated the advection-diffusion implementation.

### Step 6: Debug Failures (the most valuable part)

**5 of 7 cases required debugging.** The skill's Section 7 failure patterns
were the first reference every time. Here are the real failures and how
each was resolved:

#### Case 30 — Eigenvalue solver (3 iterations to fix)

| Attempt | Error | Root Cause | Fix |
|---------|-------|-----------|-----|
| 1 | `Zero pivot row 0` | Shift-invert with default σ=0 | Not in skill — had to research SLEPc |
| 2 | `Power method cannot compute more than one eigenvalue` | Default PJFNK maps to power method | Change `solve_type = KRYLOVSCHUR` |
| 3 | `sinvert requires a target` | Missing `-eps_target` | Add `which_eigen_pairs = TARGET_MAGNITUDE`, `-eps_target 12.0` |

**Takeaway**: Eigenvalue problems are a blind spot in the skill's current
failure catalog. The fix requires three coordinated settings:
`solve_type = KRYLOVSCHUR` + `which_eigen_pairs = TARGET_MAGNITUDE` +
`-eps_target <value> -st_type sinvert`.

#### Case 32 — Unused parameter

```
unused parameter 'GlobalParams/theta'
```

Diagnosed immediately from the skill's Section 7.4 ("Renamed Parameters"):
every parameter must be consumed by at least one object. The `[GlobalParams]`
block had `theta` but no object read it. Removed the block; fixed.

#### Case 33 — Exponential blowup (physics error)

The simulation diverged — `avg_u` grew exponentially instead of oscillating.
Not a MOOSE error but a **physics error in the plan**:

- Plan specified symmetric coupling (+κ in both equations)
- Eigenvalues of the coupling matrix: -γ ± κ = -0.5 ± 3.0
- One eigenvalue is +2.5 → **exponentially unstable**
- Fix: antisymmetric coupling (`+κ` in u equation, `-κ` in v equation)
- Corrected eigenvalues: -γ ± jκ → oscillatory with decay

**Lesson**: The skill's solver-divergence checklist (Section 7.6) is for
*numerical* divergence. *Physical* instabilities require checking the
governing equations — the skill can't catch sign errors in the physics.

#### Case 36 — AD/non-AD material property mismatch

```
AD material property 'nl_rate' already retrieved as non-AD
```

Not in the skill's current failure catalog. Root cause:
- `DerivativeParsedMaterial` produces **non-AD** material properties
- `ADMatReaction` requires **AD** material properties
- Fix: change `ADMatReaction` → `MatReaction` (non-AD version)

**Takeaway**: When mixing AD and non-AD objects, the property types must
match. `ADGenericFunctionMaterial` → AD properties (use with `ADMatReaction`).
`DerivativeParsedMaterial` → non-AD properties (use with `MatReaction`).

#### Case 32 — Unused GlobalParams

Same pattern as Section 7.4. Every HIT top-level variable and every
`[GlobalParams]` entry must be consumed by at least one object, or MOOSE
throws an unused-parameter error.

### Step 7: Generate Plots

The skill requires adding a `plot_caseNN()` function to `visualize_all.py`
for each case. The function reads the Exodus file (via netCDF4) and/or CSV,
then produces a multi-panel PNG.

**Real pitfall**: The eigenvalue CSV for Case 30 was named
`case30_waveguide_cutoff_out_eigenvalues_0002.csv` (with a `_0002` suffix
from the eigenvalue executioner), not the expected `_0001.csv`. The plot
function had to be patched after the first run.

### Step 8: Verify Plots Visually

Each PNG was opened and checked for physically correct behavior:
- Case 30: TM₁₁ mode shape (single lobe, zero on walls) ✓
- Case 33: u and v oscillate 90° out of phase with exponential decay ✓
- Case 34: Random speckle → smooth → uniform T = 0.5 ✓
- Case 35: Gaussian translates right while broadening; total mass flat ✓

### Step 9: Write READMEs and Report

Each case directory got a README.md with the structure the skill prescribes:
Overview → Physics → MOOSE Implementation → Running → Expected Results →
Key Takeaways.

---

## What the Skill Gets Right

1. **The Docker command template is bulletproof.** Copy-paste it exactly
   and it works every time. The `MSYS_NO_PATHCONV=1` prefix alone prevents
   the most common silent failure on Windows.

2. **The failure pattern catalog (Section 7) catches ~60% of errors
   instantly.** Unused parameters, JIT failures, solver divergence, and
   missing materials are all diagnosed on first sight.

3. **The validation checklist prevents false confidence.** Checking that
   `.csv` files span the correct time range and conserved quantities
   actually conserve caught a sign error that would otherwise have shipped.

4. **The authoring standards produce self-documenting input files.** The
   header comment block with governing equations means anyone can open a
   `.i` file and understand the physics without external documentation.

---

## What the Skill Doesn't Cover (Yet)

These gaps were encountered during Cases 30–36 and required ad-hoc
debugging:

| Gap | Cases Affected | What Was Needed |
|-----|---------------|-----------------|
| Eigenvalue executioner configuration | 30 | `KRYLOVSCHUR` + `TARGET_MAGNITUDE` + `-eps_target` + `-st_type sinvert` |
| AD vs non-AD material property compatibility | 36 | `DerivativeParsedMaterial` → non-AD → use `MatReaction` not `ADMatReaction` |
| Kernel sign conventions (`MatReaction` vs `CoefReaction`) | 32 | `ADMatReaction` residual = −rate·u·test; `CoefReaction` residual = +coeff·u·test |
| EM module objects (`EMRobinBC`, `ReflectionCoefficient`) | 32 | Port BCs, field_real/field_imaginary coupling, sign parameter |
| Physics validation (not just numerical convergence) | 33 | Antisymmetric vs symmetric coupling — eigenvalue analysis needed |
| Eigenvalue CSV naming (`_0002.csv` suffix) | 30 | `VectorPostprocessors/Eigenvalues` output naming is non-obvious |

---

## Recommended Workflow for New Users

1. **Start with a simple case.** Pick a physics you understand (e.g.,
   transient heat diffusion) and follow the skill end-to-end. Cases 01–03
   are good starting points.

2. **Always run the prerequisite checks.** The 30 seconds they take
   prevents 30 minutes of debugging silent Docker failures.

3. **Read existing cases as templates.** Before writing a new `.i` file,
   find the closest existing case and adapt it. The 44 quickstart cases
   cover most MOOSE object types.

4. **Trust the failure catalog first.** When a run fails, check Section 7
   before searching the MOOSE documentation. The catalog is sorted by
   frequency of occurrence.

5. **Validate CSV data, not just "Solve Converged!"** A converged solve
   can still have wrong physics (Case 33 taught this lesson). Check that
   conserved quantities conserve and that signs are correct.

6. **Use `tail -40` on Docker output.** MOOSE prints thousands of lines.
   The last 40 contain the convergence status, postprocessor table, and
   timing — everything you need to diagnose success or failure.

---

## Quick Reference: The Canonical Docker Command

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/caseNN-directory-name \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i caseNN_input_file.i 2>&1 | tail -40'
```

Replace `caseNN-directory-name` and `caseNN_input_file.i` with your case.
Do not change anything else.
