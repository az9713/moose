# Case 28: Two-Way Joule Heating — Temperature-Dependent Conductivity

## Overview

This case extends **Case 17 (Joule Heating)** by making the electrical
conductivity σ a function of temperature. In Case 17 σ is constant — the
electric field drives a heat source but the resulting temperature has no
effect on σ (one-way coupling). Here the coupling runs in both directions:

```
T changes → σ(T) changes → V distribution changes → Q changes → T changes → ...
```

This two-way coupling is physically realistic for metallic conductors. Metals
follow a well-known rule: their resistivity ρ_e = 1/σ rises roughly linearly
with temperature (phonon scattering increases as atoms vibrate more). The model
captured here is:

```
σ(T) = σ₀ / (1 + α·(T - T_ref))
```

where σ₀ is the reference conductivity at T_ref and α is the temperature
coefficient of resistivity. As temperature rises, the denominator grows, so
σ decreases. Lower conductivity means less current for the same voltage, which
means less Joule heating — a **negative feedback** that self-limits the
temperature rise.

This contrasts with **semiconductors**, where conductivity increases with
temperature (more thermally excited carriers). That positive feedback can lead
to **thermal runaway**: higher T increases σ, which increases Q, which raises T
further. The metallic negative-feedback case in this simulation is inherently
stable.

The reference physics is Melcher, *Continuum Electromechanics* (MIT Press,
1981), Chapter 10, §10.1–10.3.

---

## The Physics

### Governing Equations

Two coupled PDEs are solved simultaneously on the domain [0,2] × [0,1]:

**Current conservation (electric potential):**

```
-div(σ(T)·grad(V)) = 0
```

This is the Laplace equation for V when σ is uniform, but with a
T-dependent σ it becomes a nonlinear equation. Changes in T alter σ at every
point, which reshapes the current density J = -σ(T)·grad(V) and therefore the
spatial distribution of V.

**Transient heat equation with Joule source:**

```
ρ·cp·∂T/∂t = div(k·grad(T)) + σ(T)·|grad(V)|²
```

The volumetric heat source Q = σ(T)·|grad(V)|² now depends on T in two ways:
through σ(T) directly, and through V (which depends on T via the first
equation). This is what makes the problem truly two-way coupled.

### The Metallic Conductivity Model

```
σ(T) = σ₀ / (1 + α·(T - T_ref))
```

Parameters used in this case:

| Symbol  | Value      | Units | Meaning                                     |
|---------|------------|-------|---------------------------------------------|
| σ₀      | 1e6        | S/m   | Conductivity at reference temperature       |
| α       | 0.004      | 1/K   | Temperature coefficient of resistivity      |
| T_ref   | 300        | K     | Reference temperature                       |

At 25 K above T_ref the denominator equals 1 + 0.004 * 25 = 1.1, so σ drops
to 91% of its reference value. At 250 K above T_ref it drops to 50%. The
model is valid when T - T_ref is small enough that the denominator stays
positive (i.e., T < T_ref + 1/α = 300 + 250 = 550 K for these parameters).

### Negative Feedback Mechanism

Consider what happens when a local hot spot develops at the centre of the
conductor (x = 1 m):

1. T rises at the centre.
2. σ drops at the centre (metallic behaviour).
3. The centre region becomes more resistive.
4. Current density J = -σ·grad(V) redistributes: current avoids the resistive
   hot spot and spreads toward cooler (more conductive) regions near the top
   and bottom boundaries.
5. The Joule heating Q = σ(T)·|grad(V)|² decreases at the hot spot because σ
   is lower there. This damps the local heat source — negative feedback.

In a voltage-driven geometry (fixed V boundaries, as here), the Joule heating
at any point is Q = σ(T)·|grad(V)|². As σ decreases in the hot region, Q
decreases there, reducing the local heat source. The system approaches a
stable steady state where the temperature distribution and conductivity
distribution are mutually consistent.

### Comparison: Constant σ (Case 17) vs. σ(T) (Case 28)

| Property                  | Case 17 (constant σ)          | Case 28 (metallic σ(T))         |
|---------------------------|-------------------------------|---------------------------------|
| Coupling direction        | One-way: V drives T           | Two-way: V drives T, T drives V |
| V distribution            | Linear in x (always)          | Slightly non-linear (T-warped)  |
| Q distribution            | Uniform (same everywhere)     | Lower in hotter central region  |
| Temperature rise rate     | Fastest (full σ always)       | Slower (σ drops as T rises)     |
| Steady state              | Very high (T grows unbounded) | Lower, finite (negative FB)     |
| Newton nonlinearity       | Mild (linear Laplace for V)   | Stronger (σ couples both PDEs)  |

---

## Input File Walkthrough

The input file is `case28_twoway_joule_heating.i`. It is identical to
`case17_joule_heating.i` except for the `[Materials]` block.

### Block: `[Variables]`

Same as Case 17: `V` (electric potential, initial 0 V) and `T` (temperature,
initial 300 K).

### Block: `[Kernels]`

Same four kernels as Case 17:

- `V_diff` — `ADHeatConduction` on V, reading `electrical_conductivity`
- `T_time` — `ADHeatConductionTimeDerivative` for ρ·cp·∂T/∂t
- `T_diff` — `ADHeatConduction` for div(k·grad(T))
- `T_joule` — `ADJouleHeatingSource` reading `electric_field_heating`

The only change is that `electrical_conductivity` is now T-dependent (provided
by the new material below). The kernels themselves are unchanged.

### Block: `[BCs]`

Same as Case 17: V = 10 on left, V = 0 on right, T = 300 K on both electrodes.

### Block: `[Materials]` — KEY DIFFERENCE

This is the only substantive change from Case 17.

**Case 17 (constant conductivity):**

```
[electrical]
  type        = ADGenericConstantMaterial
  prop_names  = 'electrical_conductivity'
  prop_values = '1e6'
[]
```

**Case 28 (temperature-dependent conductivity via tabulated interpolation):**

```
[sigma_of_T]
  type         = ADPiecewiseLinearInterpolationMaterial
  property     = electrical_conductivity
  variable     = T
  x            = '250  300   350     400     450     500     550     600     700'
  y            = '1.25e6  1e6 833333  714286  625000  555556  500000  454545  384615'
  extrapolation = true
[]
```

`ADPiecewiseLinearInterpolationMaterial` evaluates a piecewise-linear function
of a coupled variable at every quadrature point. The tabulated (x, y) pairs
were precomputed from the analytic formula σ₀ / (1 + α·(T - T_ref)):

| T (K) | σ (S/m)  | Formula check                              |
|-------|----------|--------------------------------------------|
| 250   | 1 250 000 | 1e6 / (1 + 0.004·(250-300)) = 1e6 / 0.8  |
| 300   | 1 000 000 | 1e6 / (1 + 0.004·0) = 1e6 / 1.0           |
| 350   | 833 333   | 1e6 / (1 + 0.004·50) = 1e6 / 1.2          |
| 400   | 714 286   | 1e6 / (1 + 0.004·100) = 1e6 / 1.4         |
| 450   | 625 000   | 1e6 / (1 + 0.004·150) = 1e6 / 1.6         |
| 500   | 555 556   | 1e6 / (1 + 0.004·200) = 1e6 / 1.8         |
| 550   | 500 000   | 1e6 / (1 + 0.004·250) = 1e6 / 2.0         |
| 600   | 454 545   | 1e6 / (1 + 0.004·300) = 1e6 / 2.2         |
| 700   | 384 615   | 1e6 / (1 + 0.004·400) = 1e6 / 2.6         |

#### Why Tabulated Interpolation Instead of ADParsedMaterial

The Docker image `idaholab/moose:latest` does not include the JIT compilation
toolchain (LLVM/Clang dev libraries) that `ADParsedMaterial` requires to
compile expression strings at runtime. Attempting to use `ADParsedMaterial`
in that environment produces a runtime error. `ADPiecewiseLinearInterpolationMaterial`
is a compiled C++ object that performs the same task — evaluating a scalar
function of one variable — without any JIT dependency. For monotone, smooth
functions sampled at sufficient resolution, the piecewise-linear approximation
is exact to plotting precision.

#### The AD Chain Through Tabulated Materials

`ADPiecewiseLinearInterpolationMaterial` is an AD-compatible object. When MOOSE
evaluates the residual of the V equation at a quadrature point:

```
R_V = integral( σ(T) · grad(V) · grad(φ) )  dΩ
```

the derivative of R_V with respect to T (the off-diagonal Jacobian block
dR_V/dT) involves dσ/dT. Because the interpolation object operates on AD dual
numbers, the derivative dσ/dT is propagated automatically via the slope of the
piecewise-linear segments. No hand-coded Jacobian is needed.

Similarly the T residual includes:

```
R_T = integral( (ρ·cp·T_dot - Q)·φ )  dΩ
   where Q = σ(T)·|grad(V)|²
```

The derivative dR_T/dT includes both the direct ρ·cp·dT_dot/dT term and the
indirect -dQ/dT = -dσ/dT·|grad(V)|² term, computed automatically by AD.

#### Why `SMP full = true` Matters More Here

With constant σ (Case 17) the off-diagonal Jacobian block dR_V/dT is zero
because V and T are coupled only through the material's T dependence (which
was absent). The block dR_T/dV is nonzero (Joule source depends on grad(V)).

With σ(T) (Case 28):

- dR_V/dT is nonzero — changing T changes σ(T), which changes the V equation.
- dR_T/dV is nonzero — changing V changes |grad(V)|², which changes Q.
- dR_T/dT is nonzero — changing T changes σ(T), which changes Q directly.

All three off-diagonal-relevant blocks are filled. `SMP full = true` ensures
all of them are included in the preconditioner, which is critical for Newton
convergence on this strongly coupled system.

### Block: `[Executioner]`

Transient NEWTON solver, dt = 0.25 s, end_time = 5.0 s (20 timesteps). With
the stronger nonlinearity from σ(T), Newton uses 2 iterations per timestep —
identical to Case 17 in practice because the temperature rise at t = 5 s is
modest (~34 K above T_ref) and the feedback is mild. The tolerances
(nl_rel_tol = 1e-8, nl_abs_tol = 1e-10) are sufficient.

### Block: `[Postprocessors]`

Same three postprocessors as Case 17: `max_T`, `avg_T`, `max_V`. These allow
a direct comparison of temperature histories between Case 17 and Case 28.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case28-twoway-joule-heating \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case28_twoway_joule_heating.i 2>&1 | tail -30'
```

The run produces:

- `case28_twoway_joule_heating_out.e` — Exodus file with V and T fields at all
  20 timesteps plus the initial condition frame
- `case28_twoway_joule_heating_out.csv` — CSV with max_T, avg_T, max_V vs time

Expected runtime: a few seconds on a single core. Newton converges in
2 iterations every timestep.

---

## Expected Results

### Temperature Rise (Comparison with Case 17)

The most important observable difference from Case 17 is a **slower temperature
rise** with a **lower peak temperature**. The negative feedback from σ(T)
reduces the Joule source as temperature climbs.

At t = 5 s, measured values:

| Postprocessor | Case 28 (σ(T)) | Notes                                       |
|---------------|----------------|---------------------------------------------|
| max_T         | 334.3 K        | At the centre of the domain (x = 1 m)       |
| avg_T         | rises linearly | Steady linear increase across all timesteps |
| max_V         | ~9.95 V        | Matches applied BC of 10 V as expected      |

The max_T of 334.3 K represents a 34.3 K rise above the 300 K electrode
boundary condition. With the metallic α = 0.004 /K, the conductivity at the
hottest point has dropped to σ₀ / (1 + 0.004 * 34.3) ≈ 0.88 σ₀ — a 12%
reduction. This is the negative feedback at work: as the centre heats, its
conductivity drops, reducing the local Joule source and slowing further
temperature rise.

The max_V of ~9.95 V (rather than exactly 10.0 V) reflects that the
`ElementExtremeValue` postprocessor samples at element quadrature points, which
do not coincide exactly with the Dirichlet boundary nodes. The applied boundary
condition is exactly 10 V.

### avg_T Behaviour

The domain-averaged temperature rises linearly with time. This is expected:
the total Joule power input to the domain changes only slightly as σ(T) evolves
(the conductivity changes are modest over a 34 K range), so the volume-averaged
heating rate is approximately constant, giving the observed linear rise.

### Non-Uniform V Distribution

In Case 17 the V field is a perfect linear ramp V(x) = 10(1 - x/2), identical
at all times and independent of y. In Case 28 the V field develops a slight
spatial non-uniformity driven by the temperature distribution:

- The central region (x ≈ 1 m) becomes hotter and therefore less conductive.
- The equipotential lines V = const, which were straight vertical lines in
  Case 17, bow slightly around the high-resistance hot spot.
- The current density J = -σ(T)·grad(V) is no longer perfectly horizontal —
  it acquires small y-components as current routes around the more resistive
  centre.

This effect is small at t = 5 s with these parameters but would be clearly
visible in a long-duration run or with a larger α.

### Joule Heating Distribution

In Case 17 Q = σ·|grad(V)|² is spatially uniform (constant σ and 1D V field).
In Case 28 Q decreases at the hot centre: the drop in σ(T) outweighs any
concentration of grad(V). The result is a Joule source that is slightly
lower at the midplane and slightly higher near the electrodes where T is
clamped at 300 K and σ remains near σ₀.

This redistribution further slows the central temperature rise and reinforces
the negative feedback.

---

## Key Takeaways

### ADPiecewiseLinearInterpolationMaterial for Temperature-Dependent Properties

`ADPiecewiseLinearInterpolationMaterial` is the Docker-compatible MOOSE tool
for defining material properties that depend on other field variables through
tabulated data. The `variable` parameter names the field variable (here T), and
the `x`/`y` arrays define the lookup table. The object is fully AD-compatible:
the slope of each linear segment provides dσ/dT automatically to the Jacobian
assembly.

This pattern applies directly to any temperature-dependent property expressed
as tabulated data: thermal conductivity k(T), viscosity μ(T), yield stress
σ_y(T), diffusivity D(T), etc. The only changes are the property name,
variable name, and table values.

When `ADParsedMaterial` is available (a full MOOSE build with JIT support),
the same physics can be expressed with the analytic formula
`'${sigma0} / (1.0 + ${alpha} * (T - ${T_ref}))'` and the result is
mathematically identical for smooth, well-sampled functions.

### Two-Way Coupling Between V and T

One-way coupling (Case 17): V solves independently; T is driven by Q(V).
Two-way coupling (Case 28): V and T are solved simultaneously, each depending
on the other. The Jacobian has nonzero blocks in all four V-T quadrants.
MOOSE handles this automatically when both variables are in the same
`[Variables]` block and the coupling flows through AD-enabled material objects.

### Negative vs. Positive Feedback

The sign of dσ/dT determines the character of the coupled system:

- **Metallic (dσ/dT < 0)**: Higher T lowers σ, reduces Q, damps T rise.
  Negative feedback — the system is self-stabilising. This is the case here.

- **Semiconducting (dσ/dT > 0)**: Higher T raises σ, increases Q, drives T
  higher still. Positive feedback — can lead to thermal runaway. The system
  may have no stable steady state if the heat removal is insufficient.

The tabulated values used here faithfully capture the metallic (stable) regime.
Reversing the table to give increasing σ with T would simulate semiconductor-
like behaviour and would require either tighter nonlinear tolerances or an
adaptive timestep to remain stable.

### AD Chain Through Tabulated Materials

The MOOSE AD framework propagates derivatives through piecewise-linear
interpolations automatically. When `ADPiecewiseLinearInterpolationMaterial`
computes σ(T), the dual-number arithmetic encodes dσ/dT = slope of the
enclosing linear segment alongside the value. Every subsequent object that
reads `electrical_conductivity` — the `ADHeatConduction` kernel for V, and the
`ElectromagneticHeatingMaterial` — inherits this derivative information,
enabling exact Jacobian assembly for the full coupled system without any
hand-coded partial derivatives.

---

## Experiments to Try

### Experiment 1: Compare max_T History with Case 17

Run both Case 17 and Case 28 and plot max_T vs. time from their respective
CSV files. The Case 28 curve should lie below Case 17 at every timestep,
with the gap growing as temperatures (and therefore the σ reduction) increase.

At t = 5 s the measured max_T for Case 28 is 334.3 K. Compare this to the
Case 17 value to quantify the effect of the negative feedback.

### Experiment 2: Increase Voltage to Amplify the Effect

Change `V_left` from 10 V to 50 V. The Joule source scales as V², so Q
increases 25-fold. The temperature rises much faster, the σ reduction becomes
large (e.g., at 100 K above T_ref, σ drops to σ₀/1.4), and the feedback
effect is dramatic. The gap between constant-σ and σ(T) predictions will be
very large at t = 5 s. Extend the table if temperatures are expected to exceed
700 K.

### Experiment 3: Simulate Semiconductor Positive Feedback

Reverse the table so that σ increases with T:

```
x = '250  300   350     400     450     500     550     600     700'
y = '750000  1e6  1250000  1500000  1750000  2e6  2250000  2500000  3e6'
```

This simulates a material whose conductivity increases with temperature.
With a high enough voltage, the Joule heating will accelerate faster than heat
can be removed, and the simulation may fail to converge — an indicator of
thermal runaway in the physics. Try reducing dt (e.g., dt = 0.05) and
nl_max_its = 50 to follow the instability further.

### Experiment 4: Vary the Slope (α Equivalent) to See Feedback Strength

Recompute the table with different α values (0.001, 0.004, 0.010, 0.020 /K)
and run a parameter study. Plot max_T at t = 5 s vs. α. As α increases, the
negative feedback strengthens and max_T decreases. For very large α the
effective conductivity near the electrodes (where T = 300 K = T_ref, so
σ = σ₀) will differ greatly from the centre, producing a strongly non-uniform
V distribution and current crowding near the electrodes.

### Experiment 5: Add an AuxVariable to Track σ(T) Spatially

Add an AuxVariable and AuxKernel to visualise how σ varies across the domain:

```
[AuxVariables]
  [sigma_field]
    order  = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [sigma_aux]
    type             = ADMaterialRealAux
    variable         = sigma_field
    property         = electrical_conductivity
    execute_on       = timestep_end
  []
[]
```

In Paraview you can then animate σ(x,y,t) alongside T(x,y,t) and watch the
conductivity depression track the growing hot spot.

### Experiment 6: Increase Table Resolution

The current table uses 9 points spanning 250–700 K. Add more points near the
operating range (300–350 K) to reduce interpolation error:

```
x = '280  290  300  305  310  315  320  325  330  335  340  350  400  500  700'
```

Verify that the postprocessor results do not change significantly — this
confirms the 9-point table was already sufficiently resolved for the 34 K
temperature rise observed in this run.
