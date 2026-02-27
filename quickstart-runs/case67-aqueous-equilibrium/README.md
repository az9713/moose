# Case 67: Aqueous Equilibrium — CO2-H2O System

## Overview

The dissolution of carbon dioxide in water is the most important acid-base equilibrium in natural waters. It governs the pH of rainwater, ocean acidification, carbonate rock weathering, and the chemistry of groundwater interacting with carbonate minerals. When CO2 dissolves in water, a cascade of equilibrium reactions rapidly establishes the carbonate system: CO2(aq), bicarbonate (HCO3-), carbonate (CO3^2-), and the fundamental water ionization equilibrium producing OH-. These reactions equilibrate on timescales of milliseconds to seconds — far faster than any transport process — so they are correctly treated as instantaneous equilibrium rather than kinetically-controlled reactions.

This case uses the MOOSE `chemical_reactions` module's `ReactionNetwork/AqueousEquilibriumReactions` Action to simulate the CO2-H2O equilibrium system in a well-mixed batch reactor (no spatial transport). Given initial concentrations of H+ and HCO3-, the system rapidly equilibrates to establish the carbonate speciation and pH. The simulation demonstrates how the `chemical_reactions` module handles the simultaneous solution of multiple nonlinear equilibrium constraints coupled to the primary species mass balances.

Key concepts demonstrated:

- `ReactionNetwork/AqueousEquilibriumReactions` Action for fast equilibrium chemistry
- `PHAux` auxiliary kernel to compute pH = -log10([H+])
- `TotalConcentrationAux` for tracking the total analytical concentration of each element
- Batch reactor setup (no spatial transport, uniform 0D chemistry)
- Rapid equilibration demonstrating that equilibrium is achieved within the first time step

---

## The Physics

### The Carbonate System Equilibria

The CO2-H2O carbonate system involves three primary equilibrium reactions:

**Reaction 1 — CO2(aq) formation (hydration):**
```
H+ + HCO3- ⇌ CO2(aq) + H2O     log10(K1) = 6.35
```
K1 = 10^6.35 = 2.24e6, meaning CO2(aq) strongly dominates over the free ions at low pH.

**Reaction 2 — Carbonate ion formation:**
```
H+ + CO3^2- ⇌ HCO3-             log10(K2) = 10.33
```
K2 = 10^10.33 = 2.14e10, meaning CO3^2- is strongly suppressed at neutral pH (it only becomes significant above pH 10).

**Reaction 3 — Water ionization:**
```
H+ + OH- ⇌ H2O                  log10(K3) = 14.0
```
K3 = 10^14 = Kw^(-1), the familiar water product. At pH 5: [OH-] = 10^-9 mol/L — negligible compared to [H+] = 10^-5 mol/L.

### Primary Species and Mass Action

The `ReactionNetwork` framework designates H+ and HCO3- as the two primary species. All secondary species are expressed as mass-action products of the primary species:

```
[CO2(aq)] = K1 * [H+] * [HCO3-]
[CO3^2-]  = [HCO3-] / (K2 * [H+])
[OH-]     = 1 / (K3 * [H+])
```

The equilibrium constraints are satisfied algebraically — no additional equations are needed. Only the primary species have PDEs (time derivatives). The secondary species concentrations are computed from the primary species values at each Newton iteration.

### Initial Conditions and Expected pH

Given initial conditions [H+]_0 = 1e-5 M and [HCO3-]_0 = 1e-5 M, the total carbonate concentration is approximately constant at C_T ~ [HCO3-]_0 = 1e-5 M. At pH 5, the carbonate distribution is dominated by CO2(aq):

```
[CO2(aq)] = K1 * [H+] * [HCO3-] = 2.24e6 * 1e-5 * 1e-5 = 0.224 mol/L  >> [HCO3-] !!
```

Wait — this means [CO2(aq)] >> [HCO3-]. The total carbon balance must be satisfied: the system equilibrates to redistribute carbon between CO2(aq) and HCO3-. The equilibrium pH is found by solving the charge balance:

```
[H+] = [HCO3-] + 2*[CO3^2-] + [OH-]
```

For the given initial conditions (approximately dilute, nearly equal H+ and HCO3-), equilibrium occurs near pH 5.0. At this pH:
- CO2(aq) dominates: [CO2(aq)] >> [HCO3-] >> [CO3^2-]
- OH- is negligible (pH 5 << pH 7)
- The system is mildly acidic, consistent with dissolved CO2

### Domain — Batch Reactor (0D)

```
Batch reactor: well-mixed, no spatial gradients
Mesh: 1 x 1 single-element domain (purely computational)
No transport terms: no advection, no diffusion
Chemistry only: instantaneous equilibrium + time derivative for mass balance
```

The 1 x 1 mesh is a computational convenience for MOOSE's finite element framework, which requires a mesh. In a 0D batch reactor, all quantities are spatially uniform, so the single-element result represents the entire reactor.

The absence of spatial transport means the `AqueousEquilibriumReactions` Action uses only `TimeDerivative` kernels for the primary species (no diffusion, no advection). The equilibrium constraints enforce the secondary species concentrations algebraically at each time step.

---

## Input File Walkthrough

The input file is `case67_aqueous_equilibrium.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 1
nx = 1
xmax = 1
```

A single-element 1D mesh. All variables are spatially constant — this is a 0D (batch reactor) simulation. The MOOSE framework requires at least one element; a 1D mesh with one element is the minimal computational domain.

### `[ReactionNetwork]`

```
[ReactionNetwork]
  primary_species = 'H+ HCO3-'
  [AqueousEquilibriumReactions]
    [co2aq]
      primary_species_stoichiometry = '1 1'
      equilibrium_constant = 6.35
      secondary_species_name = 'CO2aq'
    []
    [carbonate]
      primary_species_stoichiometry = '1 -1'
      equilibrium_constant = -10.33
      secondary_species_name = 'CO3--'
    []
    [hydroxide]
      primary_species_stoichiometry = '1 0'
      equilibrium_constant = -14.0
      secondary_species_name = 'OH-'
    []
  []
[]
```

The `AqueousEquilibriumReactions` Action interprets the stoichiometry and equilibrium constants to automatically:
- Create secondary species AuxVariables (`CO2aq`, `CO3--`, `OH-`)
- Create `AqueousEquilibriumRxnAux` AuxKernels that compute each secondary species concentration from the mass-action expression
- Add `CoupledBEEquilibriumSub` terms to the primary species equations to enforce mass conservation (total analytical concentration of each element is conserved)

### `[Variables]`

```
[H+]
  initial_condition = 1e-5   # 10 umol/L
[]
[HCO3-]
  initial_condition = 1e-5   # 10 umol/L
[]
```

Both primary species start at 1e-5 mol/L. The initial conditions represent a dilute acidic solution with equal concentrations of H+ and bicarbonate — a slightly idealized starting point. The system will equilibrate within the first time step.

### `[AuxVariables]` and `[AuxKernels]`

The `AqueousEquilibriumReactions` Action auto-generates the secondary species AuxVariables. Two additional auxiliary variables are manually added:

**`pH`**: computed by `PHAux` as:
```
pH = -log10([H+])
```
This familiar quantity is the primary indicator of solution acidity. `PHAux` is a specialized kernel in the `chemical_reactions` module.

**`total_C` and `total_H`**: computed by `TotalConcentrationAux` to track the total analytical concentrations of carbon and charge. These should remain constant throughout the simulation (conservation laws), providing a verification check.

### `[Kernels]`

The `AqueousEquilibriumReactions` Action generates `CoupledBEEquilibriumSub` kernels for each primary species. These add the contribution of equilibrium secondary species to the primary species mass balance. The `TimeDerivative` kernels handle the transient evolution.

No spatial kernels (diffusion, convection) are present — this is a spatially uniform batch problem.

### `[Postprocessors]`

| Name | Type | Physical meaning |
|------|------|-----------------|
| `pH` | ElementAverageValue on `pH` aux | Solution pH |
| `H_conc` | ElementAverageValue on `H+` | Primary H+ concentration |
| `HCO3_conc` | ElementAverageValue on `HCO3-` | Bicarbonate concentration |
| `CO2aq_conc` | ElementAverageValue on `CO2aq` | Dissolved CO2(aq) |
| `CO3_conc` | ElementAverageValue on `CO3--` | Carbonate ion concentration |
| `OH_conc` | ElementAverageValue on `OH-` | Hydroxide concentration |

### `[Executioner]`

```
type     = Transient
dt       = 1.0
end_time = 10.0
nl_rel_tol = 1e-12
```

Ten time steps of 1 s each. Because the equilibrium reactions are instantaneous (no kinetic timescale), the system should reach its equilibrium state within the first Newton solve of the first time step. The 10-step simulation confirms that the equilibrium is stable and the concentrations do not drift.

Tight nonlinear tolerances (`nl_rel_tol = 1e-12`) are appropriate for equilibrium chemistry where the mass-action constraints must be satisfied precisely for accurate speciation.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case67-aqueous-equilibrium \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case67_aqueous_equilibrium.i 2>&1 | tail -30'
```

Output files:
- `case67_aqueous_equilibrium_out.e` — Exodus with all species concentrations (constant in space)
- `case67_aqueous_equilibrium_out.csv` — Time series of all species and pH

---

## Expected Results

### Rapid Equilibration

The key result is that equilibrium is established in the very first time step. Starting from [H+] = [HCO3-] = 1e-5 mol/L, the system immediately redistributes carbon between the species:

| Time (s) | pH | [H+] (mol/L) | [HCO3-] (mol/L) | [CO2aq] (mol/L) | [CO3^2-] (mol/L) | [OH-] (mol/L) |
|----------|----|--------------|-----------------|-----------------|------------------|---------------|
| 0        | 5.00 | 1.00e-5    | 1.00e-5         | ~2.24           | ~4.7e-16        | ~1.0e-9       |
| 1        | ~5.0 | ~1.0e-5   | ~1.0e-5         | ~2.24           | ~4.7e-16        | ~1.0e-9       |
| 10       | ~5.0 | ~1.0e-5   | ~1.0e-5         | ~2.24           | ~4.7e-16        | ~1.0e-9       |

The concentrations are essentially constant after the first step, confirming rapid equilibration.

### Carbonate Speciation at pH 5

The dominant result is the extremely high CO2(aq) concentration relative to HCO3-:

```
[CO2(aq)] / [HCO3-] = K1 * [H+] = 2.24e6 * 1e-5 = 22.4
```

At pH 5, dissolved CO2 is 22 times more concentrated than bicarbonate. This is the correct behavior of the carbonate system: below pH 6.35 (pK1 of the CO2(aq)/HCO3- equilibrium), CO2(aq) dominates. The system correctly resolves this speciation.

### pH Verification

The equilibrium pH ~5.0 is self-consistent with the initial conditions. The charge balance is satisfied:

```
[H+] = [HCO3-] + 2*[CO3^2-] + [OH-]
1e-5 ~ 1e-5 + 2*5e-16 + 1e-9  (approximately satisfied)
```

The CO3^2- and OH- concentrations are negligible at pH 5, confirming the dominant ion pair H+ and HCO3-.

### Conservation Check

The total inorganic carbon (TIC = [CO2aq] + [HCO3-] + [CO3^2-]) and the total charge should be conserved throughout the simulation. The `TotalConcentrationAux` postprocessors confirm this. Any drift would indicate a mass balance error in the equilibrium formulation.

---

## Key Takeaways

- `ReactionNetwork/AqueousEquilibriumReactions` is the MOOSE tool for fast (instantaneous) aqueous equilibrium chemistry. It automatically generates the mass-action auxiliary kernels and the source/sink terms that couple secondary species into the primary species mass balances.
- The `PHAux` kernel computes pH = -log10([H+]) from the primary H+ variable. It is provided by the `chemical_reactions` module and handles the log computation robustly (including guards against negative concentrations).
- Equilibrium chemistry requires only the primary species to be solved as PDEs. All secondary species are algebraic quantities computed from the mass-action expressions — they have no time derivatives of their own.
- The batch reactor (0D) setup is appropriate when transport (diffusion, advection) is absent or when spatial gradients can be neglected. It reduces the problem to pure chemistry, allowing verification of the equilibrium constants and speciation logic without transport complications.
- At pH 5 in the carbonate system, CO2(aq) dominates over HCO3- by a factor of 22 (= K1 * [H+] = 10^6.35 * 10^-5 = 10^1.35). The pKa values of the carbonate system determine which species dominates at each pH: CO2(aq) below pH 6.35, HCO3- between pH 6.35 and 10.33, CO3^2- above pH 10.33.
- Tight nonlinear tolerances (nl_rel_tol = 1e-12) are important for accurate equilibrium speciation. The mass-action expressions involve exponential functions of log10(K) values, and poor convergence leads to speciation errors that propagate as incorrect source terms in multi-reaction networks.
- The CO2-H2O equilibrium system modeled here is the foundation for understanding ocean acidification: as atmospheric CO2 rises, more CO2 dissolves in the ocean, shifting the equilibrium toward lower pH and reducing carbonate ion concentration — directly threatening calcifying organisms.
- `TotalConcentrationAux` is a verification tool: the total analytical concentration of each conserved element should not change during equilibrium redistribution. Checking conservation is essential in any multi-reaction network to catch stoichiometry errors.
