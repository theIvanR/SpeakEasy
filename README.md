# Speakeasy 🔈

**AGPLv3 powerhouse for electroacoustic modelling – scalable, physics‑first, and open.**

[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-orange)](https://www.mathworks.com)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)

Speakeasy is a **research‑grade toolbox** for loudspeaker modelling, measurement, and parameter inversion.  
It combines a calibrated soundcard‑based impedance analyser, a full‑wave Boundary Element Method (BEM) radiation model, and a **Bayesian inverse solver** that gives you confidence intervals on Thiele‑Small parameters – all in a modular, AGPLv3‑licensed codebase.



## ✨ What can it do?

### 1. Calibrated soundcard Bode plotter (Python) – *release*
Turn any PC with a stereo line‑in/out into an **accurate impedance measurement rig**.

- 3‑point calibration (channel mismatch, open, short) stored in a `.cal` file.
- Coherence‑weighted estimation of transfer function.
- Automatically saves `.csv` with magnitude, phase, and coherence.
- **Example:** measure a loudspeaker’s electrical impedance from 20 Hz to 20 kHz with a simple voltage divider (one sense resistor).

```yaml
                                                          --> soundcard line in (channel 2)
                                                         |
Soundcard line out --+--> series resistor (e.g., 22 Ω) --+--> driver --> gnd
                     |
                      --> soundcard line in (channel 1)
```

### 2. Forward Thiele‑Small solver (MATLAB) – *release*
Predict **SPL, impedance, cone displacement, and efficiency** from a complete set of driver parameters.

- Arbitrary baffle geometry via BEM (circle, square, or any `.stl` mesh).
- Enclosure models: infinite baffle, sealed box (horn planned soon, ported & arbitrary wave eq enclosure planned later).
- Outputs: electrical impedance `Ze`, on‑axis pressure, velocity, and true acoustic efficiency.
- **Example:** design a sealed‑box subwoofer and visualise its displacement‑limited power handling.

### 3. Inverse problem solver (MATLAB) – *experimental*
Recover **all T/S parameters from a *single* electrical impedance measurement and their confidence intervals**.

- Two‑pass Bayesian MAP fit with **SVD‑shaped priors** that handle the notorious `Rms` ⇔ `Cms` ambiguity.
- Automatically estimates `fs0` from the measured phase zero‑crossing.
- Outputs: fitted `Re, Le, Bl, Rms, Mms, Cms` plus derived `Qms, Qes, Qts, Vas` and their uncertainties.
- **Example:** measure a raw driver’s impedance, run `main_inverse.m`, and get a complete parameter set with honest error bars.

**NOTE on identifiability**
Let F: P → O be the forward map from parameter space P ⊂ ℝ⁶ (e.g., Re, Le, Bl, Rms, Mms, Cms) to the observed electrical impedance Ze(ω) ∈ O, where O is a function space. For a single measurement, the preimage F⁻¹({Ze_obs}) is not a singleton – it is a submanifold of P of positive dimension. Equivalently, the intersection of the measurement manifold M = F(P) with the observed data is a set of positive measure in parameter space.

The solver finds the intersection of this measurement manifold with the prior-informed feasible region. The problem is ill‑posed in the sense of Hadamard: solutions are not unique, and small perturbations can lead to large parameter variations.

To recover a unique solution (i.e., to collapse the feasible set to the empty set of alternative parameter combinations), we lift the observation space by adding independent observables. Each new observable – e.g., added‑mass impedance, near‑field pressure, or laser velocimetry – defines an additional forward map F_j: P → O_j. The joint forward map (F, F₁, …, F_k): P → O × O₁ × … × O_k has an associated adjoint operator that, under ideal conditions (sufficiently many informative measurements), reduces the preimage to a single point. In Bayesian terms, the posterior covariance shrinks, and the confidence intervals tighten.

For practical use: A single impedance measurement yields a ballpark estimate suitable for most engineering tasks (enclosure design, crossover simulation). To eliminate the remaining ambiguity (e.g., for research or high‑precision driver characterisation), the framework scales seamlessly to multiple observables – just append their residuals, and the same Bayesian machinery computes the joint posterior.

## 🔬 What makes it special?

This is **not** yet another T/S calculator. Speakeasy was built on a strict physical contract:

 This solver operates under the following physical assumptions:

 1. ACOUSTIC MEDIUM
    - Linear, inviscid, isentropic fluid (Helmholtz regime)
    - Constant sound speed c and density rho
    - No thermal or viscous boundary layer losses

 2. RADIATION MODEL
    - Surface Boundary Element Method (BEM)
    - Arbitrary geometry supported via discretized mesh:
        geom = {xc, A, aeq, Np}
    - Radiation is fully geometry-dependent

 3. STRUCTURAL MODEL
    - Single-degree-of-freedom rigid piston (current version)
    - No cone breakup modes included (future extension planned)

 4. ENCLOSURE MODEL
    - Lumped mechanical impedance ONLY
    - IB (infinite baffle) = Zbox = 0
    - Sealed box = compliance + optional damping
    - No spatial cavity modes included

 5. COUPLING
    - Electromechanical coupling via Bl product
    - Acoustic loading enters as mechanical impedance

 6. IMPORTANT DESIGN RULE
    - Geometry affects ONLY radiation operators
    - Enclosure affects ONLY lumped mechanical impedance
    - No hidden coupling paths are allowed

## Contributing:

Contributions are welcome! Please open an issue or submit a pull request.


Enjoy!
