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
Recover **all T/S parameters from a single electrical impedance measurement** – with **95% confidence intervals**.

- Two‑pass Bayesian MAP fit with **SVD‑shaped priors** that handle the notorious `Rms` ⇔ `Cms` ambiguity.
- Automatically estimates `fs0` from the measured phase zero‑crossing.
- Outputs: fitted `Re, Le, Bl, Rms, Mms, Cms` plus derived `Qms, Qes, Qts, Vas` and their uncertainties.
- **Example:** measure a raw driver’s impedance, run `main_inverse.m`, and get a complete parameter set with honest error bars.



## 🔬 What makes it special?

This is **not** yet another T/S calculator. Speakeasy was built on a strict physical contract:

1. ACOUSTIC MEDIUM – Linear, inviscid, isentropic (Helmholtz)
2. RADIATION MODEL – Surface BEM, arbitrary geometry
3. STRUCTURAL MODEL – Rigid piston (SDOF) – extensible
4. ENCLOSURE MODEL – Lumped mechanical impedance only
5. COUPLING – Bl product, acoustic loading as mechanical impedance
6. IMPORTANT RULE:
      - Geometry → radiation operators ONLY
      - Enclosure → lumped mechanical impedance ONLY
      - No hidden coupling paths


## Contributing:

Contributions are welcome! Please open an issue or submit a pull request.


Enjoy!
