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
Recover **all T/S parameters from a single electrical impedance measurement and their confidence intervals**.

- Two‑pass Bayesian MAP fit with **SVD‑shaped priors** that handle the notorious `Rms` ⇔ `Cms` ambiguity.
- Automatically estimates `fs0` from the measured phase zero‑crossing.
- Outputs: fitted `Re, Le, Bl, Rms, Mms, Cms` plus derived `Qms, Qes, Qts, Vas` and their uncertainties.
- **Example:** measure a raw driver’s impedance, run `main_inverse.m`, and get a complete parameter set with honest error bars.

```yaml
[SPEAKER FIT] Pass 1 fs0 = 72.612 Hz | refined fs = 72.6028 Hz
[SPEAKER FIT] cond(J_data) = 7.830e+02 | cond(F_used) = 5.567e+01
[SPEAKER FIT] Fit complete. Exitflag = 3
====================================================
SMART FIT RESULTS
====================================================
--- Fitted Parameters ---
      c: 343
    rho: 1.2000
     Sd: 0.0140
     Re: 5.6804
     Le: 4.3705e-04
     Bl: 11.4295
    Rms: 15.0450
    Mms: 0.1273
    Cms: 3.7741e-05
     fs: 72.6028
    Qms: 3.8606
    Qes: 2.5257
    Qts: 1.5268

--- Derived ---
Fs (Hz): 72.602812
Qms    : 3.860614
Qes    : 2.525651
Qts    : 1.526802
Vas (m^3): 1.044341e-03
--- 95% Confidence Intervals (Parameters) ---
Re: [4.367, 7.389]
Le: [0.0003144, 0.0006076]
Bl: [8.678, 15.05]
Rms: [8.322, 27.2]
Mms: [0.09227, 0.1757]
Cms: [2.705e-05, 5.266e-05]
--- 95% Confidence Intervals (Derived) ---
     fs: [65.8056 80.1021]
    Cms: [2.7049e-05 5.2660e-05]
    Qms: [1.5818 9.4222]
    Qes: [1.0919 5.8423]
    Qts: [0.7138 3.2658]
    Vas: [7.4847e-04 0.0015]
    Rms: [8.3222 27.1988]

cond(J_data) = 7.830e+02
```



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
