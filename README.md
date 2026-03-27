# SpeakEasy: A modular, open‑source loudspeaker simulation and parameter extraction toolkit.

## What can SpeakEasy do?
- Forward simulation: given T/S parameters (plus extended parameters like voice‑coil inductance and frequency‑dependent radiation), predict SPL, impedance, displacement, and efficiency.
- Inverse parameter extraction: from a simple electrical impedance measurement (soundcard + series resistor), recover all T/S parameters in one shot.

All code is released under the AGPLv3 license – you can use, modify, and share it freely, but any derivative work or service must also be open source.


## Prerequisites:

- MATLAB (R2018b or newer recommended)
- Python 3.8+ with the following packages:
  - numpy
  - scipy
  - matplotlib
  - sounddevice

Install them via pip:
pip install numpy scipy matplotlib sounddevice

- Soundcard with at least one line output and one line input (a standard USB audio interface works perfectly).
⚠️ Do not use a microphone‑only input, two input channels are needed. 


## Getting Started:

1. Clone the repository
   git clone https://github.com/theIvanR/speakeasy.git
   cd speakeasy

2. Read the white paper
   The detailed theory and derivation are documented in whitepaper.pdf.
   It is strongly recommended to understand the basics before using the tools.

3. Set up MATLAB path
   In MATLAB, navigate to the repository folder and add the subfolders to the path:
   addpath(genpath(pwd));


## Forward Problem:

Goal: Predict SPL, impedance, and other performance metrics from a set of T/S (plus) parameters.

1. Open MATLAB and run:
   forwards_ts_to_plots

2. A dialog will appear – select the driver parameter file or enter parameters manually.

3. Choose simulation options:
   - Configuration: single driver, isobaric series, parallel, or series‑parallel.
   - Enclosure: infinite baffle or sealed box.
   - Radiation model: constant, bessel, or bem.

4. The script will produce:
   - SPL (dB) vs frequency
   - Cone displacement
   - Electrical impedance magnitude & phase
   - Efficiency and power plots

All results are also saved to the results/ folder.


## Inverse Problem:

Goal: Determine T/S parameters from a real driver using only an electrical impedance measurement.

Step 1: Measure the impedance with bode_plotter.py

1. Connect the measurement circuit:
Soundcard line out --> series resistor (e.g., 10 Ω) --> driver --> soundcard line in (channel 2)
                 |
                 --> soundcard line in (channel 1)

2. Run the measurement script:
   python bode_plotter.py

Step 2: Fit parameters with inverse_speaker_to_ts.m

1. In MATLAB, run:
   inverse_speaker_to_ts

2. Select the CSV file from Step 1.

3. Enter an initial guess for the parameters.

4. The Levenberg‑Marquardt optimizer will fit the model to your measured impedance, producing a complete set of T/S parameters.

5. Results are displayed and saved as a JSON file for later use in the forward problem.


## Contributing:

Contributions are welcome! Please open an issue or submit a pull request.


Enjoy!
