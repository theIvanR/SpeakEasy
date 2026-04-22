import numpy as np
import pandas as pd
import sounddevice as sd
from scipy import signal as sig
import matplotlib.pyplot as plt
from typing import Optional, Tuple
import warnings

"""
A Bode plotter using a sound card to measure impedance via a voltage divider.

Wire it as:
    OUT -> Shunt Resistor -> node A -> DUT -> GND
    L_in measures node before the shunt resistor
    R_in measures node A (DUT node)
      
Adjust: 
    Output_rms as needed, but dont go too high, start with 0.2 or so.
"""

class BodePlotter:
    
    def __init__(
        self,
        input_device: int,      # SoundDevice input index
        output_device: int,     # SoundDevice output index
        fs: int,                # Sampling rate (Hz)
        duration_s: float,     # Stimulus duration (s)
        output_rms: float,     # Output RMS level (sigma!)
        r_sense: float,        # Sense resistor (Ohms)
        nperseg: int,          # Welch segment length
        f_min: float,          # Lower frequency bound (Hz)
        f_max: float,          # Upper frequency bound (Hz)
        noverlap: Optional[int] = None,  # Overlap (None → 50%)
    ):
        
        self.input_device = input_device
        self.output_device = output_device
        self.fs = fs
        self.duration_s = duration_s
        self.output_rms = output_rms
        self.r_sense = r_sense
        self.nperseg = nperseg
        self.noverlap = noverlap if noverlap is not None else nperseg // 2
        self.f_min = f_min
        self.f_max = f_max

        # Internal state
        self.stimulus = None
        self.recording = None
        self.f = None
        self.H_div = None
        self.coherence = None
        self.Z = None

    def _check_devices(self):
        """Verify that selected devices exist and support required channels."""
        try:
            sd.query_devices(self.input_device)
            sd.query_devices(self.output_device)
        except sd.PortAudioError as e:
            raise RuntimeError(f"Invalid device index: {e}")

        # Get device info
        in_dev = sd.query_devices(self.input_device, 'input')
        out_dev = sd.query_devices(self.output_device, 'output')
        if in_dev['max_input_channels'] < 2:
            raise ValueError(f"Input device {self.input_device} has only {in_dev['max_input_channels']} input channels; need at least 2.")
        if out_dev['max_output_channels'] < 1:
            raise ValueError(f"Output device {self.output_device} has no output channels.")

    def _generate_stimulus(self, seed: int = 12345, fade_s: float = 0.05) -> np.ndarray:
        """
        Generate band‑limited white noise stimulus.
        The noise is filtered to the range [f_min, f_max] to concentrate energy
        in the audio band.
        """
        n = int(self.duration_s * self.fs)
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)

        # Set RMS
        x = x / (np.sqrt(np.mean(x**2)) + 1e-15) * self.output_rms

        # Apply fade in/out
        fade_n = int(fade_s * self.fs)
        if fade_n > 0 and 2 * fade_n <= n:
            fade = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, fade_n))
            x[:fade_n] *= fade
            x[-fade_n:] *= fade[::-1]

        # Band‑limit to [f_min, f_max]
        # Design a bandpass filter using a Butterworth
        nyquist = 0.5 * self.fs
        low = self.f_min / nyquist
        high = self.f_max / nyquist
        # Ensure bounds are within (0,1)
        low = max(0.001, min(0.999, low))
        high = max(0.001, min(0.999, high))
        b, a = sig.butter(4, [low, high], btype='band')
        x = sig.filtfilt(b, a, x)

        return x

    def _record(self):
        """Play stimulus and record two input channels."""
        # Set default device and samplerate
        sd.default.device = (self.input_device, self.output_device)
        sd.default.samplerate = self.fs

        # Generate stimulus
        self.stimulus = self._generate_stimulus()

        # Record
        rec = sd.playrec(self.stimulus.reshape(-1, 1), samplerate=self.fs, channels=2, dtype='float64')
        sd.wait()

        self.recording = rec

    def _estimate_transfer(self):
        """Estimate transfer function H = V_dut / V_src and coherence."""
        if self.recording is None:
            raise RuntimeError("No recording available. Call _record() first.")

        v_src = self.recording[:, 0]   # source side
        v_dut = self.recording[:, 1]   # DUT side

        # Remove DC offset
        v_src = v_src - np.mean(v_src)
        v_dut = v_dut - np.mean(v_dut)

        # Compute cross‑spectral density and power spectra
        f, P_yx = sig.csd(v_dut, v_src,
                          fs=self.fs,
                          window='hann',
                          nperseg=self.nperseg,
                          noverlap=self.noverlap,
                          detrend=False,
                          scaling='density')

        _, P_xx = sig.welch(v_src,
                            fs=self.fs,
                            window='hann',
                            nperseg=self.nperseg,
                            noverlap=self.noverlap,
                            detrend=False,
                            scaling='density')

        _, coh = sig.coherence(v_dut, v_src,
                               fs=self.fs,
                               window='hann',
                               nperseg=self.nperseg,
                               noverlap=self.noverlap,
                               detrend=False)

        # Transfer function (FIXED as the csd does the conjugate!)
        #H = P_yx / (P_xx + 1e-30) WRONG!!! (FLIPS PHASE)
        H = np.conj(P_yx) / (P_xx + 1e-30)

        # Limit frequency range
        mask = (f >= self.f_min) & (f <= self.f_max)
        self.f = f[mask]
        self.H_div = H[mask]
        self.coherence = coh[mask]

    def _compute_impedance(self):
        """Compute impedance from voltage divider equation."""
        if self.H_div is None:
            raise RuntimeError("No transfer function available. Call _estimate_transfer() first.")
    
        den = 1.0 - self.H_div
    
        # Compute Z for all bins (no masking)
        Z = self.r_sense * self.H_div / den
    
        self.Z = Z

    def run(self):
        """Run the whole measurement process."""
        self._check_devices()
        self._record()
        self._estimate_transfer()
        self._compute_impedance()

    def save(self, filename: str = "impedance_bode", save_csv: bool = True, save_plot: bool = True):
        """Save results as CSV and/or plot."""
    
        if self.Z is None:
            raise RuntimeError("No impedance data available. Call run() first.")
    
        # --- Prepare data ---
        freq = self.f
        Z_mag = np.abs(self.Z)
        Z_phase = np.rad2deg(np.angle(self.Z))
        coh = self.coherence
    
        # --- Save CSV ---
        if save_csv:
            df = pd.DataFrame({
                'frequency_hz': freq,
                'Z_mag_ohm': Z_mag,
                'Z_phase_deg': Z_phase,
                'Coherence': coh
            })
            csv_name = f"{filename}.csv"
            df.to_csv(csv_name, index=False)
            print(f"Saved CSV → {csv_name}")
    
        # --- Save Plot ---
        if save_plot:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        
            # Magnitude
            ax1.semilogx(freq, np.abs(self.Z))
            ax1.set_ylabel('Impedance (Ω)')
            ax1.grid(True, which='both', linestyle=':')
            ax1.set_title('Measured Impedance')
        
            # Phase
            ax2.semilogx(freq, np.rad2deg(np.angle(self.Z)))
            ax2.set_ylabel('Phase (deg)')
            ax2.grid(True, which='both', linestyle=':')
        
            # Coherence
            ax3.semilogx(freq, self.coherence)
            ax3.set_ylabel('Coherence')
            ax3.set_xlabel('Frequency (Hz)')
            ax3.set_ylim(0, 1.05)
            ax3.grid(True, which='both', linestyle=':')
        
            plt.tight_layout()
        
            plot_name = f"{filename}.png"
            plt.savefig(plot_name, dpi=300, bbox_inches='tight')
            plt.close(fig)

        print(f"Saved Plot → {plot_name}")

if __name__ == "__main__":
    # Example usage – adapt device indices to your system
    INPUT_DEVICE = 1   # Realtek Line In
    OUTPUT_DEVICE = 6  # Realtek Line Out

    plotter = BodePlotter(
        input_device=INPUT_DEVICE,
        output_device=OUTPUT_DEVICE,
        fs=48000,
        duration_s=30.0,
        output_rms=0.2,
        r_sense=22.060,
        nperseg=32768,
        f_min=20.0,
        f_max=20000.0,
    )

    plotter.run()
    plotter.save("impedance_bode_plot")