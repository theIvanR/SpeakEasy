import numpy as np
import pandas as pd
import sounddevice as sd
from scipy import signal as sig
import matplotlib.pyplot as plt
from typing import Optional, Tuple
import warnings

class BodePlotter:
    """
    A Bode plotter using a sound card to measure impedance via a voltage divider.
    """
    def __init__(
        self,
        input_device: int,
        output_device: int,
        fs: int = 44100,
        duration_s: float = 10.0,
        output_rms: float = 0.1,
        r_sense: float = 100.0,
        nperseg: int = 16384,
        noverlap: Optional[int] = None,
        f_min: float = 20.0,
        f_max: float = 20000.0,
    ):
        """
        Parameters
        ----------
        input_device : int
            SoundDevice index for input.
        output_device : int
            SoundDevice index for output.
        fs : int
            Sampling rate in Hz.
        duration_s : float
            Duration of stimulus in seconds.
        output_rms : float
            Desired RMS level of the output stimulus.
        r_sense : float
            Sense resistor value in Ohms.
        nperseg : int
            Length of each segment for Welch's method.
        noverlap : int or None
            Number of points to overlap. If None, defaults to nperseg // 2.
        coh_min : float
            Minimum coherence to consider a frequency bin valid.
        f_min : float
            Lower bound of frequency range to retain in results (Hz).
        f_max : float
            Upper bound of frequency range to retain in results (Hz).
        """
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

    def save_to_csv(self, filename: str = "impedance_bode.csv"):
        """Save results to a CSV file."""
        if self.Z is None:
            raise RuntimeError("No impedance data available. Call run() first.")

        df = pd.DataFrame({
            'frequency_hz': self.f,
            'Z_mag_ohm': np.abs(self.Z),
            'Z_phase_deg': np.rad2deg(np.angle(self.Z)),
            'Coherence': self.coherence
        })

        df.to_csv(filename, index=False)
        print(f"Saved results to {filename}")

    def plot(self):
        """Create a Bode plot of impedance magnitude and phase, plus coherence."""
        if self.Z is None:
            raise RuntimeError("No impedance data available. Call run() first.")

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        # Magnitude
        ax1.semilogx(self.f, np.abs(self.Z))
        ax1.set_ylabel('Impedance (Ω)')
        ax1.grid(True, which='both', linestyle=':')
        ax1.set_title('Measured Impedance')

        # Phase
        ax2.semilogx(self.f, np.rad2deg(np.angle(self.Z)))
        ax2.set_ylabel('Phase (deg)')
        ax2.grid(True, which='both', linestyle=':')

        # Coherence
        ax3.semilogx(self.f, self.coherence)
        ax3.set_ylabel('Coherence')
        ax3.set_xlabel('Frequency (Hz)')
        ax3.set_ylim(0, 1.05)
        ax3.grid(True, which='both', linestyle=':')

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # Example usage – adapt device indices to your system
    INPUT_DEVICE = 1   # Realtek Line In
    OUTPUT_DEVICE = 4  # Sharkoon DAC (WASAPI)

    plotter = BodePlotter(
        input_device=INPUT_DEVICE,
        output_device=OUTPUT_DEVICE,
        fs=44100,
        duration_s=90.0,
        output_rms=0.1,
        r_sense=22.060,
        nperseg=16384,
        f_min=20.0,
        f_max=20000.0,
    )

    input(
        "\nWire it as:\n"
        "  OUT -> Shunt Resistor -> node A -> DUT -> GND\n"
        "  L_in measures node before the shunt resistor\n"
        "  R_in measures node A (DUT node)\n"
        "Press Enter to start..."
    )

    plotter.run()
    plotter.save_to_csv("impedance_bode.csv")
    plotter.plot()