import numpy as np
import sounddevice as sd
import scipy.signal as sig

# =========================
# Configuration
# =========================
INPUT_DEVICE = 1   # Line In
OUTPUT_DEVICE = 4   # DAC

FS = 44100
DURATION_S = 5.0
OUTPUT_RMS = 0.05

CHANNEL_L = 0
CHANNEL_R = 1

# =========================
# Helpers
# =========================
def rms(x):
    x = np.asarray(x)
    return np.sqrt(np.mean(x**2))

def make_white_noise(n, rms_target=0.1, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(n)
    x = x / (rms(x) + 1e-15) * rms_target
    return x

# =========================
# Setup devices
# =========================
print("Available devices:\n")
print(sd.query_devices())

sd.default.device = (INPUT_DEVICE, OUTPUT_DEVICE)
sd.default.samplerate = FS

# =========================
# Stimulus
# =========================
n = int(DURATION_S * FS)
x = make_white_noise(n, rms_target=OUTPUT_RMS)

input("\nLoop L_out -> splitter -> L_in and R_in.\nPress Enter to start...")

# =========================
# Acquisition
# =========================
rec = sd.playrec(x.reshape(-1, 1), samplerate=FS, channels=2, dtype='float64')
sd.wait()

xL = rec[:, CHANNEL_L]
xR = rec[:, CHANNEL_R]

# =========================
# DC removal
# =========================
xL = xL - np.mean(xL)
xR = xR - np.mean(xR)

# =========================
# Gain estimation (RMS-based)
# =========================
gain_L = rms(xL)
gain_R = rms(xR)

gain_ratio = gain_R / (gain_L + 1e-15)

print("\n=== Gain Estimation ===")
print(f"RMS L_in: {gain_L:.6f}")
print(f"RMS R_in: {gain_R:.6f}")
print(f"Gain ratio (R / L): {gain_ratio:.6f}")

# =========================
# Delay estimation (cross-correlation)
# =========================
corr = sig.correlate(xR, xL, mode='full')
lags = sig.correlation_lags(len(xR), len(xL), mode='full')

lag = lags[np.argmax(corr)]
delay_seconds = lag / FS

print("\n=== Time Delay Estimation ===")
print(f"Sample lag (R relative to L): {lag} samples")
print(f"Time delay: {delay_seconds * 1e6:.3f} µs")

# =========================
# Optional: sanity check plot
# =========================
import matplotlib.pyplot as plt

plt.figure()
plt.plot(lags, corr)
plt.title("Cross-correlation")
plt.xlabel("Lag (samples)")
plt.ylabel("Correlation")
plt.grid(True)
plt.show()