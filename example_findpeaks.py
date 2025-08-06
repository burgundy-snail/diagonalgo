import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

signal = np.array([1,2,3,4,5,6,4,3,2,2,1])
peaks = find_peaks(signal)
print(peaks)

x = np.sin(2*np.pi*(2**np.linspace(2,10,1000))*np.arange(1000)/48000) + np.random.normal(0, 1, 1000) * 0.15

#play around with these parameters.
peaks, _ = find_peaks(x, distance=35)
peaks2, _ = find_peaks(x, prominence=0.65)      # BEST!
peaks3, _ = find_peaks(x, width=8)
peaks4, _ = find_peaks(x, threshold=0.4)     # Required vertical distance to its direct neighbouring samples, pretty useless
plt.subplot(2, 2, 1)
plt.plot(peaks, x[peaks], "xr"); plt.plot(x); plt.legend(['distance'])
plt.subplot(2, 2, 2)
plt.plot(peaks2, x[peaks2], "ob"); plt.plot(x); plt.legend(['prominence'])
plt.subplot(2, 2, 3)
plt.plot(peaks3, x[peaks3], "vg"); plt.plot(x); plt.legend(['width'])
plt.subplot(2, 2, 4)
plt.plot(peaks4, x[peaks4], "xk"); plt.plot(x); plt.legend(['threshold'])
plt.show()

# -----PEAKDETECT METHOD-----

# Import library
from findpeaks import findpeaks
# Initialize with appropriate lookahead for large dataset
fp = findpeaks(method='peakdetect', lookahead=200, interpolate=None)

# Example 1d-vector with noise
i = 10000
xs = np.linspace(0,3.7*np.pi,i)
X = (0.3*np.sin(xs) + np.sin(1.3 * xs) + 0.9 * np.sin(4.2 * xs) + 0.06 * np.random.randn(i))

# Fit peakdetect method on the noisy 1d-vector
results = fp.fit(X)
# Plot
fp.plot()