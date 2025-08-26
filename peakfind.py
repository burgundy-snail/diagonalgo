import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks, savgol_filter, convolve
from scipy.ndimage import gaussian_filter1d as gaussian
from findpeaks import findpeaks

from diag import DiagCool

def debug_pd(sig, diagonal, lookahead=50, **kwargs):
    print(f"\n=== Debugging diagonal d={diagonal} ===")
    print(f"Signal length: {len(sig)}, min={sig.min()}, max={sig.max()}")

    # Run raw SciPy
    peaks, props = find_peaks(sig, distance=lookahead)
    print(f"SciPy returned {len(peaks)} raw peaks")

    # Peek at first few peaks
    if len(peaks) > 0:
        print("First peaks (x,y):", [(int(p), float(sig[p])) for p in peaks[:5]])

    # If your wrapper does more filtering, check before/after
    return peaks, props

def preprocess_signal(sig: np.ndarray, method: str = None, **kwargs) -> np.ndarray:
    """Optional preprocessing of 1D signal.
    
    Only accepts denoising-related kwargs - pass via 'denoise_kwargs' from peakfinder. Interpolation is handled by findpeaks.
    """
    if method is None:
        return sig
    
    # NaN handling
    arr = np.array(sig, dtype=float).copy()
    if np.isnan(arr).any():
        fill = np.nanmedian(arr)
        arr = np.where(np.isnan(arr), fill, arr)
    
    elif method == 'savgol':
        # expected kwargs: window_length, polyorder
        win = int(kwargs.get('window_length', 11))
        if win % 2 == 0:
            win += 1
            print(f"Warning: 'window_length' must be odd for Savitzky-Golay, increased to {win}.")
        poly = int(kwargs.get("polyorder", 3))
        if poly >= win:
            raise ValueError('Savgol: polyorder must be less than window_length')
        return savgol_filter(arr, window_length=win, polyorder=poly, mode='interp')
    elif method == "gaussian":
        # expected kwargs: std, avoids explicit window
        std = float(kwargs.get("std", 2.0))
        return gaussian(arr, sigma=std, mode="reflect")
    else:
        raise ValueError(f"Unknown method: {method} is not a valid method for 1D denoising.")

def assign_segments_hybrid(n, valley_idx=None, res_bp=10_000, fallback_span_bp=500_000):
    labels = np.zeros(n, dtype=int)
    if valley_idx is not None and len(valley_idx) > 0:
        cuts = np.sort(np.asarray(valley_idx, dtype=int))
        labels[:] = 1
        for i, c in enumerate(cuts, start=1):
            labels[c:] = i + 1
    else:
        bin_size_bins = max(1, int(round(fallback_span_bp / res_bp)))
        n_bins = int(np.ceil(n / bin_size_bins))
        for b in range(n_bins):
            start = b * bin_size_bins
            end = min((b + 1) * bin_size_bins, n)
            labels[start:end] = b + 1
    return labels

def peakfinder(df: pandas.DataFrame,*, method: str = 'peakdetect', denoise: str = None, denoise_kwargs: dict = None, **kwargs):
    """Extracts rows of data (diagonals in the original contact matrix) and finds peaks in each row using the methods in findpeaks package.
    
    Parameters
    ----------
    
    df: DataFrame
        Expects format generated from diag.get_aligned.

    method: str, optional
        Determines findpeaks method:
        - 'peakdetect': simple method for detecting peaks one scale at a time.
        - 'topology': method using persistent homology to find most significant peaks.
        - 'caerus': takes longer, looks across multiple scales.
    
    denoise: str, optional
        Preprocessing method for signals ("savgol", "gaussian", or None).
        
    denoise_kwargs: dict, optional
        Arguments for the preprocessing method.
        - savgol: window_length, polyorder
        - gaussian: window, std

    **kwargs: optional
        Arguments passed to findpeaks depending on method:
        - lookahead: int, determines scale of peaks found with peakdetect.
        - interpolate: str, applies interpolation/smoothing to valid methods.

    Returns
    -------
    out_frame: DataFrame
        Lists coordinates of detected peaks
    signals: dict
        stores HiC diagonals for reference, with preprocessing if included
    """
    # TODO: make sure OUTPUT STANDARDIZED

    # kwargs reserved keys
    fallback_span_bp = kwargs.get('fallback_span_bp', 500_000)
    res_bp = kwargs.get('res_bp', 10_000)

    fp = findpeaks(method=method, **(denoise_kwargs or {}))
    peaks = []
    signals = {}
    
    for d, row in df.iterrows():
        sig = np.array(row, dtype = float)
        # peaks_idx, props = debug_pd(sig, diagonal=d, lookahead=50)

        sig = preprocess_signal(sig, method=denoise, **(denoise_kwargs or {}))

        # for catching lack of valleys
        valley_idx, _ = find_peaks(-sig)
        labx = assign_segments_hybrid(
        n=len(sig),
        valley_idx=valley_idx if len(valley_idx) > 0 else None,
        res_bp=kwargs.get('res_bp', 10_000),
        fallback_span_bp=kwargs.get('fallback_span_bp', 500_000)
        )

        if method == 'peakdetect':
            look = kwargs.get('lookahead', None)
            if look is None:
                raise ValueError("Method 'peakdetect' requires 'lookahead' argument.")
            if len(sig) < (2 * look + 1):
                print(f'Diagonal {d} skipped: not enough points for lookahead {look}.')
                continue
        elif method in {'topology', 'caerus'}:
            if len(sig) < 5:
                print(f'Diagonal {d} skipped: {len(sig)} is not enough points for {method}.')

        result = fp.fit(sig)
        signals[d] = sig

        print(result['df'].head())

        peaks_df = result['df'].loc[result['df']['peak'] == 1, ['x', 'y']]
        peaks_df['labx'] = labx[peaks_df['x'].astype(int).values]

        if kwargs.get('interpolate'):
            factor = kwargs['interpolate']
            peaks_df['peak_x_interp'] = peaks_df['x'] # true interpolated coordinate
            peaks_df['peak_x'] = peaks_df['x'] / factor # map back to original scale, but keep float
        else:
            peaks_df['peak_x'] = peaks_df['x']

        for _, peak in peaks_df.iterrows():
            peaks.append({
                'diagonal': int(str(d).replace("d=", "")),
                'peak_x': int(peak['x']),
                'peak_y': float(peak['y']),
                'labx': int(peak['labx'])
            })
        
        print(f"Diagonal {d}: {len(peaks_df)} peaks found (signal length of {len(sig)})")
    
    out_frame = pandas.DataFrame(peaks)
    print(f'Output dataframe of size: {out_frame.size}')
    
    return out_frame, signals


def peak_linegraph(df: pandas.DataFrame, signals: dict, *, start: int, end: int):
    for d, sig in signals.items():
        sig_x = np.arange(start, min(end, len(sig)))
        sig_y = np.array(sig)[start:min(end, len(sig))]
        plt.plot(sig_x, sig_y, label=d) #may have to change label
        
        peaks_df = df[df['diagonal'] == d]
        peaks_in_range = peaks_df[(peaks_df['peak_x'] >= start) & (peaks_df['peak_x'] < end)]

        if not peaks_in_range.empty:
            plt.scatter(peaks_in_range['peak_x'], peaks_in_range['peak_y'], color='red', s=15, zorder=3)
    
    plt.xlabel('Position along chromosome')
    plt.ylabel('No. of contacts')
    plt.title(f'Diagonals with peaks {start}-{end}')
    plt.legend(loc = 'best', fontsize = 'small')
    plt.show()

def peak_dotplot(df: pandas.DataFrame):
    """Simple function to create a dotplot of peaks by position versus offset"""
    #df['diagonal'] = df['diagonal'].str.replace("d=", "").astype(int)

    plt.scatter(df['peak_x'], df['diagonal'])
    plt.xlabel('Position along chromosome')
    plt.ylabel('Offset')
    plt.title('Peak locations by offset')
    plt.show()

