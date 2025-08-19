import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from findpeaks import findpeaks

from diag import DiagCool

def preprocess_signal(sig: np.ndarray, method: str = None, **kwargs) -> np.ndarray:
    """Optional preprocessing of 1D signal."""
    if method is None:
        return sig
    elif method == "savgol":
        # kwargs: window_length, polyorder
        return findpeaks.savgol_filter(sig, **kwargs)
    elif method == "gaussian":
        # kwargs: window, std
        window = findpeaks.gaussian(kwargs.get("window", 11), kwargs.get("std", 2))
        return findpeaks.convolve(sig, window / window.sum(), mode="same")
    else:
        raise ValueError(f"Unknown method: {method} is not a valid method for 1D denoising.")

def peakfinder(df: pandas.DataFrame, method: str = 'peakdetect', **kwargs):
    """Extracts rows of data (diagonals in the original contact matrix) and finds peaks in each row using the methods in findpeaks package.
    
    Parameters
    ----------
    df: DataFrame, expects format generated from diag.get_aligned]
    method: str, optional.
        Determines findpeaks method:
        - 'peakdetect': simple method for detecting peaks one scale at a time.
        - 'topology': method using persistent homology to find most significant peaks.
        - 'caerus': takes longer, looks across multiple scales.

    **kwargs: arguments passed to findpeaks depending on method
    - lookahead: int, determines scale of peaks found with peakdetect.
    - interpolate: str, applies interpolation/smoothing to valid methods.

    Returns
    -------
    out_frame: DataFrame listing coordinates of detected peaks
    signals: dictionary storing HiC diagonals for reference
    """
    # TODO: allow preprocessing
    fp = findpeaks(method=method, **kwargs)
    peaks = []
    signals = {}
    
    for d, row in df.iterrows():
        sig = np.array(row, dtype = float)

        if method == 'peakdetect':
            look = kwargs.get('lookahead', None)
            if look is None:
                raise ValueError("Method 'peakdetect' requires 'lookahead' argument.")
            if len(sig) < (2 * look + 1):
                print(f'Diagonal {d} skipped: not enough points for lookahead {look}.')
                continue
        elif method == 'topology':
            if len(sig) < 5:
                print(f'Diagonal {d} skipped: {len(sig)} is not enough points for topology.')
        elif method == 'caerus':
            if len(sig) < 5:
                print(f'Diagonal {d} skipped: {len(sig)} is not enough points for caerus.')

        result = fp.fit(sig)
        signals[d] = sig

        peaks_df = result['df'].loc[result['df']['peak'] == 1, ['x', 'y']] # just stores the coordinates of peaks.

        for _, peak in peaks_df.iterrows():
            peaks.append({
                'diagonal': d,
                'peak_x': int(peak['x']),
                'peak_y': float(peak['y'])
            })
    
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
    plt.scatter(df[:0], df[:1])
    plt.xlabel('Position along chromosome')
    plt.ylabel('Offset')
    plt.title('Peak locations by offset')
    plt.show()

