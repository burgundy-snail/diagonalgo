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

clr10 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000') # C:/Users/hzhan/OneDrive/Documents/Curie_internship/data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000
clr25 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/25000')
clr50 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/50000')

#TODO: test preprocessing methods (interpolation, denoising, etc) and use of topology + caerus methods, generate dotplots
chr6_highres = clr10.get_aligned('6', 16, 50)
chr6 = clr25.get_aligned('6', 8, 26)
chr6_lowres = clr50.get_aligned('6', 4, 13)
chr16 = clr10.get_aligned('16', 4, 12)
chry = clr50.get_aligned('Y', 6, )

# selecting same areas as figures shown in Liu et al., 2024 to compare detection
peaks_16c, sig_16c = peakfinder(chr16, 'caerus')
peak_linegraph(peaks_16c, sig_16c, start=7750, end=8500)
peaks_16c.to_csv('/output/chr16_caerus.tsv', sep='\t')
peak6c, sig6c = peakfinder(chr6, 'caerus')
peak_linegraph(peak6c, sig6c, start=1860, end=2020)
peak6c.to_csv('/output/chr6_caerus.tsv', sep='\t')
#peak6, sig6 = peakfinder(chr6_lowres, 10)
#peak_linegraph(peak6, sig6, start=930, end=1010)

