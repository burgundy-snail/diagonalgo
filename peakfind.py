import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from findpeaks import findpeaks

from diag import MyCool

def peakfinder(df: pandas.DataFrame, look: float):
    fp = findpeaks(method = 'peakdetect', lookahead = look, interpolate = None)
    peaks = []
    signals = {}
    
    for d, row in df.iterrows():
        sig = np.array(row, dtype = float)

        if len(sig) < (2 * look + 1):
            print(f'Not enough points for lookahead {look}, row skipped.')
            continue

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
# TODO: graphing code that limits range, makes use of signals

def peakgraph(df: pandas.DataFrame, look: float):
    fp = findpeaks(method = 'peakdetect', lookahead = look, interpolate = None)
    
    for i, row in df.iterrows():
        sig = np.array(row[1:], dtype = float)
        result = fp.fit(sig)
        peaks_df = result['df'].loc[result['df']['peak'] == 1, ['x', 'y']]

        plt.plot(range(len(sig)), sig, label=f'd={i}')
        plt.scatter(peaks_df['x'], peaks_df['y'], color='red', s=40, zorder=3)

    plt.xlabel('Position along chromosome')
    plt.ylabel('No. of contacts')
    plt.title('All diagonals with peaks')
    plt.legend(loc = 'best', fontsize = 'small')
    plt.show()

def dotplot(df: pandas.DataFrame):
    plt.scatter(df[:0], df[:1])
    plt.xlabel('Position along chromosome')
    plt.ylabel('Offset')
    plt.title('Peak locations by offset')
    plt.show()

clr = MyCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000') # C:/Users/hzhan/OneDrive/Documents/Curie_internship/data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000

diags = clr.get_aligned('1', 2, 20)
print(diags.shape)
print(diags.head())
# print([len(diags[col]) for col in diags.columns])

peaks150, sig150 = peakfinder(diags, 150)
peakgraph(diags, 150)
peaks150.to_csv('./output/out_peaks150.tsv', sep = '\t')

peaks100, sig100 = peakfinder(diags, 100)
peakgraph(diags, 100)
peaks100.to_csv('./output/out_peaks100.tsv', sep = '\t')

