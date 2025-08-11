import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from findpeaks import findpeaks

from diag import MyCool

def peakfinder(df: pandas.DataFrame, look: float):
    fp = findpeaks(method='peakdetect', lookahead=look, interpolate=None)
    peaks = {}
    
    for row in df.itertuples(index = True):
        d = row[0]
        sig = np.array(row[1:], dtype = float)
        result = fp.fit(sig)
        peaks_df = result['df'].loc[result['df']['peak'] == 1, ['x', 'y']]

        for _, peak in peaks_df.iterrows():
            peaks.append({
                'diagonal': d,
                'peak_x': int(peak['x']),
                'peak_y': float(peak['y'])
            })
    
    out_frame = pandas.DataFrame(peaks_df) # may need additional modification we'll see
    print(f"Output dataframe of size: {out_frame.size}")
    
    return out_frame
# TODO: write code for dotplot, test & debug more!

def peakgraph(df: pandas.DataFrame, look: float):
    fp = findpeaks(method='peakdetect', lookahead=look, interpolate=None)
    peaks = {}
    
    for i, row in df.iterrows():
        sig = np.array(row[1:], dtype = float)
        result = fp.fit(sig)
        peaks_df = result['df'].loc[result['df']['peak'] == 1, ['x', 'y']]

        plt.plot(range(len(sig)), sig, label=f'd={i}')
        plt.scatter(peaks_df['x'], peaks_df['y'], color='red', s=40, zorder=3)

        plt.xlabel('Position along chromosome')
    plt.ylabel('Value')
    plt.title('All diagonals with peaks')
    plt.legend(loc = "best")
    plt.show()

clr = MyCool("/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000") # C:/Users/hzhan/OneDrive/Documents/Curie_internship/data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000

peaks15 = peakfinder(clr.get_aligned('1', 2, 20), 15)
peakgraph(clr.get_aligned('1', 2, 20), 15)
peaks15.to_csv("./output/out_peaks15.tsv", sep = '\t')

peaks10 = peakfinder(clr.get_aligned('1', 2, 20), 10)
peakgraph(clr.get_aligned('1', 2, 20), 10)
peaks10.to_csv("./output/out_peaks10.tsv", sep = '\t')

