import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from findpeaks import findpeaks

from diag import MyCool

def peakfinder(df: pandas.DataFrame):
    fp = findpeaks(method='peakdetect', lookahead=150, interpolate=None)
    peaks = {}
    
    for row in df.itertuples(index = True):
        d = list(row)
        result = fp.fit(d) # this is a dict
        peaks[row.Index] = result
    
    print(result.keys())
    out_frame = pandas.DataFrame(peaks) # may need additional modification we'll see
    
    return out_frame
# TODO: fix peakfinder. write code to make line graph and dotplot

clr = MyCool("/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000") # C:/Users/hzhan/OneDrive/Documents/Curie_internship/data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000

aaa = peakfinder(clr.get_aligned('1', 4, 50))
print(aaa)

