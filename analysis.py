import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from findpeaks import findpeaks

from diag import DiagCool
from peakfind import peakfinder, peak_linegraph

# clr = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000')

# clr.graph_diags('1', 1000, 1200, 2, 12) # the above example. not aligned but compare numerical results

# these are all comparable to diagrams in the paper
# clr.graph_diags('16', 775, 850, 3, 15, 3)
# clr.graph_diags('16', 788, 801, 2, 12)
# clr.graph_diags('20', 172, 183, 2, 10, 1)

clr10 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000') # C:/Users/hzhan/OneDrive/Documents/Curie_internship/data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000
clr25 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/25000')
clr50 = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/50000')

#TODO: test preprocessing methods (interpolation, denoising, etc) and use of topology + caerus methods, generate dotplots
chr6_highres = clr10.get_aligned('6', 16, 50)
chr6 = clr25.get_aligned('6', 8, 26)
#chr6_lowres = clr50.get_aligned('6', 4, 13)
chr16 = clr10.get_aligned('16', 4, 12)
#chry = clr50.get_aligned('Y', 1000, 1400) #small dimensions for testing errors

# selecting same areas as figures shown in Liu et al., 2024 to compare detection
# these do take a while to run, TODO figure out if process can be streamlined
peaks_16c, sig_16c = peakfinder(chr16, interpolate=50, lookahead=100)
peak_linegraph(peaks_16c, sig_16c, start=7750, end=8500)
peaks_16c.to_csv('./output/chr16_pd_intpol50.tsv', sep='\t')
peak6c, sig6c = peakfinder(chr6, interpolate=50, lookahead=100)
peak_linegraph(peak6c, sig6c, start=1860, end=2020)
peak6c.to_csv('./output/chr6_pd_intpol50.tsv', sep='\t')
#peak6, sig6 = peakfinder(chr6_lowres, 10)
#peak_linegraph(peak6, sig6, start=930, end=1010)
# peaks_y, sig_y = peakfinder(chry, 'caerus')
# peak_linegraph(peaks_y, sig_y, start=930, end=1010)
# peaks_y.to_csv('./output/chry_caerus.tsv', sep='\t')
