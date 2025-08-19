import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

clr = cooler.Cooler('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000')

# original method (manually change d and plot each time)
HiC_map = clr.matrix(balance=False).fetch('1')
d = 2
p = HiC_map.diagonal(d)
plt.plot(p[1000:1200])
d = 4
p = HiC_map.diagonal(d)
plt.plot(p[1000:1200])
d = 6
p = HiC_map.diagonal(d)
plt.plot(p[1000:1200])
d = 8
p = HiC_map.diagonal(d)
plt.plot(p[1000:1200])
plt.show()