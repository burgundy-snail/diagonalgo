import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

class MyCool(cooler.Cooler):
    def __init__(self, uri): # created the same way as a normal cooler
        super().__init__(uri)
    
    # need some way to align diags for both of these methods

    def get_diags(self, chrom: str, start: int, end: int, min_off: int, max_off: int):
        mat = self.matrix(balance = False).fetch(chrom)
        data = []
        data1 = []
        # ind = 0
        
        for d in range(min_off, max_off + 1, 2):
            p = mat.diagonal(d)

            if end > len(p):
                continue

            data.append(p[start:end])
            #rewrite this to align and trim
            

            # ind += 1

        for d in range(min_off+1, max_off + 2, 2):
            p = mat.diagonal(d)

            if end > len(p):
                continue

            data1.append(p[start:end])
        
        # trim to align, put data and data1 together somehow
        return data
    
    def get_aligned(self, chrom: str, min_off: int, max_off: int):
        mat = self.matrix(balance = False).fetch(chrom)
        size = mat.shape[0]

        data = {}

        for d in range(min_off, max_off + 1):
            antidiag = [
                mat[i, d-i]
                for i in range(max(0, d - size + 1), min(size, d + 1))
                if 0 <= d - i < size
            ]

            print(len(antidiag))
            if len(antidiag) > 0 and not np.isnan(antidiag).any():
                data[f'd={d}'] = antidiag

        if not data:
            print("No valid data found.")
            return pandas.DataFrame()
        
        '''min_len = min(len(v) for v in data.values())
        trimmed = {k: v[:min_len] for k, v in data.items()}'''

        # SYMMETRICAL TRIMMING THAT IS ACTUALLY GOOD

        df = pandas.DataFrame(trimmed)
        return df

    # looping version of the original code
    def graph_diags(self, chrom: str, start: int, end: int, min_off: int, max_off: int):
        mat = self.matrix(balance = False).fetch(chrom)

        for d in range(min_off, max_off + 1):

            p = mat.diagonal(d)

            plt.plot(p[start:end], label=f'd={d}')

        plt.legend()
        plt.show() # could also add labels etc
    
    def graph_aligned(self, chrom: str, start: int, end: int, min_off: int, max_off: int):
        mat = self.matrix(balance = False).fetch(chrom)
        # just do a get_aligned, then graph a segment of the data?
        # this depends on how easy the format of the data is to graph...
        # fix the goddamn get_aligned first bro
        # then think about efficiency, if get_aligned takes too long j first slice for visualization

clr = MyCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000')

# original method (manually change d and plot each time)
'''HiC_map = clr.matrix(balance=False).fetch('1')
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
plt.show()'''

clr.graph_diags('1', 1000, 1200, 2, 12) # the above example. not aligned but compare numerical results

# these are all from the diagram
# clr.graph_diags('16', 775, 850, 3, 15, 3)
# clr.graph_diags('16', 788, 801, 2, 12)
# clr.graph_diags('20', 172, 183, 2, 10, 1)

#mat1 = pandas.DataFrame(clr.get_diags('16', 788, 801, 2, 12, 1))
#mat1.to_csv("./output/out.tsv", sep = '\t')

mat2 = clr.get_aligned('1', 4, 50)
print(mat2.shape)
print(mat2.head())
mat2.to_csv("./output/out_aligned.tsv", sep = '\t')

