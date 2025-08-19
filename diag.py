import cooler
import numpy as np
import pandas
import matplotlib.pyplot as plt

class DiagCool(cooler.Cooler):
    """Wrapper for 'cooler' library with custom utilities for working with Hi-C contact matrices.

    Internally holds a 'cooler.Cooler' object and provides methods for analysis of diagonals. Made for the 'diagonalgo' tool for identification of chromatin fountains and other structures.

    Attributes
    ----------
    cooler: cooler.Cooler - the underlying cooler object. Can be accessed directly if needed.

    Methods
    -------
    matrix(balance=False): access the raw contact matrix.
    fetch(chrom): fetch a chromosome-specific submatrix.
    get_diags(chrom, start, end, min_off, max_off): returns a list of the raw diagonals for a specified range of indices and offsets.
    get_aligned(chrom, min_off, max_off): returns the diagonals of an intrachromosomal contact matrix for a specified range of offsets, symmetrically trimmed to the shortest diagonal to mantain alignment.
    graph_diags(chrom, start, end, min_off, max_off): reads raw diagonals and graphs the signals.
    graph_aligned(chrom, start, end, min_off, max_off): incomplete. Meant to replace graph_diags for aligned and trimmed diagonals.
    """
    def __init__(self, uri):
        """Calls the cooler.Cooler init."""
        super().__init__(uri)
    
    def cooler(self):
        """Return underlying cooler object."""
        return self._cool
    
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
        print(f'Processing contact matrix of shape {mat.shape}')

        data = {}

        for d in range(min_off, max_off + 1, 2): # switch back to only even offsets for now, deal with odd later
            if d < size:
                antidiag = [mat[i, i + d] for i in range(size - d)]
            else:
                antidiag = []

            print(f'd={d}, length={len(antidiag)}')
            if len(antidiag) > 0 and not np.isnan(antidiag).any():
                data[f'd={d}'] = antidiag

        if not data:
            print('No valid data found.')
            return pandas.DataFrame()
        
        min_len = min(len(v) for v in data.values())

        # symmetrically trimmed to shortest length
        def trim_centered(lst, target_len):
            total_trim = len(lst) - target_len
            start = total_trim // 2
            end = start + target_len
            return lst[start:end]

        trimmed = {k: trim_centered(v, min_len) for k, v in data.items()}

        df = pandas.DataFrame.from_dict(trimmed, orient = 'index')
        df.index.name = 'offset'
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
        # intended to graph output of get_aligned before processing

# clr = DiagCool('/Users/hzhang/repli-HiC_data/Repli-HiC_K562_WT_totalS.mcool::resolutions/10000')

# clr.graph_diags('1', 1000, 1200, 2, 12) # the above example. not aligned but compare numerical results

# these are all comparable to diagrams in the paper
# clr.graph_diags('16', 775, 850, 3, 15, 3)
# clr.graph_diags('16', 788, 801, 2, 12)
# clr.graph_diags('20', 172, 183, 2, 10, 1)

