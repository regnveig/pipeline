from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cooler
import scipy
import seaborn as sns
import statsmodels.stats.multitest
import pybedtools
import sys
from cooltools.lib.numutils import set_diag
import pickle


def CTALE_norm_multiplicate(mtx, ROI_start, ROI_end, resolution, func=scipy.stats.gmean, mult=1.54):
    """mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    func - function for calculating mean, can be numpy.mean,numpy.median and etc. Default: scipy.stats.gmean
    mult-coefficient of multiplication around ROI, default=1.54
    returns normalized matrix"""
    normalized = np.zeros(shape=mtx.shape)
    start_bin = int(ROI_start / resolution)
    end_bin = int(ROI_end / resolution)
    # fill main diagonal with zeros
    np.fill_diagonal(mtx, 0)
    # first we multiply around region
    new_mtx = mtx * mult
    new_mtx[start_bin:end_bin, start_bin:end_bin] = mtx[start_bin:end_bin, start_bin:end_bin]
    for i in range(start_bin, end_bin):
        # left=new_mtx[i,i:]
        # up=new_mtx[:i,i]
        normalized[i, i:] = new_mtx[i, i:] / (np.sum(new_mtx[i, i:]) + np.sum(new_mtx[:i, i]))  # left
        normalized[:i, i] = new_mtx[:i, i] / (np.sum(new_mtx[:i, i]) + np.sum(new_mtx[i, i:]))  # up
    for i in range(start_bin, end_bin):
        for j in range(i + 1, end_bin):
            normalized[i, j] = new_mtx[i, j] / func([np.sum(new_mtx[i, :]), np.sum(new_mtx[:, j])])

    i_lower = np.tril_indices(normalized.shape[0], -1)  # creates symmetric matrix
    normalized[i_lower] = normalized.T[i_lower]
    return (normalized)


def CTALE_norm(mtx, ROI_start, ROI_end, resolution, func=scipy.stats.gmean):
    """Same as CTALE_norm_multiplicate function but without multiplication.
    mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    func - function for calculating mean, can be numpy.mean,numpy.median and etc. Default scipy.stats.gmean
    returns normalized matrix"""
    normalized = np.zeros(shape=mtx.shape)
    start_bin = int(ROI_start / resolution)
    end_bin = int(ROI_end / resolution)
    # fill main diagonal with zeros
    np.fill_diagonal(mtx, 0)
    for i in range(start_bin, end_bin):
        # left=mtx[i,i:]
        # up=mtx[:i,i]
        normalized[i, i:] = mtx[i, i:] / (np.sum(mtx[i, i:]) + np.sum(mtx[:i, i]))  # left
        assert np.isclose(np.sum(mtx[i, i:]) + np.sum(mtx[:i, i]), np.sum(mtx[i]),
                                                       atol=np.average(mtx)/100000)
        normalized[:i, i] = mtx[:i, i] / (np.sum(mtx[:i, i]) + np.sum(mtx[i, i:]))  # up
    for i in range(start_bin, end_bin):
        for j in range(i + 1, end_bin):
            normalized[i, j] = mtx[i, j] / func([np.sum(mtx[i, :]), np.sum(mtx[:, j])])
    i_lower = np.tril_indices(normalized.shape[0], -1)  # creates symmetric matrix
    normalized[i_lower] = normalized.T[i_lower]
    return (normalized)


def multiplicate(mtx, ROI_start, ROI_end, resolution, mult=1.54):
    """mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    mult-coefficient of multiplication around ROI, default=1.54
    Function perform multiplicating of matrix around ROI for selected coefficient"""
    start_bin = int(ROI_start / resolution)
    end_bin = int(ROI_end / resolution)
    # fill main diagonal with zeros
    np.fill_diagonal(mtx, 0)
    # first we multiply around region
    new_mtx = mtx * mult
    new_mtx[start_bin:end_bin, start_bin:end_bin] = mtx[start_bin:end_bin, start_bin:end_bin]
    return (new_mtx)

def CTALE_norm_iterative(mtx,ROI_start,ROI_end,resolution,func=scipy.stats.gmean,mult=1.54,steps=20,tolerance=1e-5):
    """Main function that perform normalization until variance>tolerance
    mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    func - function for calculating mean, can be numpy.mean,numpy.median and etc. Default: scipy.stats.gmean
    mult-coefficient of multiplication around ROI, default=1.54
    steps-number of iterations, by default=20
    tolerance-when variance<tolerance algorithm stops.
    returns normalized matrix"""
    start_bin=int(ROI_start/resolution)
    end_bin=int(ROI_end/resolution)
    out=multiplicate(mtx=mtx,ROI_start=ROI_start,ROI_end=ROI_end,resolution=resolution,mult=mult)
    for s in range(steps):
        out=CTALE_norm(out,ROI_start,ROI_end,resolution,func=func)
        var=np.var(np.sum(out[start_bin:end_bin,:],axis=1))
        print('Variance is: ',var)
        if var<tolerance:
            print('Variance < ',tolerance)
            break
    return(out)
def Save_coolfile(coolfile,mtx,output_coolfile,genome):
    """Function change raw HiC matrix of cool file to user selected (normalized) and write it to new file.
    Because function rewrite data of original cool, later you should load it with balance=False flag.
    coolfile - original HiC file
    mtx - matrix to write
    output_coolfile - name of new cool file
    genome - genome assembly id"""
    #create bins
    bins=coolfile.bins()[0:]
    #Create sparse matrix
    mtx_upper_diag = np.triu(mtx, k=0)
    smtx=scipy.sparse.csr_matrix(mtx_upper_diag)
    sc=smtx.tocoo(copy=False)
    smtx_pixels=pd.DataFrame({'bin1_id': sc.row, 'bin2_id': sc.col, 'count': sc.data})
    cooler.io.create(output_coolfile,bins,smtx_pixels,assembly=genome,dtype={'count':float})
    return('Saved')

def draw_graph(mat, diagonal_offset, start_bin, end_bin, name):
    clipval = np.nanmedian(np.diag(mat, diagonal_offset))
    print(clipval)
    for i in range(-diagonal_offset + 1, diagonal_offset):
        set_diag(mat, clipval, i)
    plt.imshow(np.log(mat[start_bin:end_bin, start_bin:end_bin]), cmap="autumn_r")  # cmap="RdBu_r")
    plt.colorbar()
    plt.savefig(name)
    plt.clf()

# n.py raw.cool[1] start_cap[2] end_cap[3] bin_size[4] chr_cap[5] norm.cool[6]

chr = sys.argv[5]
start = int(sys.argv[2])
end = int(sys.argv[3])
bin_size = int(sys.argv[4])
start_bin = int(start/bin_size)
end_bin = int(end/bin_size)

#load data
raw = cooler.Cooler(sys.argv[1])
mtx_raw = raw.matrix(balance=False).fetch(chr)

#Perform normalization
mtx_normalized=CTALE_norm_iterative(mtx_raw, start, end, bin_size, steps=20, mult=1)
mtx = open('rawInv.pkl', 'wb')
pickle.dump(mtx_raw,mtx)
mtx.close()
mtx = open('normalizedInv.pkl','wb')
pickle.dump(mtx_normalized,mtx)

#Save_coolfile
Save_coolfile(raw,mtx_normalized,sys.argv[6],raw.info[u'genome-assembly'])

#plot maps
draw_graph(mtx_raw, 2, 10804, 10892, 'raw_Inv')
draw_graph(mtx_normalized, 2, 10804, 10892, 'normalized_Inv')
