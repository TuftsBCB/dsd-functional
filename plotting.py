# just combining the various plotting routines into a single callable module

import expt
import numpy as np
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt

def dsd_overlap_pairs(infile, dsdMat=None, overlapMat=None, randomize=False):
    """
    plots dsd against (normalized):
        1. summed overlap
        2. summed pairs
    """
    # need to decide on this
    LMSetSize = 50

    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    LMSet = expt.get_LMset(G, nodeList, LMSetSize)

    D = expt.dsd_matrix(G, nodeList, LMSet, dsdMat)
    # if randomize:
    #     K = np.random.random_integers(0, 1, np.shape(D))
    #     K = (K + np.transpose(K)) / 2
    # else:
    K = expt.overlap_matrix(nodeList, GOfile, overlapMat, randomize=randomize)

    # flatten D and K, sort by distance in increasing order
    #DFlat = np.ravel(D)
    #KFlat = np.ravel(K)

    DFlat = expt.triu_ravel(D)
    KFlat = expt.triu_ravel(K)
    DSorted = np.argsort(DFlat)
    numPairs = len(DSorted)

    # running sum of overlap
    overlapRSum = np.zeros(numPairs)
    overlapSum = 0.0
    for n in range(numPairs):
        index = DSorted[n]
        overlapSum += KFlat[index]
        overlapRSum[n] = overlapSum

    # dsd distances in increasing order
    distances = [DFlat[k] for k in DSorted]

    # plot dsd vs rsum, then numpairs
    plt.plot(distances, overlapRSum/overlapRSum[-1], "b-", label="cumulative overlap")

    plt.plot(distances, np.array(range(numPairs))/float(numPairs), "r-", label="protein pairs")
    plt.xlabel("DSD")
    plt.legend()


def pairs_summed_overlap(infile, dsdMat=None, overlapMat=None, randomize=False):
    """
    plots pairs against running sum of overlap
    """
    # need to decide on this
    LMSetSize = 50

    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    LMSet = expt.get_LMset(G, nodeList, LMSetSize)

    D = expt.dsd_matrix(G, nodeList, LMSet, dsdMat)
    K = expt.overlap_matrix(nodeList, GOfile, overlapMat, randomize=randomize)

    # flatten D and K, sort by distance in increasing order
    #DFlat = np.ravel(D)
    #KFlat = np.ravel(K)

    DFlat = expt.triu_ravel(D)
    KFlat = expt.triu_ravel(K)
    DSorted = np.argsort(DFlat)
    numPairs = len(DSorted)

    # running sum of overlap
    overlapRSum = np.zeros(numPairs)
    overlapSum = 0.0
    for n in range(numPairs):
        index = DSorted[n]
        overlapSum += KFlat[index]
        overlapRSum[n] = overlapSum

    plt.plot(range(numPairs), overlapRSum)
    plt.plot([0, numPairs], [overlapRSum[0], overlapRSum[-1]], "k--")
    plt.xlabel("Protein pairs")
    plt.ylabel("Running sum of function overlap")


def dsd_density(infile, dsdMat=None, overlapMat=None, randomize=False):
    """
    plots dsd against overlap density
    """
    # need to decide on this
    LMSetSize = 50

    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    LMSet = expt.get_LMset(G, nodeList, LMSetSize)

    D = expt.dsd_matrix(G, nodeList, LMSet, dsdMat)
    K = expt.overlap_matrix(nodeList, GOfile, overlapMat, randomize=randomize)

    # flatten D and K, sort by distance in increasing order
    #DFlat = np.ravel(D)
    #KFlat = np.ravel(K)

    DFlat = expt.triu_ravel(D)
    KFlat = expt.triu_ravel(K)
    DSorted = np.argsort(DFlat)
    numPairs = len(DSorted)

    # list of (summed overlap):(number of pairs) ratios
    overlapRatios = np.zeros(numPairs)
    overlapSum = 0.0
    posCI = np.zeros(numPairs)
    negCI = np.zeros(numPairs)
    for n in range(numPairs):
        index = DSorted[n]
        overlapSum += KFlat[index]
        ratio = overlapSum/(n+1)
        overlapRatios[n] = ratio
        # 95% CI
        std = np.sqrt(ratio * (1-ratio) / (n+1))
        posCI[n] = ratio + (1.96*std)
        negCI[n] = ratio - (1.96*std)

    # dsd distances in increasing order
    distances = [DFlat[k] for k in DSorted]

    plt.plot(distances, overlapRatios, 'k-', distances, negCI, 'b--', distances, posCI, 'b--')
    # set lim to min and max of CI, ignoring the first 1% of pairs
    plt.ylim((np.min(negCI[(numPairs/100):])-0.01, np.max(posCI[(numPairs/100):])+0.01))
    plt.xlabel("DSD")
    plt.ylabel("Density of function overlap")

