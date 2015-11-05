# just combining the various plotting routines into a single callable module
# plotting functions have optional matrix params -- these are paths to .npy binary files (matrices)
# for all functions:
#    if sp is true, use shortest-path distance as metric instead of DSD.
#    -- in this case, assume dsdMat, if specified, points to path of sp matrix.

import expt
import numpy as np
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt

def dsd_overlap_pairs(infile, dsdMat=None, overlapMat=None, randomize=False, sp=False):
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

    if sp:
        D = expt.sp_matrix(G, nodeList, dsdMat)
        metric = "Shortest-path distance"
    else:
        D = expt.dsd_matrix(G, nodeList, dsdMat)
        metric = "DSD"

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
    plt.xlabel(metric)
    plt.legend()


def pairs_summed_overlap(infile, dsdMat=None, overlapMat=None, randomize=False):
    """
    plots pairs against running sum of overlap
    """
    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    D = expt.dsd_matrix(G, nodeList, dsdMat)
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


def dsd_density(infile, dsdMat=None, overlapMat=None, randomize=False, sp=False):
    """
    plots dsd against overlap density
    """
    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    if sp:
        D = expt.sp_matrix(G, nodeList, dsdMat)
        metric = "Shortest-path distance"
    else:
        D = expt.dsd_matrix(G, nodeList, dsdMat)
        metric = "DSD"

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
    # set lim to min and max of CI, ignoring the first 1000 pairs
    plt.ylim((np.min(negCI[1000:])-0.01, np.max(posCI[1000:])+0.01))
    plt.xlabel(metric)
    plt.ylabel("Density of function overlap")


def dsd_density_res(infile, calc, dsdMat=None, resMat=None, randomize=False, sp=False):
    """
    plots dsd against resnik score density
    takes in .ppi file (infile) and SemSimCalculator instance (calc)
    """
    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    if sp:
        D = expt.sp_matrix(G, nodeList, dsdMat)
        metric = "Shortest-path distance"
    else:
        D = expt.dsd_matrix(G, nodeList, dsdMat)
        metric = "DSD"
    R = expt.resnik_matrix(nodeList, GOfile, calc, resMat, randomize=randomize)

    # flatten D and K, sort by distance in increasing order
    #DFlat = np.ravel(D)
    #KFlat = np.ravel(K)

    DFlat = expt.triu_ravel(D)
    RFlat = expt.triu_ravel(R)
    DSorted = np.argsort(DFlat)
    numPairs = len(DSorted)

    # list of (summed overlap):(number of pairs) ratios
    overlapRatios = np.zeros(numPairs)
    overlapSum = 0.0
    posCI = np.zeros(numPairs)
    negCI = np.zeros(numPairs)
    for n in range(numPairs):
        index = DSorted[n]
        overlapSum += RFlat[index]
        ratio = overlapSum/(n+1)
        overlapRatios[n] = ratio
        # 95% CI
        std = np.sqrt(ratio * (1-ratio) / (n+1))
        posCI[n] = ratio + (1.96*std)
        negCI[n] = ratio - (1.96*std)

    # dsd distances in increasing order
    distances = [DFlat[k] for k in DSorted]

    plt.plot(distances, overlapRatios, 'k-', distances, negCI, 'b--', distances, posCI, 'b--')
    # set lim to min and max of CI, ignoring the first 1000 pairs
    plt.ylim((np.min(negCI[1000:])-0.01, np.max(posCI[1000:])+0.01))
    plt.xlabel(metric)
    plt.ylabel("Resnik similarity score density")


ALL_DISTANCE_METRICS = {'DSD 50LM': None, 'DSD 50LM (random)' : None, 'DSD 500LM':None,\
                        'DSD 500LM (random)': None, 'DSD 200LM':None, 'DSD 200LM (random)': None,\
                        'SPD': None, 'DFD': None, 'DSD 500LM (eigen)': None, 'DSD 500LM (katz)': None,\
                        'DSD 500LM (closeness)': None, 'DSD 500LM (betweenness)': None}

def compute_all_distance_pairs_density(infile, distanceMetrics=ALL_DISTANCE_METRICS, overlapMat=None, randomize=False):
    """
    distanceMetrics is a dictionary containing distance metrics as keys and
    paths to their respective NumPy matrices, or None, as values.
    e.g. -- {'dsd': 'NumPy files/rat_dsd.npy'}
    Default is all supported metrics, with no saved matrices.

    Sample function call:
    compute_all_distance_pairs_density('PPIs and GO/rat.ppi', distanceMetrics={'DSD': 'NumPy files/rat_dsd.npy',\
     'SPD': 'NumPy files/rat_spd.npy', 'DFD': 'NumPy files/rat_dfd.npy'}, \
    overlapMat='NumPy files/rat_overlap.npy')

    return a dictionary of (DFlat, DSorted, overlapRatios) for each metric
    """
    metricDict = {}

    # mapping from metric to function that returns its matrix.
    # Assume each function can be called with identical parameters.
    metricMatrices = {'DSD 50LM' : expt.dsd_matrix, \
                'DSD 50LM (random)' : lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=50, centralityFunc='random'), \
                'SPD' : expt.sp_matrix, \
                'DFD': expt.diffusion_matrix, \
                'DSD 500LM': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize = 500),\
                'DSD 500LM (random)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=500, centralityFunc='random'), \
                'DSD 200LM': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize = 200),\
                'DSD 200LM (random)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=200, centralityFunc='random'), \
                'DSD 500LM (eigen)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=500, centralityFunc='eigenvector'), \
                'DSD 500LM (katz)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=500, centralityFunc='katz'), \
                'DSD 500LM (closeness)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=500, centralityFunc='closeness'), \
                'DSD 500LM (betweenness)': lambda G, nodeList, npyFile: expt.dsd_matrix(G, nodeList, npyFile, LMsetSize=500, centralityFunc='betweenness')}


    # assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
    GOfile = infile[:-4]
    GOfile += "_NCBI_to_GO"

    G = expt.setup_graph(infile)

    # capture canonical node order
    nodeList = G.nodes()

    K = expt.overlap_matrix(nodeList, GOfile, overlapMat, randomize=randomize)
    KFlat = expt.triu_ravel(K)
    numPairs = len(KFlat)

    for metric in distanceMetrics:
        if metric not in metricMatrices:
            print(metric + 'not in metricMatrices')
            continue
        D = metricMatrices[metric](G, nodeList, distanceMetrics[metric])
        DFlat = expt.triu_ravel(D)
        DSorted = np.argsort(DFlat)
        # list of (summed overlap):(number of pairs) ratios
        overlapRatios = np.zeros(numPairs)
        overlapSum = 0.0
        for n in range(numPairs):
            index = DSorted[n]
            overlapSum += KFlat[index]
            ratio = overlapSum/(n+1)
            overlapRatios[n] = ratio

        metricDict[metric] = (DFlat, DSorted, overlapRatios)

    return metricDict


def all_distance_pairs_density(infile, distanceMetrics=ALL_DISTANCE_METRICS, overlapMat=None, randomize=False):
    """
    Plots sorted pairs against density for various distance metrics.
    Currently, only DSD ('DSD') and SPD ('SPD') are supported.
    Experimenting with DFD (diffusion distance)

    distanceMetrics is a dictionary containing distance metrics as keys and
    paths to their respective NumPy matrices, or None, as values.
    e.g. -- {'dsd': 'NumPy files/rat_dsd.npy'}
    Default is all supported metrics, with no saved matrices.

    Sample function call:
    all_distance_pairs_density('PPIs and GO/rat.ppi', distanceMetrics={'DSD': 'NumPy files/rat_dsd.npy',\
     'SPD': 'NumPy files/rat_spd.npy', 'DFD': 'NumPy files/rat_dfd.npy'}, \
    overlapMat='NumPy files/rat_overlap.npy')
    """
    metricDict = compute_all_distance_pairs_density(infile, distanceMetrics, overlapMat, randomize)

    lower = 1.0
    upper = 0.0
    for metric in distanceMetrics:
        try:
            DFlat, DSorted, overlapRatios = metricDict[metric]
        except KeyError:
            continue

        plt.plot(range(len(overlapRatios)), overlapRatios, label=metric)
        lower = min(lower, np.min(overlapRatios[1000:]))
        upper = max(upper, np.max(overlapRatios[1000:]))

    plt.ylim(lower-0.01, upper+0.01)

    plt.xlabel("Sorted pairs in order of increasing distance")
    plt.ylabel("Density of function overlap")
    plt.legend()
