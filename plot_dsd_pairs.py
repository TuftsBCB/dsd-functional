from matplotlib import pyplot as plt
import expt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("infile")

# optional: specify npy files containing DSD matrix previously computed
# if file does not exist, compute matrix and save to file
parser.add_argument("--dsdMat")

# not using these atm
# parser.add_argument("trialNum", type=int)
# tNum = options.trialNum

options = parser.parse_args()

# need to decide on this
LMSetSize = 50

# assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
GOfile = options.infile[:-4]
GOfile += "_NCBI_to_GO"

G = expt.setupGraph(options.infile)

# capture canonical node order
nodeList = G.nodes()

LMSet = expt.getLMSet(G, nodeList, LMSetSize)

D = expt.getDSD(G, nodeList, LMSet, options.dsdMat)

# flatten D, sort by distance in increasing order
DFlat = np.ravel(D)
DSorted = np.argsort(DFlat)

numPairs = len(DSorted)

# dsd distances in increasing order
distances = [DFlat[k] for k in DSorted]

# make another file to plot numpairs
plt.plot(distances, range(numPairs))
plt.xlabel("DSD")
plt.ylabel("Running sum of protein pairs")
plt.show()