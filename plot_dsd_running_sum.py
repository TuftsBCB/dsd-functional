from matplotlib import pyplot as plt
import expt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("infile")

# optional: specify npy files containing DSD and overlap matrices previously computed
# if files do not exist, compute matrices and save them to file
parser.add_argument("--dsdMat")
parser.add_argument("--overlapMat")

# if specified, write LaTeX code into file as text
parser.add_argument("--texOutfile")

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
K = expt.overlapMatrix(nodeList, GOfile, options.overlapMat)

# flatten D and K, sort by distance in increasing order
DFlat = np.ravel(D)
DSorted = np.argsort(DFlat)

KFlat = np.ravel(K)
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

# make another file to plot numpairs
plt.plot(distances, overlapRSum)
plt.xlabel("DSD")
plt.ylabel("Running sum of function overlap")

# get organism name from GOfile
outfile = GOfile[:-11].split('\\')
outfile = outfile[-1]

outfile = "plots/" + outfile + "_dsd_running_sum"

print "saving plot.."
plt.savefig(outfile)

# append code to file -- assume there is already existing code that we dont want to overwrite
if options.texOutfile:
    with open(options.texOutfile, "a") as f:
        f.write("%DSD vs Running Sum Plot:%\n")
        f.write("\\subfloat[DSD vs Running Sum]{\n"
                "\\includegraphics[width=0.5\\textwidth]{" + outfile + ".png}\n"
                "}\n \n")