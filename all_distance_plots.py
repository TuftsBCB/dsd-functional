# these lines needed for working on ELF -- feel free to comment out
import sys
sys.path.append("/usr/lib64/python2.7/site-packages")
sys.path.append("/usr/lib/python2.7/site-packages/decorator-3.4.0-py2.7.egg")
sys.path.append("/usr/lib/python2.7/site-packages")

import plotting
import argparse
from matplotlib import pyplot as plt
import os

plt.rc('axes.formatter', use_mathtext=True)

parser = argparse.ArgumentParser()

# path to folder containing ppi and GO files
# organisms is a string of organisms we want to generate plots for, separated by spaces
parser.add_argument("directory")
parser.add_argument("organisms")

# # file to write LaTeX code to
# parser.add_argument("outfile")

options = parser.parse_args()

organisms = options.organisms.split(' ')

# ppi file extension
ppiExt = ".ppi"

# folder we want to put plots in
plotsDir = "plots"

# assume that any npy matrices will be in this folder
# if folder doesnt exist, create
npyDir = "NumPy files"
try:
    os.makedirs(npyDir)
except OSError:
    if not os.path.isdir(npyDir):
        raise

distanceMetrics = {'DSD': None,  'SPD': None}

figNum = 1
for organism in organisms:
    # setup
    print("generating plots for: " + organism)
    infile = options.directory + "/" + organism + ppiExt
    npyPath = npyDir + "/" + organism

    distanceMetrics['DSD'] = npyPath + "_dsd.npy"
    distanceMetrics['SPD'] = npyPath + "_spd.npy"

    plt.figure(figNum)
    plotting.all_distance_pairs_density(infile, distanceMetrics=distanceMetrics, overlapMat=npyPath + "_overlap.npy")

    # save plot
    figFile = plotsDir + "/" + organism + "_alld_pairs_density"
    figFile += ".png"
    plt.savefig(figFile)
    figNum += 1
