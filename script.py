import numpy as np
import argparse
import random
import matplotlib as mp
from matplotlib import pyplot as plt
import sys
import dsd
import networkx as nx
import time
import os.path

def announce(str):
    print time.strftime('%H:%M'),str
    sys.stdout.flush()

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("trialNum", type=int)

# optional: specify npy files containing DSD and overlap matrices previously computed
# if files do not exist, compute matrices and save them to file
parser.add_argument("--dsdMat")
parser.add_argument("--overlapMat")

options = parser.parse_args()

tNum = options.trialNum

# need to decide on this
LMSetSize = 50

# test all nodes
# nSamps = 0

# get input graph as a networkx graph
G = nx.read_adjlist(options.infile)

# remove selfloop edges
G.remove_edges_from(G.selfloop_edges())

# we can only work with a connected graph
# DONT Do this: dsd.make_connected(G)
# instead, select the largest connected component
G = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]

# capture canonical node order
nodeList = G.nodes()

# initial bookkeeping
n = G.number_of_nodes()
degree = np.array([G.degree(i) for i in nodeList],dtype='float')

# sort nodes by degree
degAscending = np.argsort(degree)
degDescending = degAscending[::-1]

# use the highest degree nodes as landmarks in each case
# this is based on the assumption that high-degree nodes will be 
# easiest to match using sequence or other side information
LMset = degDescending[0:LMSetSize]

# make sure sample set and landmark set don't have nodes in common
# put this back in if/when we are doing cross-species because matching landmarks is misleading
# nonLMset = [i for i in xrange(n) if i not in LMset]
# if (nSamps != 0):
#    sampleSet = random.sample(nonLMset,nSamps)
# else:

# modified for testing ALL nodes
# not using this right now
# sampleSet = nodeList
# nSamps = len(sampleSet)

# if npy path not entered, or file does not exist, compute D
if not options.dsdMat or not os.path.isfile(options.dsdMat):
    #
    # construct hemat
    #
    HEmatrix = dsd.hematrix(np.array(nx.adjacency_matrix(G,nodeList).todense()))

    # construct DSD
    D = dsd.DSD(HEmatrix,LMset)

    if options.dsdMat:
        np.save(options.dsdMat, D)
else:
    D = np.load(options.dsdMat)


def GOToDict(nodeList, GOfile):
    """
    given set of nodes and name of GO file,
    return dictionary with node(gene id) as key, list of labels as value
    """
    node_labels = {}

    with open(GOfile) as f:
        for line in f:
            # assume first key in each line is either gene id or empty (a tab char)
            # only consider lines with gene ids
            if line[0] == '\t':
                continue
            part = line.strip().split('\t', 1)

            # if no label, add empty list
            if len(part) == 2:
                labels = part[1].split('; ')
            else:
                labels = []
            node_labels[part[0]] = labels

            # consider any node not listed in GO file as having no labels
            for node in nodeList:
                if node not in node_labels:
                    node_labels[node] = []

    return node_labels

def setOverlap(set1, set2):
    return not set1.isdisjoint(set2)

# assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
GOfile = options.infile[:-4]
GOfile += "_NCBI_to_GO"
node_labels = GOToDict(nodeList, GOfile)

# if npy path not entered, or file does not exist, compute K
if not options.overlapMat or not os.path.isfile(options.overlapMat):
    # construct K = overlap matrix
    K = np.zeros((n, n), dtype=int)

    # if there is function overlap, set K[ij] to 1
    for i in xrange(n):
        for j in xrange(i+1,n):
            K[i][j] = K[j][i] = setOverlap(set(nodeList[i]), set(nodeList[j]))

    if options.overlapMat:
        np.save(options.overlapMat, K)
else:
    K = np.load(options.overlapMat)


#
# Code for initial plot
#

# flatten D and K, sort by distance in increasing order
DFlat = np.ravel(D)
DSorted = np.argsort(DFlat)

KFlat = np.ravel(K)
numPairs = len(DSorted)
# list of (summed overlap):(number of pairs) ratios
overlapRatios = np.zeros(numPairs)
overlapSum = 0.0
for n in range(numPairs):
    index = DSorted[n]
    overlapSum += KFlat[index]
    overlapRatios[n] = overlapSum/(n+1)

# dsd distances in increasing order
distances = [DFlat[k] for k in DSorted]

plt.plot(distances, overlapRatios)
plt.xlabel("DSD")
plt.ylabel("Cumulative function overlap")
plt.show()
