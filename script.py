import numpy as np
import argparse
import random
import matplotlib as mp
import sys
import dsd
import networkx as nx
import time

def announce(str):
    print time.strftime('%H:%M'),str
    sys.stdout.flush()

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("trialNum", type=int)

# optional: specify npy files containing DSD and overlap matrices previously computed
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


if not options.dsdMat:
    #
    # construct hemat
    #
    HEmatrix = dsd.hematrix(np.array(nx.adjacency_matrix(G,nodeList).todense()))

    # construct DSD
    D = dsd.DSD(HEmatrix,LMset)

    npyOutfile = raw_input("DSD matrix computed. Enter file path to save matrix for future use, "
                           "or just press Enter to continue without saving:")
    if not npyOutfile:
        print "fine, be that way!"
    else:
        np.save(npyOutfile, D)
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

# assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
GOfile = options.infile[:-4]
GOfile += "_NCBI_to_GO"
node_labels = GOToDict(nodeList, GOfile)

if not options.overlapMat:
    # construct K = overlap matrix
    K = np.zeros((n, n), dtype=int)

    # if there is function overlap, set K[ij] to 1
    for i in range(1,n):
        nd = nodeList[i]
        for j in range(i+1,n):
            for label in node_labels[nd]:
                if label in node_labels[nodeList[j]]:
                    K[i][j] = 1
                    K[j][i] = 1
                    break
    npyOutfile = raw_input("overlap matrix computed. Enter file path to save matrix for future use, "
                           "or just press Enter to continue without saving:")
    if not npyOutfile:
        print "fine, be that way!"
    else:
        np.save(npyOutfile, K)
else:
    K = np.load(options.overlapMat)
