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
parser.add_argument("trialNum",type=int)
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
#
# construct hemat
#
HEmatrix = dsd.hematrix(np.array(nx.adjacency_matrix(G,nodeList).todense()))

# construct DSD
D = dsd.DSD(HEmatrix,LMset)

# get list of nodes in LMSet
LMnodes = [nodeList[index] for index in LMset]

def GOToDict(LMnodes, GOfile):
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

            # if gene id in our landmark set, add labels to dict
            # if no label, add empty list
            if part[0] not in LMnodes:
                continue
            if len(part) == 2:
                labels = part[1].split('; ')
            else:
                labels = []
            node_labels[part[0]] = labels

    # if node not found in GO file, consider it as having no labels
    for node in LMnodes:
        if node not in node_labels:
            node_labels[node] = []

    return node_labels

# assuming GOfile is in same directory as ppi file, replace .ppi extension with NCBI_to_GO
GOfile = options.infile[:-4]
GOfile += "_NCBI_to_GO"
node_labels = GOToDict(LMnodes, GOfile)

#for node in node_labels:
#    print "{}: ".format(node), node_labels[node]

# construct K = overlap matrix
num_nodes = len(node_labels)
K = np.zeros((num_nodes, num_nodes), dtype=int)

# if there is function overlap, set K[ij] and K[ji] to 1
for i in xrange(1,num_nodes):
    nd = LMnodes[i]
    for j in xrange(i+1,num_nodes):
        for label in node_labels[nd]:
            if label in node_labels[LMnodes[j]]:
                K[i][j] = 1
                K[j][i] = 1
                break

for i in range(num_nodes):
    print K[i]