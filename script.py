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
nSamps = 0

# get input graph as a networkx graph
G = nx.read_adjlist(options.infile)

# remove selfloop edges
G.remove_edges_from(G.selfloop_edges())

# make connected if not so
dsd.make_connected(G)

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
nonLMset = [i for i in xrange(n) if i not in LMset]
if (nSamps != 0):
    sampleSet = random.sample(nonLMset,nSamps)
else:
    # modified for testing ALL nodes
    sampleSet = nonLMset
    nSamps = len(sampleSet)
#
#
# construct hemat
#
HEmatrix = dsd.hematrix(np.array(nx.adjacency_matrix(G,nodeList).todense()))

# construct DSD
D = dsd.DSD(HEmatrix,LMset)


