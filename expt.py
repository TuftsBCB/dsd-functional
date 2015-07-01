#
# Just a container for functions we will probably be using often,
# and to separate plotting scripts from the more basic routines
#

import numpy as np
import sys
import dsd
import networkx as nx
import time
import os.path


def setup_graph(infile):
    """
    creates networkx graph and does necessary processing
    in order to generate matrices or do further work.
    Returns graph.
    """
    # get input graph as a networkx graph
    G = nx.read_adjlist(infile)

    # remove selfloop edges
    G.remove_edges_from(G.selfloop_edges())

    # we can only work with a connected graph
    # DONT Do this: dsd.make_connected(G)
    # instead, select the largest connected component
    G = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]

    return G


def get_LMset(G, nodeList, LMSetSize):
    """
    given a graph, its corresponding nodelist and LM set size, computes and returns LM set.
    """
    degree = np.array([G.degree(i) for i in nodeList],dtype='float')

    # sort nodes by degree
    degAscending = np.argsort(degree)
    degDescending = degAscending[::-1]

    # use the highest degree nodes as landmarks in each case
    # this is based on the assumption that high-degree nodes will be
    # easiest to match using sequence or other side information
    return degDescending[0:LMSetSize]


def dsd_matrix(G, nodeList, LMset, npyFile):
    # if npy path not entered, or file does not exist, compute D
    if not npyFile or not os.path.isfile(npyFile):
        #
        # construct hemat
        #
        HEmatrix = dsd.hematrix(np.array(nx.adjacency_matrix(G,nodeList).todense()))

        # construct DSD
        D = dsd.DSD(HEmatrix,LMset)

        if npyFile:
            np.save(npyFile, D)
    # otherwise just load and return it
    else:
        D = np.load(npyFile)

    return D


def GO_to_dict(nodeList, GOfile):
    """
    given set of nodes and name of GO file,
    return dictionary with node(gene id) as key, list of labels as value
    """
    nodeLabels = {}

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
            nodeLabels[part[0]] = labels

            # consider any node not listed in GO file as having no labels
            for node in nodeList:
                if node not in nodeLabels:
                    nodeLabels[node] = []

    return nodeLabels


def set_overlap(set1, set2):
    return not set1.isdisjoint(set2)


def overlap_matrix(nodeList, GOfile, npyFile, randomize=False):
    """
    given list of nodes and path to GO file(with labels),
    construct overlap matrix
    also added option to randomly shuffle label sets for comparison
    """
    # if npy path not entered, or file does not exist, compute K
    # if randomizing, DON'T LOAD and DON'T SAVE
    if not npyFile or not os.path.isfile(npyFile) or randomize:
        n = len(nodeList)
        # construct K = overlap matrix, dictionary of node:labels
        K = np.zeros((n, n), dtype=int)
        # get dictionary and reconstruct to only include nodes in nodeList
        nodeLabels = GO_to_dict(nodeList, GOfile)
        nodeLabels = {key: nodeLabels[key] for key in nodeList}

        if randomize:
            nodeLabels = dict(zip(np.random.permutation(nodeLabels.keys()), nodeLabels.values()))

        # if there is function overlap, set K[ij] to 1
        for i in xrange(n):
            for j in xrange(i+1,n):
                K[i][j] = K[j][i] = set_overlap(set(nodeLabels[nodeList[i]]), set(nodeLabels[nodeList[j]]))

        if npyFile and not randomize:
            np.save(npyFile, K)
    else:
        K = np.load(npyFile)

    return K