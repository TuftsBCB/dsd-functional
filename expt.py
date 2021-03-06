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
sys.path.append("./semsimcalc")
import semsimcalc


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


def get_LMset(G, nodeList, LMsetSize, centralityFunc='degree', **kwargs):
    """
    given a graph, its corresponding nodelist and LM set size, computes and returns LM set.
    """
    centralityFuncs = {'degree': nx.degree_centrality, \
                       'closeness': nx.closeness_centrality, \
                       'betweenness': nx.betweenness_centrality, \
                       'eigenvector': nx.eigenvector_centrality_numpy, \
                       'katz': nx.katz_centrality_numpy}

    f = centralityFunc
    if type(f) == str:
        if f == 'random':
            return np.random.choice(len(nodeList), LMsetSize, replace=False)
        else:
            f = centralityFuncs[f]

    centrality = f(G, **kwargs)
    scores = np.array([centrality[node] for node in nodeList])

    # sort nodes by score
    scoreAscending = np.argsort(scores)
    scoreDescending = scoreAscending[::-1]

    return scoreDescending[0:LMsetSize]


def dsd_matrix(G, nodeList, npyFile, LMsetSize=50, centralityFunc='degree', **kwargs):
    """
    any kwargs, if specified, will be passed into centrality function call.
    """
    # if npy path not entered, or file does not exist, compute D
    if not npyFile or not os.path.isfile(npyFile):
        #
        # construct hemat
        #
        adjMatrix = np.array(nx.adjacency_matrix(G,nodeList))
        if np.shape(adjMatrix) == ():
            adjMatrix = np.array(nx.adjacency_matrix(G,nodeList).todense())
        HEmatrix = dsd.hematrix(adjMatrix)

        # construct DSD
        LMset = get_LMset(G, nodeList, LMsetSize, centralityFunc, **kwargs)
        D = dsd.DSD(HEmatrix,LMset)

        if npyFile:
            try:
                np.save(npyFile, D)
            except IOError:
                os.makedirs(npyFile[:npyFile.rfind('/')])
                np.save(npyFile, D)
    # otherwise just load and return it
    else:
        D = np.load(npyFile)
    return D


def sp_matrix(G, nodeList, npyFile):
    # same as dsd_matrix, but for shortest-path distance
    if not npyFile or not os.path.isfile(npyFile):
        S = dsd.spmatrix(G, nx.shortest_path(G), nodeList)

        if npyFile:
            np.save(npyFile, S)
    else:
        S = np.load(npyFile)

    return S


def diffusion_matrix(G, nodeList, npyFile, gamma=8):
    """
    compute inverse of Laplacian matrix and symmetrize. (default gamma is arbitrary)
    the higher the influence, the closer we expect function to be, so let distances be reciprocal.
    """
    if not npyFile or not os.path.isfile(npyFile):

        L = np.array(nx.laplacian_matrix(G,nodeList))
        # depending on version of networkx, might get sparse matrix instead. if so, do this:
        if np.shape(L) == ():
            L = np.array(nx.laplacian_matrix(G,nodeList).todense())
        m, n = np.shape(L)
        L = L + (np.eye(m, n)*gamma)
        D = np.linalg.inv(L)
        n = len(nodeList)
        for i in xrange(n):
            for j in xrange(i+1,n):
                D[i][j] = D[j][i] = 1/(min(D[i][j], D[j][i]))

        if npyFile:
            np.save(npyFile, D)
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
            #nodeLabels = dict(zip(np.random.permutation(nodeLabels.keys()), nodeLabels.values()))
            keys = nodeLabels.keys()
            np.random.shuffle(keys)
            nodeLabels = dict(zip(keys, nodeLabels.values()))

        # if there is function overlap, set K[ij] to 1
        for i in xrange(n):
            for j in xrange(i+1,n):
                K[i][j] = K[j][i] = set_overlap(set(nodeLabels[nodeList[i]]), set(nodeLabels[nodeList[j]]))

        if npyFile and not randomize:
            np.save(npyFile, K)
    else:
        K = np.load(npyFile)

    return K


def triu_ravel(A):
    """
    Given a matrix, returns the upper triangle (not including main diagonal)
    as a 1D array. Row-major.
    """
    m, n = np.shape(A)
    flat = np.array([])
    for i in range(m):
        if i >= n-1:
            break
        flat = np.append(flat, A[i][i+1:])
    return flat


def resnik_matrix(nodeList, GOfile, calc, npyFile, randomize=False):
    """
    like overlap matrix, this time using semsimcalc's resnik similarity measure
    calc refers to SemSimCalculator instance
    """
    if not npyFile or not os.path.isfile(npyFile) or randomize:
        n = len(nodeList)
        # construct K = overlap matrix, dictionary of node:labels
        R = np.zeros((n, n), dtype=float)
        # get dictionary and reconstruct to only include nodes in nodeList
        nodeLabels = GO_to_dict(nodeList, GOfile)
        nodeLabels = {key: nodeLabels[key] for key in nodeList}

        if randomize:
            #nodeLabels = dict(zip(np.random.permutation(nodeLabels.keys()), nodeLabels.values()))
            keys = nodeLabels.keys()
            np.random.shuffle(keys)
            nodeLabels = dict(zip(keys, nodeLabels.values()))

        empty_count = 0
        scoreless_count = 0
        for i in xrange(n):
            for j in xrange(i+1,n):
                lefts = nodeLabels[nodeList[i]]
                rights = nodeLabels[nodeList[j]]
                if not lefts or not rights:
                    R[i][j] = R[j][i] = 0
                    empty_count += 1
                else:
                    score = calc.pairwise_average_term_comp(lefts, rights, calc.simRes)
                    if score == None:
                        R[i][j] = R[j][i] = 0
                        scoreless_count += 1
                    else:
                        R[i][j] = R[j][i] = score

        print "number of pairs comprising protein(s) with empty label set(s): {}".format(empty_count)
        print "number of non-empty pairs without proper resnik score: {}".format(scoreless_count)

        if npyFile and not randomize:
            np.save(npyFile, R)
    else:
        R = np.load(npyFile)

    return R
