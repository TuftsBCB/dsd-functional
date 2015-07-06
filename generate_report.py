# script for generating part (hopefully most) of a LaTeX report.
# generates plots using routines from plotting.py, saves them,
# and generates the corresponding LaTeX markup
# LaTeX packages used: graphicx, float, subfig, geometry(margins set to .75in)
# run in command line with something like:
# python generate_report.py "PPIs and GO" "rat mouse" "tex.txt"

import plotting
import argparse
from matplotlib import pyplot as plt
import os

parser = argparse.ArgumentParser()

# path to folder containing ppi and GO files
# organisms is a string of organisms we want to generate plots for, separated by spaces
parser.add_argument("directory")
parser.add_argument("organisms")

# file to write LaTeX code to
parser.add_argument("outfile")

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

# 3 plots currently:    dsd against normalized overlap and pairs,
#                       pairs against summed overlap
#                       dsd against density
width = str(1.0/3)

# append, don't overwrite
with open(options.outfile, "a") as f:
    figNum = 1
    for organism in organisms:
        # setup
        print("generating plots for: " + organism)
        infile = options.directory + "/" + organism + ppiExt
        npyPath = npyDir + "/" + organism
        f.write("\\section{" + organism + "}\n")

        f.write("\\begin{figure}[H]\n\\centerline{\n")
        #################
        # dsd : density #
        #################
        f.write("\\subfloat[DSD vs Density]{\n")
        plt.figure(figNum)
        plotting.dsd_density(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy")

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_density"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1

        ########################
        # dsd : pairs, overlap #
        ########################
        f.write("\\subfloat[DSD vs Overlap, Pairs]{\n")
        plt.figure(figNum)
        plotting.dsd_overlap_pairs(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy")

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_overlap_pairs"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1

        ##############################
        # pairs : cumulative overlap #
        ##############################
        f.write("\\subfloat[Pairs vs Cumulative Overlap]{\n")
        plt.figure(figNum)
        plotting.pairs_summed_overlap(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy")

        # save plot
        figFile = plotsDir + "/" + organism + "_pairs_overlap"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1
        f.write("}\n\\end{figure}\n\n")

        # RANDOMIZE
        f.write("\\subsection{Plots based on random permutation of label sets}\n")
        f.write("\\begin{figure}[H]\n\\centerline{\n")
        #################
        # dsd : density #
        #################
        f.write("\\subfloat[DSD vs Density]{\n")
        plt.figure(figNum)
        plotting.dsd_density(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy", randomize=True)

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_density_randomized"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1

        ########################
        # dsd : pairs, overlap #
        ########################
        f.write("\\subfloat[DSD vs Overlap, Pairs]{\n")
        plt.figure(figNum)
        plotting.dsd_overlap_pairs(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy", randomize=True)

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_overlap_pairs_randomized"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1

        ##############################
        # pairs : cumulative overlap #
        ##############################
        f.write("\\subfloat[Pairs vs Cumulative Overlap]{\n")
        plt.figure(figNum)
        plotting.pairs_summed_overlap(infile, npyPath + "_dsd.npy", npyPath + "_overlap.npy", randomize=True)

        # save plot
        figFile = plotsDir + "/" + organism + "_pairs_overlap_randomized"
        plt.savefig(figFile)
        figFile += ".png"

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1
        f.write("}\n\\end{figure}\n\n")

        print("done!")
