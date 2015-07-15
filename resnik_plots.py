# script for generating dsd-resnik plots and LaTex markup
# generates plots using routines from plotting.py, saves them,
# and generates the corresponding LaTeX markup
# LaTeX packages used: graphicx, float, subfig, geometry(margins set to .75in)
# run in command line with something like:
# python resnik_plots.py "PPIs and GO" "rat mouse" "go.obo" "sprot.stripped" "tex.txt"

# these lines needed for working on ELF -- feel free to comment out
import sys
sys.path.append("/usr/lib64/python2.7/site-packages")
sys.path.append("/usr/lib/python2.7/site-packages/decorator-3.4.0-py2.7.egg")
sys.path.append("/usr/lib/python2.7/site-packages")

import plotting
import argparse
from matplotlib import pyplot as plt
import os
sys.path.append("./semsimcalc")
import semsimcalc

plt.rc('axes.formatter', use_mathtext=True)

parser = argparse.ArgumentParser()

# path to folder containing ppi and NCBI_to_GO files
# organisms is a string of organisms we want to generate plots for, separated by spaces
parser.add_argument("directory")
parser.add_argument("organisms")

# paths to gene ontology (.obo) and annotation corpus (.stripped) files needed to
# construct SemSimCalculator
parser.add_argument("oboFile")
parser.add_argument("acStrippedFile")

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

# initialize semsimcalc
calc = semsimcalc.SemSimCalculator(options.oboFile, options.acStrippedFile)

# for resnik scores, only plot dsd:density (normal and randomized)
width = str(0.5)

# append, don't overwrite
with open(options.outfile, "a") as f:
    figNum = 1
    f.write("\\section{Plots based on average pairwise Resnik scores}\n")
    for organism in organisms:
        # setup
        print("generating plots for: " + organism)
        infile = options.directory + "/" + organism + ppiExt
        npyPath = npyDir + "/" + organism

        f.write("\\begin{figure}[H]\n\\caption{Plots for " + organism + "}\n\\centerline{\n")
        #################
        # dsd : density #
        #################
        f.write("\\subfloat[DSD vs Density]{\n")
        plt.figure(figNum)
        plotting.dsd_density_res(infile, calc, npyPath + "_dsd.npy", npyPath + "_res.npy")

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_density_res"
        figFile += ".png"
        plt.savefig(figFile)

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1

        ##############################
        # dsd : density (randomized) #
        ##############################
        f.write("\\subfloat[DSD vs Density (with label sets randomized)]{\n")
        plt.figure(figNum)
        plotting.dsd_density_res(infile, calc, npyPath + "_dsd.npy", npyPath + "_res.npy", randomize=True)

        # save plot
        figFile = plotsDir + "/" + organism + "_dsd_density_res_randomized"
        figFile += ".png"
        plt.savefig(figFile)

        f.write("\\includegraphics[width=" + width + "\\textwidth]")
        f.write("{" + figFile + "}\n}\n")
        figNum += 1
        f.write("}\n\\end{figure}\n\n")

        print("done!")
