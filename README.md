# dsd-functional
Experiments testing ability of DSD to capture functional relationships. Required (additional) modules include: NumPy, NetworkX, and matplotlib.

This repository's design revolves around 3 modules:
* `dsd.py`
  * contains basic functions for manipulating graphs and for obtaining various distance matrices -- you shouldn't need to call functions from this module directly.
* `expt.py`
  * serves as a wrapper for the dsd module, adding functions to compute matrices for functional overlap/similarity, as well as the option to specify numpy binary files for previously-computed matrices.
* `plotting.py`
  * contains routines that makes use of functions from the expt module, as well as numpy and matplotlib, to generate plots. Each function generates a single plot -- refer to each docstring for more information.

as well as `semsimcalc.py` (from [semsimcalc](https://github.com/TuftsBCB/semsimcalc)), which is used as a submodule for obtaining Resnik similarity scores.

Scripts making use of functions and routines from these modules are used to conduct our experiments. Currently, there are 2 such scripts: `overlap_plots.py` and `resnik_plots.py`. These scripts generate and save the necessary plots (as well as any computed (non-randomized) matrices), along with TeX markup, which you can see in `report.tex`.

Run these scripts from the command line as such:
* `python overlap_plots.py "PPIs and GO" "<organisms>" "tex.txt"`
* `python resnik_plots.py "PPIs and GO" "<organisms>" "go.obo" "sprot.stripped" "tex.txt"`

, replacing `<organisms>` with the organisms you wish to compute plots for (their corresponding files have to be present in the PPIs and GO directory), separated by spaces. `tex.txt` can be any text file; the TeX markup will be written to this file. Plots are saved to the `plots` directory, and numpy matrices are saved to the `NumPy files` directory.

Note: these scripts currently do not give you the option of *not* saving the computed matrices, which can be quite large. If this is a problem, feel free to raise an issue; alternatively, you can not include paths to numpy matrices when calling routines from the plotting module.
