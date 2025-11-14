# TRIPbasic

## System Requirements
- **Essential software and packages**
	- Python 3.12.5 with `numpy`, `pandas`, `networkx`, `graph-tool`, `community` (python-louvain), `torch`, `matplotlib`, and `seaborn`.
	- R 4.4.1 with `InteractionSet`, `igraph`, `tidyverse`, `shiny`, `ggplot2`, and `patchwork`.
	- FitHiChIP (loop caller; see the installation link below).
- **Language versions**
	- Python: 3.12.5
	- R: 4.4.1
- **Hardware**
	- For optimal performance, we recommend a computer with dual Intel(R) Xeon(R) Platinum CPUs (152 hardware threads, up to 3.6 GHz) and 1.0 TiB RAM to comfortably process the large sequencing data. A standard computer with enough RAM can be used to run ACA/TLP/TLPloops with the provided testdata. 

## Installation Guide
- Install the Python dependencies with your preferred package manager (pip or conda) and obtain the R packages via BiocManager or install.packages as appropriate.
- Follow the official FitHiChIP installation instructions: https://ay-lab.github.io/FitHiChIP/html/usage/installation.html.

## Instructions
0. Extract files from testdata.tar.gz into `testdata` directory.
1. To run the ACA Python notebook, replace ACA input path with the .aca file in the `testdata` directory or your project-specific data. A mcool file from actual data is also required.
2. To run the TLP identification workflow, execute `TLP/run_identify_loop_clusters.sh` with the 2D loop files swapped for the examples in `testdata/loops` or your own loops.
3. To call TRIPloops, run `bash filter_2d_loops.sh original_loop composite_ENCODE_peaks.bed 0.05`; demo FitHiChIP-called original loops and the `composite_sorted_unmerged.bed` reference live under `testdata` for quick validation.
