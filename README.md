# RNA-sequencing course HS2023

## Module: Ribosome profiling

[![GitHub tag](https://img.shields.io/github/tag/JeremyRotzetter/ribosome-profiling?include_prereleases=&sort=semver&color=blue)](https://github.com/JeremyRotzetter/ribosome-profiling/releases/)
[![License](https://img.shields.io/badge/License-GPLv3-blue)](#license "View license summary")
[![Issues - ribosome-profiling](https://img.shields.io/github/issues/JeremyRotzetter/ribosome-profiling)](https://github.com/JeremyRotzetter/ribosome-profiling/issues "View open issues")
[![Made with R](https://img.shields.io/badge/R-4.3.2-blue?logo=r&logoColor=white)](https://cran.r-project.org/ "Go to CRAN homepage")
[![Made with R](https://img.shields.io/badge/RStudio-2023.09.1_Build_494-blue?logo=rstudio&logoColor=white)](https://posit.co/products/open-source/rstudio/ "Go to RSTUDIO IDE homepage")
[![Shell Script](https://img.shields.io/badge/Shell_Script-blue?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/ "Go to Bash homepage")

The code and scripts contained within this repository were written during the fall semester 2023 for the *RNA-sequencing* course at the University of Bern.

## Setup for R
In order to use the functionalities provided by the R source code here, some preliminary steps might have to be taken. Please note that, since the project is based inside an R project, any considered path within Rscripts is relative to the working directory of the project itself and needs to be changed if used otherwise.

1. Install [R 4.3.2](https://www.r-project.org/) for your corresponding OS.
2. Install [RStudio](https://posit.co/download/rstudio-desktop/).
3. Download the repository (either download and extract zip or git clone).
4. Open project file [ribosome-profiling.Rproj](https://github.com/JeremyRotzetter/ribosome-profiling/blob/main/ribosome-profiling.Rproj). renv should automatically be installed and activated. If not `source(renv\activate.R)`.
5. renv will now check the *renv.lock* file and ask if it should install all dependencies of the project. Select `Y`. renv will then download and install required R packages from their respective sources and in the correct version into the project library.
6. Should the error message `Error in system(paste(MAKE, p1(paste("-f", shQuote(makefiles))), "compilers"),  : 
  'make' not found` appear during the restoration of the project library, then this means that one of the corresponding package
versions has been archived and the binary is no longer available. In this case, the package must be compiled and built from
source for the respective OS. The 'make' utility is commonly used in R packages to compile and build software from source code.
Since the 'make' utility is not installed on your system, there are different ways to install it, depending on your OS:

    - For Windows: Windows users should install the appropriate iteration of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for their R version. Rtools is a toolchain bundle used for building R packages from source (those that need compilation of C/C++ or Fortran code) and for build R itself.

    - For Linux: You can install 'make' using the package manager specific to your Linux distribution. For example, on Ubuntu or Debian-based systems, you can run the following command in the terminal: `apt-get install build-essential`

    - For macOS: In order to compile R for macOS, you will need Xcode developer tools from Apple and GNU Fortran compiler. More information can be found on the [r-project website for mac](https://mac.r-project.org/tools/).

    - Alternatively it might be possible to also install binaries from [Posit Public Package Manager](https://packagemanager.posit.co/) snapshots, which often has binary packages available.

## Repository structure

- [renv](renv) folder: contains renv and the files needed to initialize the project library.
- [scripts](scripts) folder: contains the scripts used in the analysis of the data. The logical order in which they should be executed is given by their numbering, with support scripts being unnumbered.
    >[!TIP]
    > All bash scripts were written with SBATCH options so that the script could be submitted to the SLURM workload manager to reduce the computational load at the head node of the IBU cluster (or any other HPC).
- [renv.lock](renv.lock) file: this file lists all dependencies used in the project, i.e. the R version, all packages, their version and their source (CRAN, Bioconductor, GitHub). renv uses this file to restore the library.
- [ribosome-profiling.Rproj](ribosome-profiling.Rproj) file: the project can be opened with this .Rproj file in RStudio, and the working directory is automatically set to the project directory where the .Rproj file is located. As a consequence, any path found within the scripts is relative to this project working directory.

## Acknowledgements
Raw RNA-seq data as well as guidance, helpful tips/hints and suggestions for the coding part were kindly provided by Dr. Puneet Sharma.

## License
Released under [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/) by [@JeremyRotzetter](https://github.com/JeremyRotzetter).

This license means:
- You can freely copy, modify, distribute and reuse this software.
- The ***original license*** must be included with copies of this software.
- Please _link back_ to this repo if you use a significant portion of the source code.
- The software is provided "as is", without warranty of any kind.
- Source code must be made available when this software is distributed.