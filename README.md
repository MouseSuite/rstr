# README
rstr: Rodent Statistics Toolbox in R

Copyright (C) 2025 The Regents of the University of California

Created by Shantanu H. Joshi, Yeun Kim, Kayla A. Schroeder, and David W. Shattuck

rstr is licensed under an GPLv2-only license (https://spdx.org/licenses/GPL-2.0-only.html).
Please see the enclosed LICENSE file for more details.

---

The Rodent Statistics toolbox in R (rstr) is a software package developed in R that performs statistical analysis of population-level rodent brain imaging data processed using [MouseSuite](https://mousesuite.org/). Specifically, it provides statistical tools for conducting tensor based morphometry, and analysis of gray matter volume and diffusion measures.

## rstr Installation
For more detailed installation instructions, usage examples, or to check for updated versions of rstr, please visit the rstr website: http://mousesuite.org/rstr/.

### Prerequisites
* rstr uses tools from BrainSuite for rendering statistical maps and showing 3D visualizations. Ensure BrainSuite is installed on your computer: <http://brainsuite.org/download> _Note: you do not need to install the MATLAB Runtime to use the BrainSuite tools with rstr_.
* Ensure that R is installed on your computer: https://cran.r-project.org/
* Install RStudio: https://www.rstudio.com/products/rstudio/#Desktop
* Windows only: Install Rtools, available at https://cran.r-project.org/bin/windows/Rtools/
* Install the R package [remotes](https://remotes.r-lib.org) - `install.packages('remotes')`

### Install from GitHub 
* Open RStudio and enter the following command: `remotes::install_github("MouseSuite/rstr")`to install the latest **rstr** version. 


### Download a specific version from [GitHub](https://github.com/MouseSuite/rstr) and install from a local directory
```
remotes::install_local('/path/to/rstr-main.zip')
```

### Windows - Install from GitHub
Note that on Windows, you will need to use double backslashes (\\) in the path because backslash is an escape character. You can also replace the backslashes with forward slashes.
Replace `C:\\path\\to\\` or `C:/path/to/` with the path to the folder where you downloaded rstr.


### Check your installation
Verify that rstr is installed by typing:
```
library(rstr)
```
This will load the library into R. Then test that rstr found your local installation of BrainSuite:
```
get_brainsuite_install_path()
```
This should display the BrainSuite installation path.


## Methods
Rstr performs statistical analysis on the outputs of the MouseSuite structural and diffusion workflows, which perform volumetric registration (rodentreg), and, optionally, processing of diffusion MRI data using the MouseSuite diffusion pipeline. Rstr is used to perform population-level statistical analysis of various neuroimaging measures produced by these components. Statistical analysis of voxel-wise data is performed in the common coordinate space of the atlas by resampling the data from subject coordinates to a the chosen atlas space.

Rstr supports the following analysis methods:

* tensor based morphometry (TBM) analysis of voxel-wise magnitudes of the 3D deformation fields of MRI images registered to the atlas
* region of interest (ROI)-based analysis of gray matter volume within cortical ROIs
* diffusion parameter maps analysis (DBA) of fractional anisotropy, mean diffusivity, radial diffusivity
* correction for multiple comparisons using false discovery rate (FDR) or permutation testing methods

Rstr is cross-platform and is available on macOS, Windows,and Linux based systems (all platforms with R support). Rstr is distributed under an open source license (GPLv2-only). Rstr  supports functionality for automated report generation to visualize statistical results using R-shiny and R markdown. The volumetric analysis report contains the cluster table, visualizations of clusters on image slices, and shows both the unadjusted and the adjusted versions of p-values and t statistics, respectively. The ROI analysis report shows the demographic spreadsheet, automatic bar plots for ANOVA and regressions, and scatter plot for correlation analyses. Rstr also exports an R markdown report that contains reproducible R commands in both the Rmd file and in the html document. This enables complete reproducibility of statistical results and only requires packaging the R markdown file along with the data.

---
