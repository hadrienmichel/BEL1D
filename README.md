# BEL1D

BEL1D is a program for the **stochastic** 1D imaging of the subsurface based on geophysical data.
It relies on the Bayesian Evidential Learning (BEL) framework called BEL 1D imaging (BEL1D in short) [(Michel et al., 2020)](https://doi.org/10.1016/j.cageo.2020.104456). The framework basically consists in 6 steps:
1. Generation of models from the prior model space and their associated data (forward modeling)
2. Reduction of the dimensionality of the dataset
3. Canonical correlation between the data and the models parameters
4. Extraction of the posterior distributions in the reduced model space
5. Sampling of models in reduced space
6. Back-transformation to original space

The different toolboxes that are inside this package are used in [(Michel et al., 2020)](https://doi.org/10.1016/j.cageo.2020.104456). They consist in 3 different toolboxes:
1. Interpretation of SNMR data (BEL1DSNMR)
2. Interpretation of dispersion curves from surface waves analyses (BEL1DSW)
3. General case (BEL1DGENERAL)

## Installation
To install the packages, just copy the main directory at the location that pleases you. In order to use the full functionalities of the codes, you **must** have:
1. A windows version (the MEX code for surface waves is compiled for Windows), or compiling the linux files and change the code to use the correct functions.
2. [MRSmatlab](https://doi.org/10.1190/geo2015-0461.1
) installed in your MATLAB environment (for the computation of sensitivity kernels in SNMR).

## Usage

To launch the GUIs, simply type ```BEL1D``` in the MATLAB command window. You will load a welcome GUI where you can choose the specific GUI that you are interested in.

For more details on the use of the toolboxes, refer to the user manual (Manual.pdf).

## Compatibility

The toolboxes have been coded under Matlab R2018a (9.4). If you encounter issues running the codes, it may be due to incompatibility of the Matlab version (functions not released yet, etc.).