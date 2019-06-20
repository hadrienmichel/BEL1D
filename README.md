# BEL1D

BEL1D is a program for the **stochastic** 1D imaging of the subsurface based on geophysical data.
It relies on the Bayesian Evidential Learning (BEL) framework called BEL 1D imaging (BEL1D in short) [(Michel et al., 2019a) - not yet available](http://www.link_to_the_paper.be). The framework basically consists in 6 steps:
1. Generation of models from the prior model space and their associated data (forward modeling)
2. Reduction of the dimensionality of the dataset
3. Canonical correlation between the data and the models parameters
4. Extraction of the posterior distributions in the reduced model space
5. Sampling of models in reduced space
6. Back-transformation to original space

The different toolboxes that are inside this package are used in [(Michel et al., 2019b) - not yet available](http://www.link_to_the_paper.be). They consist in 3 different toolboxes:
1. Interpretation of SNMR data (BEL1DSNMR)
2. Interpretation of dispersion curves from surface waves analyses (BEL1DSW)
3. General case (BEL1DGENERAL)

## Installation
To install the packages, just copy the main directory at the location that pleases you. In order to use the full functionalities of the codes, you **must** have:
1. [Geopsy](http://www.geopsy.org) installed on your machine and [**added to your PATH**](https://www.howtogeek.com/118594/how-to-edit-your-system-path-for-easy-command-line-access/), to access the gpdc command line tool for the forward modeling of dispersion curves.
2. [MRSmatlab](https://doi.org/10.1190/geo2015-0461.1
) installed in your MATLAB environment (for the computation of sensitivity kernels in SNMR).

## Usage

To launch the GUIs, simply type ```BEL1D``` in the MATLAB command window. You will load a welcome GUI where you can choose the specific GUI that you are interested in.

For more details on the use of the toolboxes, refer to the user manual (Manual.pdf).
