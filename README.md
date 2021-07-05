# lunar1Dheat

This repository contains a variety of source code used to model heat transfer in the upper 2 meters of the lunar surface. 

Each folder 


The codes are split up such that the one-dimensional heat flow equation may be solved with standard thermophysical parameters, or 




The codes presented here numerically solve the one-dimensional heat flow equation using a forward differencing technique such that time-dependent temperature variations can be determined at the surface and at depth. A complete description of this technique can be found in Appendix A of Hayne et al. (2017). The standard models presented here contain the thermophysical parameters and models described in Appendix A of Hayne et al. (2017), however, we incorporate the bolometric Bond albedo model developed in Feng et al. (2020). Though this repository fully incorporates the Feng albedo model, the original source code used in Feng et al. (2020) may be downloaded from https://zenodo.org/record/3579654. 

The updated models presented here incorporate 
