# lunar1Dheat

This repository contains a variety of source code used to model heat transfer in the upper 2 meters of the lunar surface, and it was used to produce the results described in "A Global Thermal Conductivity Model for Lunar Regolith at Low Temperatures" by A. Martinez and M. Siegler (2021). The codes presented here numerically solve the one-dimensional heat flow equation using a forward differencing technique such that time-dependent temperature variations may be determined at the surface and at depth. A complete description of this technique may be found in Appendix A of Hayne et al. (2017). 

The repository is organized such that each folder (excluding `1DFunctions`) contains either a standard thermal model or an updated thermal model. The standard thermal models contain the thermophysical properties described in Feng et al. (2020). Though this repository fully incorporates the Feng model (with a couple of minor modifications to the grid parameters), original source code used in Feng et al. (2020) may be downloaded from https://zenodo.org/record/3579654. The updated thermal models incorporate a new low-temperature model for the thermal conductivity of lunar regolith. These folders contain the core programs to execute. 

`Global`  - Contains the main model used to model global lunar surface temperatures. Valid for flat and relatively rock-free regions. 

Execute [temperature,P,totalsteps,z] = heat1D(latitude,surface) in the command line to run the thermal model with the updated thermal conductivity model.

Execute [temperature,P,totalsteps,z,D] = heat1DcraterStandard(latitude,D) to run the thermal model with standard thermophysical parameters. 

`Craters` - Models temperatures in near-polar bowl-shaped craters  
(Similar execution as global model, but substitute the appriopriate function) 

`PSRShoemaker` - Models temperature in a permanently shadowed region located in Shoemaker crater 

`1DFunctions` - Contains all the functions needed for core programs 
