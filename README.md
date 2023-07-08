# 2D Finite Element Code for Structural Analysis
The purpose of this project is to use 3-noded triangular elements to construct a 2-D finite 
element code to determine nodal displacements, element stresses, and element strains for an arbitrary 
solid structure. The solid structure is a 2-D square plate with a circular hole in the center.

## Files
`hw10_struct_a.m` is the actual MATLAB code.  
`Structural Analysis Finite Element Code Report.pdf` is the report writeup.  

The report takes into account 4 different data sets, but I've only uploaded one data set, the R set:  

`displacementsR.txt` contains data for the initial displacements in the structure.  
`forcesR.txt` contains data for the external forces on the structure.  
`nodesR.txt` contains data about the locations of each node.  
`elementsR.txt` contains data about the 3-noded elements, including the number of elements, the value for Young’s modulus, and the value for 
Poisson’s ratio.
