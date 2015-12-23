The purpose of this project is to demonstrate a simplified process I undertook for my master's thesis titled "Modelling the Out-of-Plane Behaviour of Paperboard" by replacing Paperboard with a fictitious material with simplified mechanics. Please note, to understand exactly what is going on, an understanding of Continuum Mechanics and Elasto-Plasticity is needed.

In order to demonstrate the process I have written a program in Matlab, and the same program Fortran, to show how to model elasto-plasticity using a simplified hardening and yield function. The example used here is similar to the von Mises yield function and demonstrates the use of the Backward-Euler Return-Mapping Algorithm. The hardening function is linear isotropic.

The material parameters are given as follows:

    - Young's modulus: E_const = 1000MPa
    - Shear modulus: G_const = 500MPa
    - Poisson's ratio: nu = 0:0
    - Yield stress: Sy0 = 1MPa
    - Hardening constant: H = 100MPa

The hardening function is linear isotropic and is given by:

    Sy = Sy0 - H * kap
where kap is the accumulated plastic strain.


The yield function is given by:

    f(Sigma, Sy) = norm(Sigma) - Sy

where Sigma is the Cauchy stress in the out of plane direction

To get an idea of the problem at hand, refer to Matlab/inputOutput.bmp, which shows the strain input and stress-strain output. The red crosses represent the 'hand' calcs. 

The list of files are as follows:

    - Matlab/input.m refers to the input file
    - Matlab/linearInterpolation.m interpolates between two values and is used to determine the input
    - Matlab/normCalcs.m shows the hand calcs for 1D compression/tension
    - Matlab/shearCalcs.m shows the hand calcs for 1D shear (not shown)
    - Matlab/plotout.m plots inputOutput.bmp
    - Matlab/symbolic.m calculates the gradients used in the return mapping algorithm
    - Matlab/umat.m is the material model (aka user material subroutine)

    - Fortran/input.90 is the same as input.m, except it also includes linear interpolation.
    - Fortran/umat.90 is the same as umat.m

umat.f90 can then be implemented into Abaqus (finite element software) in order to simulate elasto-plasticity for a fictitious material.
