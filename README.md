### Modelling the Out-of-Plane Behaviour of Paperboard

The purpose of this project is to demonstrate a simplified process I undertook for my master's thesis by replacing paperboard with a fictitious material with simplified mechanics. Knowledge is expected in Continuum and Material Mechanics in order to understand exactly what is going on.

In order to demonstrate the process I have written a simplified elasto-plastic model in Matlab and in Fortran, which takes strain as input and determines the resulting stress. The example used here is similar to the von Mises yield function and demonstrates the use of the Backward-Euler Return-Mapping Algorithm. The hardening function is assumed to be linear isotropic.

The material parameters are given as follows:

- Young's Modulus: `E_const = 1000MPa`
- Shear Modulus: `G_const = 500MPa`
- Poisson's Ratio: `nu = 0.0`
- Yield Stress: `Sy0 = 1MPa`
- Hardening Constant: `H = 100MPa`
- Cauchy Stress components: `sigma_11 = 0`, `sigma_12 = 0`, `sigma_22 = 0`

The hardening function is assumed to be linear isotropic and is given by:

    Sy = Sy0 - H * kap
    
where `kap` is the accumulated plastic strain.


The yield function is assumed to be:

    f(sigma, Sy) = norm(sigma) - Sy

The list of files are as follows:

- `Matlab/inputOutput.bmp` shows the assumed strain input and resulting stress-strain output
- `Matlab/input.m` determines the linear strain input shown in `Matlab/inputOutput.bmp`
- `Matlab/linearInterpolation.m` splits the linear strain input up into smaller intervals
- `Matlab/normCalcs.m` analytically calculates 1D compression/tension output, shown as red crosses in `Matlab/inputOutput.bmp`
- `Matlab/shearCalcs.m` analytically calculates 1D shear output (not used)
- `Matlab/plotout.m` plots `Matlab/inputOutput.bmp`
- `Matlab/symbolic.m` calculates the gradients used in the Return-Mapping Algorithm
- `Matlab/umat.m` demonstrates the material model (a.k.a. the user material subroutine)  
- `Fortran/input.90` is the same as input.m (linear interpolation is included)
- `Fortran/umat.90` is the same as umat.m

`Fortran/umat.90` can then be implemented into finite element software (eg. Abaqus) as a user material subroutine (UMAT) in order to simulate elasto-plasticity for this material.
