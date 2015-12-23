The purpose of this project is to demonstrate the process of modelling elasto-plasticity using a simplified hardening and yield function. The example used here is similar to the von Mises yield function. The difference is, only the out-of-plane directions are considered. The hardening function is linear isotropic.

First the task is defined. This is followed by the implementation into MATLAB, which is split up into an input and a function. The results from MATLAB are then compared to the analytical solution for one-dimensional behaviour in tension-compression. Finally, the model is written in FORTRAN, which is the language ABAQUS uses to interpret material models. The results from ABAQUS are then compared to the results from MATLAB as well as the analytical solution.

TODO


