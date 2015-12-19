#Taylor series based integrator for contact hamiltonian flows

We provide the code that supports the numerical simulation reported in the manuscript [A thermostat algorithm generating target ensembles](http://arxiv.org/abs/1510.03942).

The code is organized as follows. 

The ``src`` folder contains the file  ``integrator.jl``. This file defines the module ContactIntegrator that exports the function ``contactintegration`` which is used as follows
```
contactintegration(nsteps, condinicial , timesampling, c, beta)
```

The function is called in the file ``script.jl`` which takes certain values for the parameters (that can be modified by the user), performs the numerical integration and generates a ``.hdf5`` file which is saved in the HDF5 folder with the name given by the user (asked by the script).

The folder ``notebooks`` contain a ``jupyter notebook`` which it is shown how to manipulate the HDF5 data to generate the figures 1-3 displayed on the paper.

###Authors

**Diego Tapias** (Facultad de Ciencias, UNAM)
**Alessandro Bravetti** (Instituto de Ciencias Nucleares, UNAM)

*2015.*







Programa para integrar flujos de contacto basado en el integrador de Taylor del paquete [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

