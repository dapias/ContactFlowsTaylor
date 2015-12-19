#Taylor series based integrator for contact hamiltonian flows

We provide the code that supports the numerical simulation reported in the manuscript [A thermostat algorithm generating target ensembles](http://arxiv.org/abs/1510.03942).

The code is based on the package  [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

It is organized as follows. 

The ``src`` folder contains the file  *ContactIntegrator.jl*. This file defines the module **ContactIntegrator** that exports the main function **contactintegration** that performs the core of the simulation.

The module is imported in the file *script.jl* which takes certain values for the parameters (see below), performs the numerical integration and generates a *.hdf5* file which is saved in the ``HDF5`` folder with the name given by the user (asked by the script).

The parameters reside in the file *parameters.yaml* and the users may modify it for their convenience.

Finally, the folder ``notebooks`` contain a *jupyter notebook* which illustrates how to manipulate the HDF5 data to generate the kind of figures displayed on the paper.

##Usage

Clone this repository with the following command in a UNIX terminal
```
~$ git clone https://github.com/dapias/ContactFlowsTaylor.git 
```

Then move into the created folder and execute the script.  To do that you may proceed in the third different following ways.

1. In a UNIX terminal type

 ```
 ~$ julia script.jl
 ```
2. In a unix terminal execute julia as
 ```
 ~$ julia
 ```
And then type the following command
 ```
 julia> include("script.jl")
 ```

3. Open a Jupyter Notebook and type in a cell
 ```
 include("script.jl")
 ```

### Requirements

#### General
Julia. It may be downloaded from its [webpage](http://julialang.org/downloads/)

#### Particular
The following packages are needed for the adequate execution of the program

- TaylorSeries
- HDF5
- YAML
- PyPlot (optional, it is used in the notebook)

To add a package inside Julia try the following
```
julia> Pkg.add("PackageName")
```
 
###Authors

**Diego Tapias** (Facultad de Ciencias, UNAM)

**Alessandro Bravetti** (Instituto de Ciencias Nucleares, UNAM)

*2015.*







