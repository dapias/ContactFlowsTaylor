#include("src/integrator.jl")
push!(LOAD_PATH,"src/")

import YAML
using ContactIntegrator
using HDF5

println("type the output filename:")
name = string(readline(STDIN))
filename = name[1:end-1]

while isfile("./HDF5/$filename.hdf5")
    println("The filename typed already exists in the HDF5 folder. Try another one:")
    filename = string(readline(STDIN))
end

parameters = YAML.load(open("parameters.yaml"))

c = parameters["c"]
beta = parameters["beta"]
n = parameters["nsampling"]
deltat = parameters["deltatsample"]
condinic = Array{Float64}(3)

try 
    condinic = parameters["condinic"]
    condinic = [condinic["q0"], condinic["p0"], condinic["S0"]] 
catch
    condinic = [rand(), rand(), rand()]
end

file = h5open("./HDF5/$filename.hdf5", "w")

attrs(file)["condinicial"] = condinic
attrs(file)["nsampling"] = n
attrs(file)["deltatsample"] = deltat
attrs(file)["beta"] = beta
attrs(file)["c"] = c

results = contactintegration(n, condinic, deltat, c, beta)
t = results[1]
q = results[2]
p = results[3]
S = results[4]

file["simulation-1/p"] = p
file["simulation-1/q"] = q
file["simulation-1/t"] = t
file["simulation-1/S"] = S

close(file)

println("The simulation was succesfully done. \nOutput in /HDF5/$filename.hdf5")
