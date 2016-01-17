#include("src/integrator.jl")
push!(LOAD_PATH,"src/")

import YAML
using ContactIntegrator
using HDF5

println("Type the output filename (without a format):")
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
nsimulations = parameters["nsimulations"]
initcond = Array{Float64}(3)

if nsimulations == 1
  try
    initcond = parameters["initcond"]
    initcond = [initcond["q0"], initcond["p0"], initcond["S0"]]
  catch
    initcond = [rand(), rand(), rand()]
  end
else
   initcond = [rand(), rand(), rand()]
end

file = h5open("./HDF5/$filename.hdf5", "w")

attrs(file)["nsampling"] = n
attrs(file)["deltatsample"] = deltat
attrs(file)["beta"] = beta
attrs(file)["c"] = c

for i in 1:nsimulations

  if i > 1
    initcond = [rand(), rand(), rand()]
  end

  file["simulation-$i/initcond"] = initcond

  results = contacthointegration!(n, initcond, deltat, c, beta)

  t = results[1]
  q = results[2]
  p = results[3]
  S = results[4]

  file["simulation-$i/p"] = p
  file["simulation-$i/q"] = q
  file["simulation-$i/t"] = t
  file["simulation-$i/S"] = S
end


close(file)

println("The simulation was succesfully done. \nOutput in /HDF5/$filename.hdf5")
