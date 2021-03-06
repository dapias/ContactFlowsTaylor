using ODE
using HDF5
import YAML

include("field.jl")

println("Type the output filename (without a format):")
name = string(readline(STDIN))
filename = name[1:end-1]

while isfile("./HDF5/$filename.hdf5")
  println("The filename typed already exists in the HDF5 folder. Try another one:")
  filename = string(readline(STDIN))
end

parameters = YAML.load(open("../parameters.yaml"))

logisticparameter = parameters["c"]
temperatureparameter = parameters["beta"]

const c = logisticparameter
const beta = temperatureparameter

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

time = collect(0.0:deltat:(n-1)*deltat);

for i in 1:nsimulations

  if i > 1
    initcond = [rand(), rand(), rand()]
  end

  file["simulation-$i/initcond"] = initcond

  push!(initcond, 0.0)   ##Initial value of the extended variable that will be used to calculate the extended invariant

  #t, results = ode45(conthofield, initcond, time; points=:specified, reltol = 1.0e-10, abstol = 1.0e-14);
  t, results = ode45(conthofield, initcond, time; points=:specified, reltol = 1.0e-16, abstol = 1.0e-20);

  q = Float64[results[i][1] for i in 1:length(t)]
  p = Float64[results[i][2] for i in 1:length(t)]
  S = Float64[results[i][3] for i in 1:length(t)]
  invariant = Float64[results[i][4] for i in 1:length(t)]
  invariant = invariant + logrho(q,p, S, beta, c)


  file["simulation-$i/p"] = p
  file["simulation-$i/q"] = q
  file["simulation-$i/t"] = t
  file["simulation-$i/S"] = S
  file["simulation-$i/invariant"] = invariant
end


close(file)


println("The simulation was succesfully done. \nOutput in /HDF5/$filename.hdf5")



