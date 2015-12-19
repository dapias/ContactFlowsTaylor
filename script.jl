include("src/integrator.jl")

import YAML
using ContactIntegrator

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

contactintegration(n, condinic, deltat, c, beta)
