#Test of the Taylor integrator for the contact Hamiltonian h = p + q + S

include("../src/ContactIntegrator.jl")

using ContactIntegrator
using TaylorSeries
using FactCheck
import YAML

const order = 28

function field{T<:Real}(vec0::Array{T,1},c::T, beta::T)

  dim = length(vec0) #Dimension of the phase space

  vec0T = [ Taylor1([vec0[i]], order) for i in 1:dim]  #Promotion of the initial array to a Taylor series array

  y = [ Taylor1(0., order) for i=1:dim ]   #Auxiliar array
  D = [ Taylor1(0., order) for i=1:dim ]   #Auxiliar array

  for k = 0:order-1
    knext = k+1
    for i in 1:dim
      y[i] = Taylor1( vec0T[i].coeffs[1:k+1], k)  #Copy in the array "y" the coefficients of the Taylor series array vecO
    end

    #Field for the contact Hamiltonian h = p + q + w

    D[1] = 1.+0.*y[1]
    D[2] = -1. + y[2]
    D[3] = y[1] + y[3]

    #Update of the coefficients
    for i in 1:dim
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

parameters = YAML.load(open("testparameters.yaml"))

beta = parameters["beta"]
n = parameters["nsampling"]
deltat = parameters["deltatsample"]
c = parameters["c"]
initcond = Array{Float64}(3)


initcond = parameters["initcond"]
initcond = [initcond["q0"], initcond["p0"], initcond["S0"]]


t, q, p, S = contacthointegration!(field, n, initcond, deltat, c, beta)
###Solución analítica

tarray = collect(0.: deltat:(n-1)*deltat)
initcond = parameters["initcond"]
initcond = [initcond["q0"], initcond["p0"], initcond["S0"]]
qarray = tarray + initcond[1]
parray = 1. + exp(tarray)*(initcond[2] - 1.)
sarray = -initcond[1] - 1. - tarray + (initcond[3] + initcond[1] + 1.0)*exp(tarray)

relativetol = 1.0e-10 ##Chosen for Runge-Kutta

#These facts confirm that the local error in Taylor stays below or in the order of the
#local error in the Runge-Kutta method.

facts("Error (tolerance) Tests") do
  @fact (abs(sarray[end]) - abs(S[end]))/abs(sarray[end]) <= relativetol*10. --> true
  @fact (abs(parray[end]) - abs(p[end]))/abs(parray[end]) <= relativetol*10. --> true
  @fact (abs(qarray[end]) - abs(q[end]))/abs(qarray[end]) <= relativetol*10. --> true
end

#plot(tarray, q, "r--", tarray, qarray, "bs")
#plot(tarray, p, "r--", tarray, parray, "bs")
#plot(tarray, S, "r--", tarray, sarray, "bs")
