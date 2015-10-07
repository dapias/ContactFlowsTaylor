#using ForwardDiff
using TaylorSeries

function fofS(S::Float64)
  f = (S + 1./S)^2
end

function potential(q::Vector, w = 1)  #Adecuar parámetros según el problema.
  V = sum(1/2*w^2*q.^2)
end

function kinetic(p::Vector, m = 1)
  K = sum(1/(2*m)*p.^2)
end

function shamiltonian(y::Vector, n::Int64, w = 1, m = 1)
  H = potential(y[1:n]) + kinetic(y[n+1:end])
end

function rho_target(H, beta = 1)   #¿Dónde entra la normalización?
  rho = exp(-beta*H)
end


function chamiltonian(x::Vector)
  f = fofS(x[end])
  n = iceil((length(x) - 1)/2)
  H = shamiltonian(x[1:end-1], n)
  rho = rho_target(H)
  h = (f/rho)^(1/(n+1))
end

#gradofH = ForwardDiff.gradient(chamiltonian)


convert(q::Type{FloatingPoint}, p::Taylor1{Float64}) = p = Taylor1{q}
function fofS(S::Taylor1{Float64})
  f = (S + 1./S)^2
end


function potential(q::Array{Taylor1{Float64},1}, w = 1)  #Adecuar parámetros según el problema.
  V = sum([1/2*w^2*(q[i])^2 for i in 1:length(q)])
end

function kinetic(p::Array{Taylor1{Float64},1}, m = 1)
  K = sum([1/(2*m)*(p[i])^2 for i in 1:length(p)])
end

#####Si los grados de libertad son 1
function potential(q::Taylor1{Float64}, w = 1)  #Adecuar parámetros según el problema.
  V = 1/2*w^2*q^2
end

function kinetic(p::Taylor1{Float64}, m = 1)
  K = 1/(2*m)*p^2
end

#################

function shamiltonian(y::Array{Taylor1{Float64},1}, n::Int64, w = 1, m = 1)
  H = potential(y[1:n]) + kinetic(y[n+1:end])
end

function rho_target(H, beta = 1)   #¿Dónde entra la normalización?
  rho = exp(-beta*H)
end


function chamiltonian(x::Array{Taylor1{Float64},1})
  f = fofS(x[end])
  n = iceil((length(x) - 1)/2)
  H = shamiltonian(x[1:end-1], n)
  rho = rho_target(H)
  h = (f/rho)^(1/(n+1))
end


