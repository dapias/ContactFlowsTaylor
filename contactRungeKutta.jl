include("src/RungeKsixthorder.jl")

function campoContacto(vector, beta)
  y = vector
  len = length(y)
  D = zeros(len)
    D[1] = (y[end]^2+1)*beta*y[2]*exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))/(2*y[end])
    D[2] = -exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))*(beta*y[1]*y[end]^3-2*y[2]*y[3]^2+beta*y[1]*y[3]+2*y[2])/(2*y[end]^2)
    D[end] = -(y[end]^2+1)*exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))*(beta*y[2]^2-2)/(2*y[end])
  return D
end


function contactIntegration(campo, nsteps, condinicial, deltat, beta=1.)

  t = 0.0
  x = condinicial

  q = Array(Float64, nsteps)
  p = Array(Float64, nsteps)
  S = Array(Float64, nsteps)
  tiempo = Array(Float64, nsteps)

  q[1] = x[1]
  p[1] = x[2]
  S[1] = x[3]
  tiempo[1] = t

  for i in 2:nsteps
          x = rungeK(x,campo,deltat, beta)
    q[i] = x[1]
    p[i] = x[2]
    S[i] = x[3]
    tiempo[i] = tiempo[i-1] + deltat
  end

  tiempo, q, p, S

end


