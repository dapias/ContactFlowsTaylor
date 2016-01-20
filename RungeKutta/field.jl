function contacthofield{T<:Real}(t::T, vec0::Array{T,1})  #The t is passed as parameter to use the function ode45 of the module
  y = vec0
  len = length(y)
  D = zeros(len)
  #Equations of motion generated by h = (exp(-beta H)/Z f(S))^(-1/2) with H = 1/2p^2 + 2q^2, Z = pi/beta and f(S) = exp(S-c)/(1 + exp(S-c))^2
  #The order in the array D (and y) is D[1] = q, D[2] = p, D[3] = S
  Z = pi/beta
  H = y[2]^2/2. + 2.*y[1]^2
  f = exp(y[3] - c)/(1 + exp(y[3] - c))^2.
  h = (exp(-beta*H)*f/Z)^(-0.5)

  D[1] = h/2.*beta*y[2]
  D[2] = h/2.*(-beta*4*y[1] + y[2]*(exp(y[3]-c) -1)/(exp(y[3]-c)+1))
  D[3] = h/2.*(-y[2]^2.*beta + 2.)

  return D
end