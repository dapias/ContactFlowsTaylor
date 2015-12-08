using TaylorSeries


const ordenTaylor = 28
const epsAbs = 1.0e-20

# Integrador
function taylorStepper{T<:Number}(campo::Function, vec0::Array{T,1}, order::Int64, epsilon::T, beta::Float64)
  #Esta función recibe un arreglo correspondiente a los datos iniciales, resuelve las ecuaciones de movimiento y evalúa la solución en el tiempo óptimo (determinado mediantela función stepsize).

  n = length( vec0 )
  vec1T = campo( vec0, beta)

  # Step-size
  hh = Inf
  for i=1:n
    h1 = stepsize( vec1T[i], epsilon )
    hh = min( hh, h1 )
  end

  # Values at t0+h
#   for i=1:n
#     vec0[i] = evaluate( vec1T[i], hh )
#   end

  return hh, vec1T
end


# Returns the maximum step size from epsilon and the last two coefficients of the x-Taylor series
function stepsize{T<:Real}(x::Taylor1{T}, epsilon::Float64)
  ord = x.order
  h = Inf
  for k in [ord-1, ord]
    kinv = 1.0/k
    aux = abs( x.coeffs[k+1] )
    h = min(h, (epsilon/aux)^kinv)
  end
  return h
end
