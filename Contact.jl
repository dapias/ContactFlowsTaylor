using TaylorSeries


const ordenTaylor = 28
const epsAbs = 1.0e-20

#n será el número de grados de libertad. La dimensión del espacio fase será
# 2n + 1.

function generarics(n)
  ics = zeros(2n+1)
  #Primeras n variables = q^i
  #Segundas n variables = p_i
  #Última variable = S
  for i in 1:length(ics)
    ics[i] = 1.0 ##Es algo que tenemos que cambiar para mostrar que la ergodicidad no depende de las condiciones iniciales.
  end
  ics
end

#Tengo la libertad de elegit f(S) y p_target(p_i, q^i)


# Integrador
function taylorStepper{T<:Number}( jetEqs::Function, vec0::Array{T,1}, order::Int64, epsilon::T)
  #Esta función recibe un arreglo correspondiente a los datos iniciales, resuelve las ecuaciones de movimiento y evalúa la solución en el tiempo óptimo (determinado mediantela función stepsize).

  n = length( vec0 )
  vec1T = jetEqs( vec0 )

  # Step-size
  hh = Inf
  for i=1:n
    h1 = stepsize( vec1T[i], epsilon )
    hh = min( hh, h1 )
  end

  # Values at t0+h
  for i=1:n
    vec0[i] = evaluate( vec1T[i], hh )
  end

  return hh, vec0
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





function campoContacto{T<:Real}(x ::Array{T,1})
  #Esta función recibe un arreglo, lo convierte en uno tipo Taylor, resuelve las ecuaciones de movimiento y regresa los coeficientes de la solución.
  order = ordenTaylor
  dim = length(x) #Dimension of the phase space
  vec0T = [ Taylor1([x[i]], order) for i in 1:dim]
  n = iceil((dim - 1)/2)

  y = [ Taylor1(0., order) for i=1:dim ]
  D = [ Taylor1(0., order) for i=1:dim ]


  for k = 0:order-1
    knext = k+1

    for i in 1:dim
      y[i] = Taylor1( vec0T[i].coeffs[1:k+1], k)
    end

    #Acá entra el Hamiltoniano (parece que si lo obtengo numéricamente con ForwardDiff no funciona el método de Taylor)
    # porque necesito una expresión analítica que me permita hacer las iteraciones, esto parece ser un problema grave
    # para las aplicaciones numéricas del método de Taylor


    #Ecuaciones de movimiento
    D[end] = -(y[end]^2+1)*exp((1/4)*y[2]^2+(1/4)*y[1]^2)*(y[2]^2-2)/(2*y[end])
    D[1] = (y[end]^2+1)*y[2]*exp((1/4)*y[2]^2+(1/4)*y[1]^2)/(2*y[end])
    D[2] = exp((1/4)*y[2]^2+(1/4)*y[1]^2)*(-q*y[end]^3+2*y[2]*y[3]^2-y[1]*y[3]-2*y[2])/(2*y[end]^2)


    #Actualización de los coeficientes
    for i in 1:dim
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

function contactIntegration(n, campo, nsteps)
  t0 = 0.0
  x = generarics(n)

  q = Array(Float64, nsteps)
  p = Array(Float64, nsteps)
  S = Array(Float64, nsteps)
  tiempo = Array(Float64, nsteps)

  q[1] = x[1]
  p[1] = x[2]
  S[1] = x[3]
  tiempo[1] = t0

  for i in 2:nsteps
    q[i] = x[1]
    p[i] = x[2]
    S[i] = x[3]
    tiempo[i] = t + tiempo[i-1]
    t, x = taylorStepper(campo, x, ordenTaylor, epsAbs )
  end

  t, q, p, S

end

x = contactIntegration(1, campoContacto, 100) ###Después se grafica