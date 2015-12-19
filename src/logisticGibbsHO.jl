module ContactIntegrator

export contactintegration

using TaylorSeries

const ordenTaylor = 28
const tolerance = 1.0e-20

# Integrador
function taylorStepper{T<:Real}(campo::Function, vec0::Array{T,1}, order::Int64, epsilon::T, c::Float64, beta::Float64)
  #Esta función recibe un arreglo correspondiente a los datos iniciales, resuelve las ecuaciones de movimiento y reporta el time-step (hh) y el Taylor solución

  #n = length( vec0 )
  #n = 3
  vec1T = campo( vec0, c, beta)

  # Step-size
  hh = Inf
  for i in eachindex(vec0)
    @inbounds h1 = stepsize( vec1T[i], epsilon )
    hh = min( hh, h1 )
  end

  # Values at t0+h
#   for i=1:n
#     vec0[i] = evaluate( vec1T[i], hh )
#   end

  return hh, vec1T
end


# Returns the maximum step size from epsilon and the last two coefficients of the x-Taylor series
function stepsize{T<:Real}(x::Taylor1{T}, tolerance::Float64)
  ord = x.order
  h = Inf
  for k in [ord-1, ord]
    kinv = 1.0/k
    aux = abs( x.coeffs[k+1] )
    h = min(h, (tolerance/aux)^kinv)
  end
  return h
end


function campoContacto{T<:Real}(vec0::Array{T,1},c::T, beta::T)
  #Esta función recibe un arreglo, lo convierte en uno tipo Taylor, resuelve las ecuaciones de movimiento y regresa los coeficientes de la solución.
  order = ordenTaylor

  dim = length(vec0) #Dimension of the phase space
  #  dim = 3
  vec0T = [ Taylor1([vec0[i]], order) for i in 1:dim]


  y = [ Taylor1(0., order) for i=1:dim ]
  D = [ Taylor1(0., order) for i=1:dim ]


  for k = 0:order-1
    knext = k+1

    for i in 1:dim
      y[i] = Taylor1( vec0T[i].coeffs[1:k+1], k)   #Copia en el arreglo "y" los coeficientes de la serie de Taylor hasta lo que va la iteración
    end

    #Acá entra el Hamiltoniano (parece que si lo obtengo numéricamente con ForwardDiff no funciona el método de Taylor)
    # porque necesito una expresión analítica que me permita hacer las iteraciones, esto parece ser un problema grave
    # para las aplicaciones numéricas del método de Taylor

    #Las entradas *D* y *y* estan en el orden q, p, S

    #Ecuaciones de movimiento

    #Campo de Gibbs de contacto con "c" libre (Logistico)

    D[1] = 1/2*pi^(1/2)*beta/(exp(y[3])/exp(c)*exp(-1/2*beta*y[2]^2)/exp(beta*y[1]^2)^2*beta/(1+exp(y[3]-c))^2)^(1/2)*y[2]
    D[2] = -(4*beta*y[1]*exp(y[3]-c)-y[2]*exp(y[3]-c)+4*beta*y[1]+y[2])*pi^(1/2)/(exp(y[3])/exp(c)*exp(-1/2*beta*y[2]^2)/exp(beta*y[1]^2)^2*beta/(1+exp(y[3]-c))^2)^(1/2)/(2+2*exp(y[3]-c))
    D[3] = -1/2*(beta*y[2]^2-2)*pi^(1/2)/(exp(y[3])/exp(c)*exp(-1/2*beta*y[2]^2)/exp(beta*y[1]^2)^2*beta/(1+exp(y[3]-c))^2)^(1/2)

    #Actualización de los coeficientes
    for i in 1:dim
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

function contactintegration{T<:Real}(nsteps::Int64, condinicial::Array{T,1}, timesampling::T, c::T, beta::T=1.)
  t::Float64 = 0.0

  x = condinicial
  #Crea los arreglos donde almacenaré la información de la simulación
  q = Array(Float64, nsteps)
  p = Array(Float64, nsteps)
  S = Array(Float64, nsteps)
  tiempo = Array(Float64, nsteps)
  #El primer paso corresponde a la condición inicial
  q[1] = x[1]
  p[1] = x[2]
  S[1] = x[3]
  #Éste es el arreglo del tiempo de muestreo.
  tiempo = [timesampling*(i-1) for i in 1:nsteps]
  temporarytime = 0.

  for i in 2:nsteps
    j = true
    while j
      #Integro para encontrar la serie de Taylor con su t de convergencia
      t, vec1T::Array{Taylor1{Real},1} = taylorStepper(campoContacto, x, ordenTaylor, tolerance, c, beta )
      #Evalúo las series a tiempo t

      for k in eachindex(x)
        @inbounds x[k] = evaluate( vec1T[k], t )
      end
      #Controlo el tiempo para saber si es momento de muestrear
      temporarytime += t

      if tiempo[i] < temporarytime
        #Si me pasé del tiempo de muestreo,evalúo la última serie a tiempo[i] (pero considerando que son deltas)
        t_last = temporarytime -t
     for k in eachindex(x)
          @inbounds x[k] = evaluate( vec1T[k], tiempo[i] - t_last)
        end

        temporarytime = tiempo[i]
        q[i] = x[1]
        p[i] = x[2]
        S[i] = x[3]
        j = false

      end
    end
  end

  tiempo, q, p, S


end

end
#condinicial
#x = [0.0, 1.0,0.001]

#Hay que checar que una vez corre la simulacion la condicion inicial cambia.

#x = contactIntegration(campoContacto, 100, x) ###Después se grafica
