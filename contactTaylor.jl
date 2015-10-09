include("./src/steeperandstepsize.jl")

function campoContacto{T<:Real}(x ::Array{T,1}, beta::Float64)
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

    #Las entradas *D* y *y* estan en el orden q, p, S

    #Ecuaciones de movimiento
    D[1] = (y[end]^2+1)*beta*y[2]*exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))/(2*y[end])
    D[2] = -exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))*(beta*y[1]*y[end]^3-2*y[2]*y[3]^2+beta*y[1]*y[3]+2*y[2])/(2*y[end]^2)
    D[end] = -(y[end]^2+1)*exp(beta*((1/4)*y[2]^2+(1/4)*y[1]^2))*(beta*y[2]^2-2)/(2*y[end])


    #Actualización de los coeficientes
    for i in 1:dim
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

function contactIntegration(campo, nsteps, condinicial, beta=1.)
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
          t, x = taylorStepper(campo, x, ordenTaylor, epsAbs, beta )
    q[i] = x[1]
    p[i] = x[2]
    S[i] = x[3]
    tiempo[i] = t + tiempo[i-1]
  end

  tiempo, q, p, S

end

#condinicial
#x = [0.0, 1.0,0.001]

#Hay que checar que una vez corre la simulacion la condicion inicial cambia.

#x = contactIntegration(campoContacto, 100, x) ###Después se grafica
