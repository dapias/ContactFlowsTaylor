include("./src/steeperandstepsize.jl")

#n será el número de grados de libertad. La dimensión del espacio fase será
# 2n + 1.

# function generarics(n)
#   ics = zeros(2n+1)
#   #Primeras n variables = q^i
#   #Segundas n variables = p_i
#   #Última variable = S
#   for i in 1:length(ics)
#     ics[i] = 1.0 ##Es algo que tenemos que cambiar para mostrar que la ergodicidad no depende de las condiciones iniciales.
#   end
#   ics
# end

#Tengo la libertad de elegit f(S) y p_target(p_i, q^i)

#k es tomada igual a 1 y entonces beta = 1/T

function campoNoseHoover{T<:Real}(x ::Array{T,1}, beta::Float64)
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

    #1,2,3,4 corresponden a q,p, eta, p_eta

    #Ecuaciones de movimiento
    D[1] = y[2]
    D[2] = - y[4]*y[2] - y[1]
    D[3] = y[4]
    D[4] = y[2]^2-1/beta

    #Actualización de los coeficientes
    for i in 1:dim
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

function NoseHooverIntegration(campo, nsteps, condinicial, beta = 1)
  t = 0.0
  x = condinicial

  q = Array(Float64, nsteps)
  p = Array(Float64, nsteps)
  eta = Array(Float64, nsteps)    #Pueden obviarse pero igual quiero ver la evolución esta vez
  p_eta = Array(Float64, nsteps)  #Pueden obviarse pero igual quiero ver la evolución esta vez
  tiempo = Array(Float64, nsteps)

  q[1] = x[1]
  p[1] = x[2]
  eta[1] = x[3]
  p_eta[1] = x[4]
  tiempo[1] = t

  for i in 2:nsteps
    t, x = taylorStepper(campo, x, ordenTaylor, epsAbs, beta )
    q[i] = x[1]
    p[i] = x[2]
    eta[i] = x[3]
    p_eta[i] = x[4]
    tiempo[i] = t + tiempo[i-1]
  end

  tiempo, q, p, eta, p_eta

end

#condinicial
#x = [0.0, 1.0,0.001]

#Hay que checar que una vez corre la simulacion la condicion inicial cambia.

#x = contactIntegration(1, campoContacto, 100, x) ###Después se grafica
