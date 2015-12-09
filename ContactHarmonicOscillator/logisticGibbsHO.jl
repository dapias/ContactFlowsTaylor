include("./steeper.jl")

function campoContacto{T<:Real}(x ::Array{T,1},c::Float64, beta::Float64)
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

function contactIntegration{T<:Number}(campo::Function, nsteps::Int64, condinicial::Array{T,1}, timestep::Float64, c::Float64, beta=1.)
  t = 0.0

  x = condinicial

  n = length( x )

  q = Array(Float64, nsteps)
  p = Array(Float64, nsteps)
  S = Array(Float64, nsteps)
  tiempo = Array(Float64, nsteps)

  q[1] = x[1]
  p[1] = x[2]
  S[1] = x[3]
  tiempo = [timestep*(i-1) for i in 1:nsteps]

  temporarytime = 0.

  for i in 2:nsteps

    j = true
    while j

      t, vec1T = taylorStepper(campo, x, ordenTaylor, epsAbs, c, beta )

      for k=1:n
          x[k] = evaluate( vec1T[k], t )
      end

      temporarytime += t
      #println("time = $temporarytime")

      if tiempo[i] < temporarytime

      for k=1:n
          x[k] = evaluate( vec1T[k], tiempo[i] - (temporarytime - t))
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

#condinicial
#x = [0.0, 1.0,0.001]

#Hay que checar que una vez corre la simulacion la condicion inicial cambia.

#x = contactIntegration(campoContacto, 100, x) ###Después se grafica
