#Tomado  de http://www.ams.org/journals/mcom/1968-22-102/S0025-5718-68-99876-1/S0025-5718-68-99876-1.pdf
#(An explicit sixth-order Runge-Kutta formula) by H.A . Luther


function rungeK(x_init,campo,deltat,beta)
  state = x_init
  k1 = deltat*campo(state,beta)
  k2 = deltat*campo(state + k1, beta)
  k3 = deltat*campo(state + (3*k1 + k2)/8, beta)
  k4 = deltat*campo(state +(8*k1 + 2*k2 + 8*k3)/27, beta)
  k5 = deltat*campo(state + (3*(3*sqrt(21)- 7)*k1 - 8*(7-sqrt(21))*k2
                        + 48*(7-sqrt(21))*k3 - 3*(21 - sqrt(21))*k4
                        )/392, beta)
  k6 = deltat*campo(state + (-5*(231 + 51*sqrt(21))*k1 - 40*(7 + sqrt(21))*k2
                        -320*sqrt(21)*k3 + 3*(21 + 121*sqrt(21))*k4
                        + 392*(6 + sqrt(21))*k5)/1960,beta)
  k7 = deltat*campo(state + (15*(22+ 7*sqrt(21))*k1 + 120*k2
               + 40*(7*sqrt(21)-5)*k3 - 63*(3*sqrt(21)-2)*k4
                        -14*(49+9*sqrt(21))*k5 + 70*(7-sqrt(21))*k6)/180, beta)

  state = state+(9*k1 +64*k3 + 49*k5 + 49*k6 + 9*k7)/180.
  return state
end
