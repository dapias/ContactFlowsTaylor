function rungeK(x_init,campo,h)
  state = x_init
  k1 = h*campo(state)
  k2 = h*campo(state + k1)
  k3 = h*campo(state + (3*k1 + k2)/8)
  k4 = h*campo(state +(8*k1 + 2*k2 + 8*k3)/27)
  k5 = h*campo(state + (3*(3*sqrt(21)- 7)*k1 - 8*(7-sqrt(21))*k2
                        + 48*(7-sqrt(21))*k3 - 3*(21 - sqrt(21))*k4
                        )/392)
  k6 = h*campo(state + (-5*(231 + 51*sqrt(21))*k1 - 40*(7 + sqrt(21))*k2
                        -320*sqrt(21)*k3 + 3*(21 + 121*sqrt(21))*k4
                        + 392*(6 + sqrt(21))*k5)/1960)
  k7 = h*campo(state + (15*(22+ 7*sqrt(21))*k1 + 120*k2
               + 40*(7*sqrt(21)-5)*k3 - 63*(3*sqrt(21)-2)*k4
                        -14*(49+9*sqrt(21))*k5 + 70*(7-sqrt(21))*k6)/180)

  state = state+(9*k1 +64*k3 + 49*k5 + 49*k6 + 9*k7)/180.
  return state
end
