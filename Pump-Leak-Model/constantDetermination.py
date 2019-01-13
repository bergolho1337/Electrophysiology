import math

F = 9.64e+04
Na_e = 145.0
K_e = 3.5
Na_i = 10.0
K_i = 130.0
g_Na = 1.22258163e-02
g_K = 6.28564508e-02
Cm = 4.16121841e-03
Gamma = 6.4103e+03
nw_chap = 1.7e-05
alpha_i0 = 1.0/1.15
Vm = -70.0
Imax = 13.0
Kk = 2.0
Kna = 7.7



# Testando equilibrio da equacao Vm
VM = (-0.001*F*(1-alpha_i0)*(Na_e+K_e))/(Cm*Gamma)
print "{:.8e}".format(VM)

# Testando o equilibrio da equacao do Na_i
Vna = 31.9*math.log(Na_e/Na_i)
BOMBA = (1+(Kk/K_e))**(-2) * (1+(Kna/Na_i))**(-3)
NAI = g_Na*(Vm-Vna) + 2*Imax*BOMBA
print "{:.8e}".format(NAI)

# Testando o equilibrio da equacao do K_i
Vk = 31.9*math.log(K_e/K_i)
KI = g_K*(Vm-Vk) - 3*Imax*BOMBA
print "{:.8e}".format(KI)

print "{:.8e}".format(alpha_i0)


