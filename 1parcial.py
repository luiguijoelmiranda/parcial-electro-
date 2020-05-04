import numpy as np
from math import sin, pi, exp
import scipy.constants as sc
from matplotlib import pyplot as plt
#variables y costantes
ke = 100
kc = 50
t0 = 40
ddx = 0.001 # discretización de las celadas
dt = ddx / 12.5e8 # tiempo de paso
freq_in = 10e9 # frecuencia minima
spread = 12
# permitividad electrica en los diferentes medios
cb = np.ones(ke)
cb = 0.5 * cb
#bordes absorventes
boundary_low = [0, 0]
boundary_high = [0, 0]
cb_start = 100
epsilon = 12 #epsolon del nuevo medio
cb[44:55] = 0.5 / epsilon #cambio de medio 
nsteps = 440 # paso del tiempo
# Impedancias intrínsecasde cada medio
n1 = (sc.mu_0/sc.epsilon_0)**(1/2)
n2 = (sc.mu_0/(12*sc.epsilon_0))**(1/2)
# Coeficiente de reflexión y transmisión de cada interfaz 
t12 =(2*n1)/(n1+n2) 
r12 =(n2-n1)/(n1+n2)
t21 =(2*n2)/(n1+n2)
r21 =(n1-n2)/(n1+n2)
s = 1/(1-(r21*r21))
#campo electrico y magnetico
ex = np.zeros((ke), dtype = float)
hy = np.zeros((ke), dtype = float)
for time_step in range(1, nsteps + 1):
# calculando los campos
    for k in range(1, ke):
        ex[k] = ex[k] + cb[k] * (hy[k - 1] - hy[k])
        # pulo seno
    if time_step < 130 :
        pulse = sin(2 * pi * freq_in * dt * time_step)
    ex[1] = pulse
    ex[0] = boundary_low.pop(0)
    boundary_low.append(ex[1])
    ex[ke - 1] = boundary_high.pop(0)
    boundary_high.append(ex[ke - 2])
    for k in range(ke - 1):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])
x0 = (np.linspace(-100,100, num = ke))*ddx
plt.plot(x0,ex[:])# graficamos
plt.ylabel("Campo Electrico (V/m)")
plt.xlabel("Posicion (m)")
plt.show()
print(time_step)
Max = ex[56:100].max()
Maxt = t12*t21
print (Max, Maxt)    
    