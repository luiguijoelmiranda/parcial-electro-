import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin, pi
from scipy import constants
import scipy.constants as sc
# declaramos variables y costantes
ke = 100
ex = np.zeros(ke)
hy = np.zeros(ke)
E0 = 1
eps = np.zeros(ke)
#tamaño de las laminas 
l1_start = 33
l1_end = 40
l2_start = 53
l2_end = 60
# Tamaño del dieléctrico
diel_start = 40
diel_end = 53
# permitividad electrica en los diferentes medios
eps[0:l1_start] = sc.epsilon_0
eps[l1_start:l1_end] = math.sqrt(12)*sc.epsilon_0
eps[diel_start:diel_end] = 12*sc.epsilon_0
eps[l2_start:l2_end] = math.sqrt(12)*sc.epsilon_0
eps[l2_end:ke] = sc.epsilon_0
mu = np.zeros(ke)
# permeabilidad magnetica en los diferentes medios
mu[0:ke] = sc.mu_0
# condiciones para el seno
vel = sc.c
dt = 0.001e-9
dx = 2*dt*vel
freq_in = 10e9
# Condiciones en los bordes
boundary_low = [0, 0]
boundary_high = [0, 0]
# Impedancias intrínsecasde cada medio
n0 = math.sqrt(sc.mu_0/sc.epsilon_0)
n1 = math.sqrt(mu[diel_start+1]/eps[diel_start+1])
n2 = math.sqrt(sc.mu_0/sc.epsilon_0)
# Coeficiente de reflexión de cada interfaz
r01 = (n1-n0)/(n1+n0)
r12 = (n2-n1)/(n2+n1)
# Coficiente de transmisión de cada interfaz
t01 = (2*n1)/(n1+n0) 
t12 = (2*n2)/(n2+n1)  
nsteps = 275
for time_step in range(1, nsteps + 1):        
    for k in range(1, ke):
        ex[k] = ex[k] + (dt/(eps[k]*dx))*(hy[k-1]-hy[k])
    if time_step < 130 : 
        pulse = sin(2 * pi * freq_in * dt * time_step)
    ex[0] = pulse
    if time_step > 400:
        ex[1] = pulse
        ex[0] = boundary_low.pop(0)
        boundary_low.append(ex[1])
        ex[ke - 1] = boundary_high.pop(0)
        boundary_high.append(ex[ke - 2])
    for k in range(ke - 1):
        hy[k] = hy[k] + (dt/(mu[k]*dx))*(ex[k] - ex[k+1])
    x0 = (np.linspace(-100,100, num = ke))*dx
    plt.plot(x0,ex[:])# graficamos
    plt.ylabel("Campo Electrico (V/m)")
    plt.xlabel("Posicion (m)")
    plt.show()
print (ex[60:100].max())# se aumento en un 20% el desempeño de la antena