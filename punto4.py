print("---------------------------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------------------------")
print("---------------------------------------------------- Ejercicio 4 ----------------------------------------------")
print("---------------------------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------------------------")

# Ejercicio4
# Tiene una serie de tiempo, donde los datos de tiempo estan almacenados en un arreglo t y los datos de amplitud en un arreglo amp.
#1) Usando los paquetes de scipy de la transformada de Fourier, haga un FILTRO de la senial que elimine las frecuencias mayores a 1000Hz en las senial original.
#2) Haga una grafica de la senial original y la senial filtarada y guardela SIN MOSTRARLA en filtro.pdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

n = 1280 # number of point in the whole interval
f = 200.0 #  frequency in Hz
dt = 1 / (f * 320 ) #320 samples per unit frequency
t = np.linspace( 0, (n-1)*dt, n)
amp = np.cos(2 * np.pi * f * t) - 0.4 * np.sin(2 * np.pi * (2*f) * t ) + np.random.random(n)

# SU FILTRO

def fun(n):
    N = len(n)
    n1= np.arange(len(n))
    F_k=[]
    for k in range(N):
        omega = 2*np.pi*n1*k/N
        F_real = x*np.cos(2*np.pi*n1*k/N)
        F_im = x*np.sin(2*np.pi*n1*k/N)
        Freal = np.sum(F_real)
        Fimag = np.sum(F_im)
        Fk= Freal -1j*Fimag
        F_k.append(Fk)
return F_k

f_sig=fftfreq(n, dt)

filtro=1000.0 
fourier=fun(f)

filtrado=[]
for i in range(len(f_sig)):
    if abs(f_sig[i]) > filtro:filtrado.append(0)
        else:
            filtrado.append(fourier[i]) 
filtrada=ifft(np.array(filtrado))

# SU GRAFICA


# Puede usar los siguientes paquetes:
#from scipy.fftpack import fft, fftfreq, ifft


plt.figure()
plt.plot(t, amp,label="Signal")
plt.plot(x_sig, filtrada.real)
plt.xlabel("t")
plt.ylabel("amp")
plt.title("Signal")
plt.savefig("filtro.pdf")
