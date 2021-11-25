#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from scipy import fft

# Transformada de Fourier
senal_f = fft(senal_Tx)

# Muestras de la señal
Nm = len(senal_Tx)

# Número de símbolos (198 x 89 x 8 x 3)
Ns = Nm // mpp

# Tiempo del símbolo = periodo de la onda portadora
Tc = 1 / fc

# Tiempo entre muestras (período de muestreo)
Tm = Tc / mpp

# Tiempo de la simulación
T = Ns * Tc

# Espacio de frecuencias
f = np.linspace(0.0, 1.0/(2.0*Tm), Nm//2)

# Gráfica
plt.plot(f, 2.0/Nm * np.power(np.abs(senal_f[0:Nm//2]), 2))
plt.xlim(0, 20000)
plt.grid()
plt.show()


# In[ ]:


# --------INICIO DE LA SOLUCIÓN DEL PROYECTO.------


# In[ ]:


# ------Asignación 4.1----------------
# primero hay que empezar modificando la función de modulacion a 8-psk
import numpy as np

def modulador1(bits, fc, mpp):
    '''Un método que simula el esquema de 
    modulación digital 8-PSK.

    :param bits: Vector unidimensional de bits
    :param fc: Frecuencia de la portadora en Hz
    :param mpp: Cantidad de muestras por periodo de onda portadora
    :return: Un vector con la señal modulada
    :return: Un valor con la potencia promedio [W]
    :return: La onda portadora c(t)
    :return: La onda cuadrada moduladora (información)
    '''
# 1. Parámetros de la 'señal' de información (bits)
    N = len(bits) # Cantidad de bits

# 2. Construyendo un periodo de la señal portadora c(t)
    Tc = 1 / fc  # periodo [s]
    t_periodo = np.linspace(0, Tc, mpp)  # mpp: muestras por período
    portadora1 = np.cos(2*np.pi*fc*t_periodo)
# se debe agregar una segunda portadora
    portadora2 = np.sin(2*np.pi*fc*t_periodo)

# 3. Inicializar la señal modulada s(t)
    t_simulacion = np.linspace(0, N*Tc, N*mpp) 
    senal_Tx = np.zeros(t_simulacion.shape)
    moduladora = np.zeros(t_simulacion.shape)  
    
# Se define el valor de h 
    h = 0.707
  
            #A1 = 1, & A2 = 0 &   b1b2b3 = 111 
            # A1 = h, & A2 = h &   b1b2b3 = 110 
            #A1 = 0, & A2 = 1 &   b1b2b3 = 010 
            #A1 =-h, & A2 = h &   b1b2b3 = 011 
            #A1 =-1, & A2 = 0 &   b1b2b3 = 001 
            #A1 =-h, & A2 =-h &   b1b2b3 = 000 
            #A1 = 0, & A2 =-1 &   b1b2b3 = 100
            #A1 = h, & A2 =-h &   b1b2b3 = 101 
# hay que tener el formato anterior para el modulador
# por ende el range se vera ido de 3 en 3 par asi poder acceder a los 3 bits que necesitamos
    for i in range (0, N, 3):
        
        #111 
        if bits[i] == 1 and bits[i+1] == 1 and bits[i+2] == 1:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * 1 + portadora2 * 0
            #la señal vendra de esa forma ya que recoremos que son un sen y un cos por ende la señal vendra con esa
            #forma como en las siguientes bits
        
        
        # 000
        elif bits[i] == 0 and bits[i+1] == 0 and bits[i+2] == 0:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * -h + portadora2 * -h
            
        
        
        # 101
        elif bits[i] == 1 and bits[i+1] == 0 and bits[i+2] == 1:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * h + portadora2 * -h
        
        
         #010
        elif bits[i] == 0 and bits[i+1] == 1 and bits[i+2] == 0:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * 0 + portadora2 * 1
            
        
        #110 
        elif bits[i] == 1 and bits[i+1] == 1 and bits[i+2] == 0:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * h + portadora2 * h
            
        
         # 100
        elif bits[i] == 1 and bits[i+1] == 0 and bits[i+2] == 0:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * 0 + portadora2 * -1
        
       
            
        
        # 011
        elif bits[i] == 0 and bits[i+1] == 1 and bits[i+2] == 1:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * -h + portadora2 * h
            
        # 001
        elif bits[i] == 0 and bits[i+1] == 0 and bits[i+2] == 1:
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1 * -1 + portadora2 * 0
            
            
            
       
            
            
        
            
        
      
           
# 5. Calcular la potencia promedio de la señal modulada
    P_senal_Tx = (1 / (N*Tc)) * np.trapz(pow(senal_Tx, 2), t_simulacion)
    
    return senal_Tx, P_senal_Tx, portadora1, portadora2, moduladora  


# In[ ]:


# segundo hay que modificar la función del demoludador para asi demule en 8-psk
import numpy as np

def demodulador(senal_Rx, portadora, portadora2, mpp):
    '''Un método que simula un bloque demodulador
    de señales, bajo un esquema 8PSK. El criterio
    de demodulación se basa en decodificación por 
    detección de energía'''

    :param senal_Rx: La señal recibida del canal
    :param portadora: La onda portadora c(t)
    :param mpp: Número de muestras por periodo
    :return: Los bits de la señal demodulada
    '''
# Cantidad de muestras en senal_Rx
    M = len(senal_Rx)

# Cantidad de bits (símbolos) en transmisión
    N = int(M / mpp)

# Vector para bits obtenidos por la demodulación
    bits_Rx = np.zeros(N)

# Vector para la señal demodulada
    senal_demodulada = np.zeros(senal_Rx.shape)

# energía  de la primera portadora
    Es = np.sum(portadora1 * portadora1)
    
# energia para la portadora 2
    E2 = np.sum(portadora2 * portadora2)

# Demodulación
    for i in range(N):
        # Producto de la portadora 1
        producto1 = senal_Rx[i*mpp : (i+1)*mpp] * portadora1
        Ep1 = np.sum(producto1) 
        
        #para la portadora 2
        producto2 = senal_Rx[i*mpp : (i+1)*mpp] * portadora2
        Ep2 = np.sum(producto2)
        
# sumamos las dos funciones para asi obtener la total
        
        senal_demodulada[i*mpp : (i+1)*mpp] = producto1 + producto2

        ''' 
           Los  umbrales de energia son los siguientes
            -(1+h)/2 > E ---- -1
            -h/2  > E > -(1+h)/2 ---- -h
            h/2 > E > -h/2 ----- 0
            (1+h)/2 > E > h/2 ---- h 
            E> (1+h)/2 -----1
            
        '''
        h = 0.707
        tuc = 1 + h
 
 # se empezaran a poner los umbrales de las 2 ondas
# tanto para la portadora 1 como la 2
# para que vaya segun el ciclo
# y vaya pon8endo el umbral entre cada bit
        if (h/2)*E2 >= Ep2 and Ep2 >= (-h/2)*E2 and Ep1 >= (tuc/2)*Es:
            bits_Rx[i]   = 1
            bits_Rx[i+1] = 1
            bits_Rx[i+2] = 1
        
            
        
        elif Ep1 >= (h/2)*Es and Ep2 <= (tuc/2)*E2 and Ep2 >= (h/2)*E2 and Ep1 <= (tuc/2)*Es:
            bits_Rx[i]   = 1
            bits_Rx[i+1] = 1
            bits_Rx[i+2] = 0
          
       
        
        elif (h/2)*Es >= Ep1 and  Ep1 >= (-h/2)*Es and Ep2>= (tuc/2)*E2 :
            bits_Rx[i]   = 0
            bits_Rx[i+1] = 1
            bits_Rx[i+2] = 0
                
       
       
                
        elif (h/2)*E2 >= Ep2 and Ep2 >= (-h/2)*E2 and (-tuc/2)*Es >= Ep1  :
            bits_Rx[i]   = 0
            bits_Rx[i+1] = 0
            bits_Rx[i+2] = 1
                
        
        elif (-h/2)*Es  >= Ep1 and Ep1 >= (-tuc/2)*Es and Ep2 <= (tuc/2)*E2 and Ep2 >= (h/2)*E2:
            bits_Rx[i]   = 0
            bits_Rx[i+1] = 1
            bits_Rx[i+2] = 1
        
        elif (-h/2)*Es >= Ep1 and Ep1 >= (-tuc/2)*Es and (-h/2)*E2 >= Ep2 and Ep2>= (-tuc/2)*E2:
            bits_Rx[i]   = 0
            bits_Rx[i+1] = 0
            bits_Rx[i+2] = 0
                
    
                
        elif (tuc/2)*Es >= Ep1 and Ep1 >= (h/2)*Es and (-h/2)*E2 >= Ep2 and Ep2 >= -tuc/2*E2:
            bits_Rx[i]   = 1
            bits_Rx[i+1] = 0
            bits_Rx[i+2] = 1
        
        elif  (h/2)*Es >= Ep1 and Ep1 >= (-h/2)*Es and (-tuc/2)*E2 >= Ep2 :
            bits_Rx[i]   = 1
            bits_Rx[i+1] = 0
            bits_Rx[i+2] = 0

    return bits_Rx.astype(int), senal_demodulada


# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import time

# Parámetros
fc = 5000  # frecuencia de la portadora
mpp = 20   # muestras por periodo de la portadora
SNR = 5  # relación señal-a-ruido del canal

# Iniciar medición del tiempo de simulación
inicio = time.time()

# 1. Importar y convertir la imagen a trasmitir
imagen_Tx = fuente_info('arenal.jpg')
dimensiones = imagen_Tx.shape

# 2. Codificar los pixeles de la imagen
bits_Tx = rgb_a_bit(imagen_Tx)

# 3. Modular la cadena de bits usando el esquema BPSK
senal_Tx, P_senal_Tx, portadora1, portadora2, moduladora = modulador1(bits_Tx, fc, mpp)

# 4. Se transmite la señal modulada, por un canal ruidoso
senal_Rx = canal_ruidoso(senal_Tx, P_senal_Tx, SNR)

# 5. Se desmodula la señal recibida del canal
bits_Rx, senal_demodulada = demodulador(senal_Rx, portadora1, portadora2, mpp)

# 6. Se visualiza la imagen recibida 
imagen_Rx = bits_a_rgb(bits_Rx, dimensiones)
Fig = plt.figure(figsize=(10,6))

# Cálculo del tiempo de simulación
print('Duración de la simulación: ', time.time() - inicio)

# 7. Calcular número de errores
errores = sum(abs(bits_Tx - bits_Rx))
BER = errores/len(bits_Tx)
print('{} errores, para un BER de {:0.4f}.'.format(errores, BER))

# Mostrar imagen transmitida
ax = Fig.add_subplot(1, 2, 1)
imgplot = plt.imshow(imagen_Tx)
ax.set_title('Transmitido')

# Mostrar imagen recuperada
ax = Fig.add_subplot(1, 2, 2)
imgplot = plt.imshow(imagen_Rx)
ax.set_title('Recuperado')
Fig.tight_layout()

plt.imshow(imagen_Rx)

#mostramos las ondas
import matplotlib.pyplot as plt

# Visualizar el cambio entre las señales
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, sharex=True, figsize=(14, 7))

# La onda cuadrada moduladora (bits de entrada)
ax1.plot(moduladora[0:600], color='r', lw=2) 
ax1.set_ylabel('$b(t)$')

# La señal modulada por BPSK
ax2.plot(senal_Tx[0:600], color='g', lw=2) 
ax2.set_ylabel('$s(t)$')

# La señal modulada al dejar el canal
ax3.plot(senal_Rx[0:600], color='b', lw=2) 
ax3.set_ylabel('$s(t) + n(t)$')

# La señal demodulada
ax4.plot(senal_demodulada[0:600], color='m', lw=2) 
ax4.set_ylabel('$b^{\prime}(t)$')
ax4.set_xlabel('$t$ / milisegundos')
fig.tight_layout()
plt.show()


# In[ ]:


#asignación 2
#Realice pruebas de estacionaridad y ergodicidad a la señal modulada senal_Tx y obtenga conclusiones sobre estas.
import numpy as np
import matplotlib.pyplot as plt

# Tiempo para la muestra
T= 100 
t_final = 10
t=np.linspace(0, t_final, T)
fc = 5000

#Posibles combinaciones
X=8

# Valor A_1
A_1=[1, 0.707, 0, -0.707, -1, -0.707, 0, 0.707]

# Valor A_2
A_2=[0, 0.707, 1, 0.707, 0, -0.707, -1, -0.707]

# Nueva figura 
plt.figure()

combinaciones= np.zeros((X, len(t)))

#Forma de onda
for i in range (len(A_2)):
    s_t = A_1[i]*(np.cos(2*np.pi*fc*t)) + A_2[i]*(np.sin(2*np.pi*fc*t))
    combinaciones[i, :] = s_t
    plt.plot(t, s_t)
                  
#X posibilidades 
media_X = [np.mean(combinaciones[:, i]) for i in range (len(t))]

#Promedio de de combinaciones
valor_real =plt.plot(t, media_X, lw=7, color='b', label = 'Valor promedio de la srealizaciones')
                  
#Valo teorico
teorico = np.mean(senal_Tx)*t
                  
valor_teorico =plt.plot(t, teorico, '--',lw= 4, color='r', label = 'Valor teórico')
                  
# 8. Mostrar las realizaciones, y su promedio calculado y teórico
plt.title('Realizaciones del proceso aleatorio $X(t)$')
plt.xlabel('$t$')
plt.ylabel('$x_i(t)$')
plt.legend()
plt.show()

'''Se sabe que un sistema es ergodico si la hipersuperficie de energía 
constante del espacio de las fases es toda la hipersuperficie de energía  y esta energía es constante, como se observa
en la gráfica la energía se mantiene constante por lo que el sistema es ergodico. Además si el promedio se mantiene t
es estacionario, de manera que su varianza es nula, cosa que también se observa en la gráfica.'''


# In[ ]:


# 4.3. - Densidad espectral de potencia
# Determine y grafique la densidad espectral de potencia para la señal modulada 'senal_Tx'.

# Importamos la librerías y módulos requeridos

import numpy as np
from scipy import fft
import matplotlib.pyplot as plt

# Transformada de Fourier de la señal t_x
senal_f = fft.fft(senal_Tx)

# Muestras de la señal
Nm = len(senal_Tx)

# Número de símbolos (198 x 89 x 8 x 3)
Ns = Nm // mpp

# Tiempo del símbolo = periodo de la onda portadora
Tc = 1 / fc

# Tiempo entre muestras (período de muestreo)
Tm = Tc / mpp

# Tiempo de la simulación
T = Ns * Tc

# Espacio de frecuencias
f = np.linspace(0.0, 1.0/(2.0*Tm), Nm//2)

# Densidad espectral de potencia Sxx = |s(w)|^2
S_xx = np.power(np.abs(senal_f), 2)

# Gráfica
plt.plot(f, 2.0/Nm * np.power(np.abs(senal_f[0:Nm//2]), 2))
plt.xlim(0, 20000)
plt.title('Densidad espectral de potencia para señal t_x')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Densidad espectral de potencia')
plt.legend() # leyenda de la gráfica.
plt.grid()   
plt.show()   # Muestra la gráfica.

