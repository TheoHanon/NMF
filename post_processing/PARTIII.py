import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

### This script plots :
#       * Spectrum of the wave packet
#       * Group velocity of the wave packet for each scheme



# Constants ----------------------------

p = 12
a = 3/5
L = 1
U = 1  


kp = 2*np.pi*p/L  
sigma = L/16  
N = 128 
x = np.linspace(-1/2, 1/2, N) 
dx = x[1] - x[0]  


# Spectrum of the wave packet ----------------------------


def initial_wave_packet(x, U, kp, sigma):
    return U * np.cos(kp * x) * np.exp(-(x**2) / sigma**2)


wave_packet = initial_wave_packet(x, U, kp, sigma)
wave_fft = np.fft.fft(wave_packet)
wave_fft_magnitude = np.abs(wave_fft)
j = np.fft.fftfreq(N, d=dx) * dx * N  



plt.figure(figsize=(10, 5))

plt.stem(j, wave_fft_magnitude, 'k', markerfmt=".", basefmt="-k", )
plt.title('Spectrum of the Wave Packet')
plt.xlabel('j')
plt.ylabel('Magnitude of Fourier Coefficients')
plt.grid(True)

plt.show()


# Group velocity of the wave packet for each scheme ----------------------------

j = np.arange(0, N/2 +1)
k = 2 * np.pi * j / L 

c_E2 = lambda kh: np.cos(kh)
c_E4  =lambda kh: 4/3 * np.cos(kh) - 1/3 * np.cos(2*kh)
c_I4 = lambda kh: 3/2 * (1/2 + np.cos(kh))/(1 + 1/2 * np.cos(kh))**2
c_ED = lambda kh: (1/3 * np.cos(2*kh) - 4/3 * np.cos(kh) ) + 1j *(2/3 *np.sin(kh) - 1/3 * np.sin(2*kh))


plt.figure()

plt.plot(j, np.abs(c_ED(k*dx)),'k:' ,label='ED Scheme')
plt.plot(j, c_E2(k*dx), 'k-', label='E2 Scheme') 
plt.plot(j, c_E4(k*dx), 'k-.', label='E4 Scheme') 
plt.plot(j, c_I4(k*dx), 'k--', label='I4 Scheme') 
plt.xlabel('j')
plt.ylabel(r'$\dfrac{c_{g}^{*}}{c}$')
plt.grid()
plt.legend()

plt.show()
