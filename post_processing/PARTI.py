import numpy as np
from numpy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt



### Plot of Question: 
#       * Discretized periodic domain 
#       * Partially decentered scheme 
#       * Stability   

                                                       


## Discretized periodic domain ----------------------------
L = 1  
U = 1
Sigma = np.array([L/4, L/16])
H = Sigma / 8

for i, (h,sigma) in enumerate(zip(H, Sigma)):

    fig, ax = plt.subplots(1, 1)
    N = int(L/h)
    j = np.arange(0, N/2 +1)
    k = 2 * np.pi * j / L 

    X = np.arange(-L/2, L/2, h)  
    u_gauss = U * np.exp(-X**2 / (sigma**2)) 

    fft_gauss = np.abs(rfft(u_gauss)) 
    fft_freqs = rfftfreq(N, d=h) 

    Q = U * np.sqrt(np.pi) * sigma
    unbounded_fft_gauss = np.abs(Q*np.exp(-sigma**2 * k**2 / 4))

    ax.semilogy(fft_freqs*h*N, fft_gauss, 'o', label='DFT', color = 'k')
    ax.semilogy(j, unbounded_fft_gauss, label=f'Continuous FT', color = 'k')
    ax.xticks = j[::2]
    ax.set_ylabel('Magnitude of the Fourier Coefficients')
    ax.grid(True)
    ax.legend(loc = 'lower left')
    ax.set_xlabel('j')
    plt.show()


## Partially decentered scheme ----------------------------

kh_star = lambda kh : 1/6 * np.sin(2*kh) - 4/3 * np.sin(kh) + 1j *(1/2 + 1/6 * np.cos(2*kh) - 2/3 * np.cos(kh))
kh = np.linspace(0.0001, np.pi, 1000)


fig, ax = plt.subplots(1,1, figsize=(8, 6))

ax.plot(kh/np.pi, np.abs(kh_star(kh)) / np.pi, label='$|k^*h|$', color = 'k')
ax.plot(kh/np.pi, np.abs(kh) / np.pi, label='reference $|kh|$', linestyle='--', color = "k", linewidth = 2)
ax.set_xlabel('$\dfrac{kh}{\pi}$', fontsize = 12)
ax.set_ylabel('$\dfrac{|k^*h|}{\pi}$', fontsize = 12)
ax.legend()
ax.grid(True)

plt.show()

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

ax.plot(kh/np.pi, np.angle(kh_star(kh)) / np.pi, label='arg($k^*h$)', color = 'k')
ax.plot(kh/np.pi, (np.angle(kh) + np.pi) / np.pi, label='reference arg($kh$)', linestyle='--', color = "k", linewidth = 2)
ax.set_xlabel('$\dfrac{kh}{\pi}$', fontsize = 12)
ax.set_ylabel('$\dfrac{arg(k^*h)}{\pi}$', fontsize = 12)
ax.legend()
ax.axis('equal')
ax.grid(True)

plt.show()


## Stability ----------------------------------------------

# Rounge-Kutta 4th order
R1 = lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24

X, Y = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100))
Z = X + 1j*Y
stability_outputs1 = R1(Z)

beta  = 14/15

# Decentered scheme stability region

kh = np.linspace(-np.pi, np.pi, 1000)
lam_ED = lambda kh: -(1/2 + 1/6 * np.exp(-1j*2*kh) - np.exp(-1j*kh) + 1/3 * np.exp(1j*kh))
lam_ED_outputs = beta*lam_ED(kh)

# E2 

kh = np.linspace(-np.pi, np.pi, 10)
lam_E2  = lambda kh: -1j*np.sin(kh)
lam_E2_outputs = beta*lam_E2(kh)

# E4  

lam_E4 = lambda kh : -1j * (4/3 * np.sin(kh) - 1/6 * np.sin(2*kh))
lam_E4_outputs = beta*lam_E4(kh)

# I4 

lam_I4 = lambda kh : -1j * (3/2 * np.sin(kh) / (1 + 1/2 * np.cos(kh)))
lam_I4_outputs = beta*lam_I4(kh)

# Contour plots

plt.figure(figsize=(8, 6))
plt.plot(lam_ED_outputs.real, lam_ED_outputs.imag, color='k', linestyle = 'dashdot', linewidth = 2, label = 'ED Scheme')
plt.scatter(lam_E2_outputs.real, lam_E2_outputs.imag, color='b', marker = "x", label = 'E2 Scheme ')
plt.scatter(lam_E4_outputs.real, lam_E4_outputs.imag, color='r', marker = "x", label = 'E4 Scheme ')
plt.scatter(lam_I4_outputs.real, lam_I4_outputs.imag, color='g', marker = "x", label = 'I4 Scheme ')
plt.contour(X, Y, np.abs(stability_outputs1), levels=[1], colors='r')
plt.xlabel('$\Re(\lambda \Delta t$)', fontsize = 12)
plt.ylabel('$\Im(\lambda \Delta t$)', fontsize = 12)
plt.title(r'Stability region $CFL = \dfrac{14}{15}$')
plt.grid(True)
plt.legend()
plt.axis('equal') 
plt.show()












