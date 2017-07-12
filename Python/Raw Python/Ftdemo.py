# ftdemo - Discrete Fourier transform demonstration program

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Initialize the sine wave time series to be transformed.
N = input('Enter the number of points: ')
freq = input('Enter frequency of the sine wave: ')
phase = input('Enter phase of the sine wave: ')
tau = 1.   # Time increment
t = np.arange(N)*tau               # t = [0, tau, 2*tau, ... ]
y = np.empty(N)
for i in range(N):              # Sine wave time series
    y[i] = np.sin(2*np.pi*t[i]*freq + phase)    
f = np.arange(N)/(N*tau)           # f = [0, 1/(N*tau), ... ] 


#* Compute the transform using desired method: direct summation
#  or fast Fourier transform (FFT) algorithm.
yt = np.zeros(N,dtype=complex)
Method = input('Compute transform by: 1) Direct summation; 2) FFT  :')

import time
startTime = time.time()
if Method == 1 :             # Direct summation
    twoPiN = -2. * np.pi * (1j) /N    # (1j) = sqrt(-1)
    for k in range(N):
        for j in range(N):
            expTerm = np.exp( twoPiN*j*k )
            yt[k] += y[j] * expTerm
else:                        # Fast Fourier transform
    yt = np.fft.fft(y)

stopTime = time.time()

print 'Elapsed time = ', stopTime - startTime, ' seconds'


#* Graph the time series and its transform.
plt.subplot(1, 2, 1) # Left plot
plt.plot(t,y)
plt.title('Original time series')
plt.xlabel('Time')

plt.subplot(1, 2, 2) # Right plot
plt.plot(f,np.real(yt),'-',f,np.imag(yt),'--')
plt.legend(['Real','Imaginary  '])
plt.title('Fourier transform')
plt.xlabel('Frequency')

plt.show()

#* Compute and graph the power spectrum of the time series
powspec = np.empty(N)
for i in range(N):
    powspec[i] = abs(yt[i])**2
plt.semilogy(f,powspec,'-')
plt.title('Power spectrum (unnormalized)')
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.show()

