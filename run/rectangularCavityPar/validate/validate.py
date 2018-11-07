#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

convertfunc = lambda x: float(x.decode('utf-8').strip('()'))

data = np.genfromtxt('../probes1/0/U', converters={1: convertfunc, 3: convertfunc})

t = data[:, 0]
v = data[:, 2]

fig, ax = plt.subplots()
ax.plot(t, v)

dt = t[1]-t[0]
nperseg = v.shape[0]//6 # six cycles are present in the data
noverlap = nperseg//2
f, S = signal.welch(v, fs=1./dt,
        nperseg=nperseg,
        noverlap=noverlap)

# transform to radians
w = f*2*np.pi
Sw = S/2/np.pi

Sw_0 = 0.008 * np.ones(w.shape)

fig1, ax1 = plt.subplots()
#ax1.semilogy(w, Sw, w, Sw_0, 'r--')
ax1.plot(w, Sw, w, Sw_0, 'r--')
ax1.set_xlim([0, 10*2*np.pi])

plt.show()
