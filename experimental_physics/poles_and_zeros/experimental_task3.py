import numpy as np
import matplotlib.pyplot as plt

#from theoretical_task3 import *

node = 'V4' # V2=T_H, V3=T_B, V4=T_L
resistor = '33' # 10, 33, 150. In k Ohm
interval = '10-20k'
filename = str(node + '_R2-' + resistor + 'k_' + interval + '.txt')
with open(str('C:/Users/sebas/.vscode/experimental_physics/poles_and_zeros/data/' + filename)) as f:
    f.readline()
    f.readline()
    f.readline()
    lines = f.readlines()

freq = []
gain = []
phase = []
for x in lines:
    freq.append(float(x.split()[0]))
    gain.append(float(x.split()[1]))
    phase.append(float(x.split()[2]))

amplitude = [np.power(10,g/20) for g in gain]
x = [r * np.cos(np.pi / 180 * f) for r,f in zip(amplitude, phase)]
y = [r * np.sin(np.pi / 180 * f) for r,f in zip(amplitude, phase)]

match node:
    case 'V2':
        trans = '$T_H$'
    case 'V3':
        trans = '$T_B$'
    case 'V4':
        trans = '$T_L$'
    case _:
        print('Incorrect node name')

plt.figure(0)
plt.semilogx(freq,gain)
plt.title(str('Bode amplitude of ' + trans + ', $R_2$ =' + resistor + ' k$\Omega$, experimental'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.grid()

plt.figure(1)
plt.semilogx(freq,phase)
plt.title(str('Bode phase of ' + trans + ', $R_2$ =' + resistor + ' k$\Omega$, experimental'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shifft (deg)')
plt.grid()

plt.figure(2)
plt.plot(x,y)
plt.plot(x[0], y[0], 'rx', label=str(interval[0] + interval[1] + ' $\Omega$'))
plt.plot(x[-1], y[-1], 'go', label=str(interval[3] + interval[4] + ' k$\Omega$'))
plt.title(str('Nyquist plot of ' + trans + ', $R_2$ =' + resistor + ' k$\Omega$, experimental'))
plt.legend()
plt.grid()

plt.show()