import numpy as np
import matplotlib.pyplot as plt

import theoretical_task3 

node = 'V3' # V2=T_H, V3=T_B, V4=T_L
resistor = '150' # 10, 33, 150. In k Ohm
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
plt.semilogx(freq,gain, label='Experimental')
plt.title(str('Bode amplitude of ' + trans + ', $R_2$ =' + resistor + ' k$\Omega$'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.legend()
#plt.grid()

plt.figure(1)
plt.semilogx(freq,phase, label='Experimental')
plt.title(str('Bode phase of ' + trans + ', $R_2$ =' + resistor + ' k$\Omega$'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shift (deg)')
plt.legend()
#plt.grid()

plt.figure(2)
plt.plot(x,y, label='Experimental')
#plt.plot(x[0], y[0], 'rx', label=str(interval[0] + interval[1] + ' Hz'))
plt.plot(x[-1], y[-1], 'ro', label=str(interval[3] + interval[4] + ' kHz'))
plt.title(str('Nyquist plot of ' + trans + ', $R_2$ =' + resistor + ' kHz'))
plt.legend()
#plt.grid()

plt.show()