import numpy as np
import matplotlib.pyplot as plt

import theoretical_task2

task = '2.4'
name = '(d)'
interval = '20-20k'
filename = str(task[0] + '_' + task[2] + '_bode_20-20k.txt')
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

plt.figure(0)
plt.semilogx(freq,gain, label='Experimental')
#plt.title(str('Bode amplitude of Circuit ' + name))
plt.title(str('Experimental gain of Circuit ' + name + ' with asymptotes'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.legend()
#plt.grid()

plt.figure(1)
plt.semilogx(freq,phase, label='Experimental')
#plt.title(str('Bode phase of Circuit ' + name))
plt.title(str('Experimental phase shift of Circuit ' + name + ' with phase guidelines'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shift (deg)')
plt.legend()
#plt.grid()

plt.figure(2)
plt.plot(x,y, label='Experimental')
plt.plot(x[-1], y[-1], 'ro', label=str(interval[3] + interval[4] + ' kHz'))
plt.title(str('Nyquist plot of Circuit ' + name))
plt.legend()
#plt.grid()


#print(np.arctan(np.min(np.array(y)/np.array(x)))*180/np.pi)

plt.show()