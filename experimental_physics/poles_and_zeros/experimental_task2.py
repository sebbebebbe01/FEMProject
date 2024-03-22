import numpy as np
import matplotlib.pyplot as plt

task = '2.4'
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
plt.semilogx(freq,gain)
plt.title(str('Bode amplitude of task ' + task + ', experimental'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.grid()

plt.figure(1)
plt.semilogx(freq,phase)
plt.title(str('Bode phase of task ' + task + ', experimental'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shifft (deg)')
plt.grid()

plt.figure(2)
plt.plot(x,y)
plt.title(str('Nyquist plot of task ' + task + ', experimental'))
plt.grid()

plt.show()