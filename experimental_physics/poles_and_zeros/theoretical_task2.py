import numpy as np
import matplotlib.pyplot as plt

import transferfuncs

task = '2.4'
xfer = getattr(transferfuncs, str('xfer_' + task[0] + '_' + task[2]))

f = np.logspace(np.log10(20),np.log10(20e3),500)
T, omega_c = xfer(f)

#Asymptotes
x_as = [20,omega_c[0]/2/np.pi,omega_c[1]/2/np.pi,20000]
y_as = [1,1,(omega_c[0]/omega_c[1]),(omega_c[0]/omega_c[1])]
print(omega_c)
#Phase Guidelines
x_ph = [0.1*omega_c[0]/2/np.pi, 0.1*omega_c[1]/2/np.pi, 10*omega_c[0]/2/np.pi, 10*omega_c[1]/2/np.pi, 20000]
y_ph = [0,-57,-57,0,0]

plt.figure(0)
#plt.semilogx(f, 20*np.log10(abs(T)), 'black', linestyle=':', label='Theoretical')
plt.semilogx(x_as, 20*np.log10((y_as)), 'black', linestyle=':', label='Asymptotes')
plt.title(str('Bode amplitude of task ' + task + ', theoretical'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.grid()
plt.figure(1)
#plt.semilogx(f, (180/np.pi)*np.unwrap(np.angle(T)), 'black', linestyle=':', label='Theoretical')
plt.semilogx(x_ph,y_ph, 'black', linestyle=':', label='Phase guideline')
plt.title(str('Bode phase of task ' + task + ', theoretical'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shift (deg)')
plt.grid()
plt.figure(2)
plt.plot(np.real(T), np.imag(T), 'black', linestyle=':', label='Theoretical')
plt.title(str('Nyquist plot of task ' + task + ', theoretical'))
plt.grid()

#print(np.arctan(np.min(np.imag(T)/np.real(T)))*180/np.pi)

#plt.show()