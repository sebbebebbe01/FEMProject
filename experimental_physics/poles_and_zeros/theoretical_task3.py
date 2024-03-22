import numpy as np
import matplotlib.pyplot as plt
import transferfuncs

trans = 'T_H'
xfer = getattr(transferfuncs, str('xfer_' + trans))

R2 = 10e3

f = np.logspace(np.log10(10),np.log10(20e3),500)
T = xfer(f, R2)

label = str('$'+trans+'$')
plt.figure(0)
plt.semilogx(f, 20*np.log10(abs(T)))
plt.title(str('Bode amplitude of task ' + label + ', theoretical'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain (dB)')
plt.grid()
plt.figure(1)
plt.semilogx(f, (180/np.pi)*np.unwrap(np.angle(T)))
plt.title(str('Bode phase of task ' + label + ', theoretical'))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase shift (deg)')
plt.grid()
plt.figure(2)
plt.plot(np.real(T), np.imag(T))
plt.title(str('Nyquist plot of task ' + label + ', theoretical'))
plt.grid()
plt.show()