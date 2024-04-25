import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import transferfuncs

node = 'V2'
interval = '50-10k'
lower = int(interval[0]+interval[1])
upper = int(interval[3]+interval[4])*1000
match node:
    case 'V2':
        trans = 'T_H'
    case 'V3':
        trans = 'T_B'
    case 'V4':
        trans = 'T_L'
    case _:
        print('Incorrect node name')
xfer = getattr(transferfuncs, str('xfer_' + trans))

### Getting 10 kOhm ###
#Experimental
resistor = '10'
filename = str(node + '_R2-' + resistor + 'k_' + interval + '.txt')
with open(str('C:/Users/sebas/.vscode/experimental_physics/poles_and_zeros/data/' + filename)) as f:
    f.readline()
    f.readline()
    f.readline()
    lines = f.readlines()

freq_10 = []
gain_10 = []
phase_10 = []
for element in lines:
    freq_10.append(float(element.split()[0]))
    gain_10.append(float(element.split()[1]))
    phase_10.append(float(element.split()[2]))

amplitude_10 = [np.power(10,g/20) for g in gain_10]
x_10 = [r * np.cos(np.pi / 180 * f) for r,f in zip(amplitude_10, phase_10)]
y_10 = [r * np.sin(np.pi / 180 * f) for r,f in zip(amplitude_10, phase_10)]

#Theoretical
R2 = 10e3
f_10 = np.logspace(np.log10(lower),np.log10(upper),500)
T_10 = xfer(f_10, R2)

### Getting 33 kOhm
#Experimental
resistor = '33' # In k Ohm
filename = str(node + '_R2-' + resistor + 'k_' + interval + '.txt')
with open(str('C:/Users/sebas/.vscode/experimental_physics/poles_and_zeros/data/' + filename)) as f:
    f.readline()
    f.readline()
    f.readline()
    lines = f.readlines()

freq_33 = []
gain_33 = []
phase_33 = []
for element in lines:
    freq_33.append(float(element.split()[0]))
    gain_33.append(float(element.split()[1]))
    phase_33.append(float(element.split()[2]))

amplitude_33 = [np.power(10,g/20) for g in gain_33]
x_33 = [r * np.cos(np.pi / 180 * f) for r,f in zip(amplitude_33, phase_33)]
y_33 = [r * np.sin(np.pi / 180 * f) for r,f in zip(amplitude_33, phase_33)]

#Theoretical
R2 = 33e3
f_33 = np.logspace(np.log10(lower),np.log10(upper),500)
T_33 = xfer(f_33, R2)


### Getting 150 kOhm
#Experimental
resistor = '150' # In k Ohm
filename = str(node + '_R2-' + resistor + 'k_' + interval + '.txt')
with open(str('C:/Users/sebas/.vscode/experimental_physics/poles_and_zeros/data/' + filename)) as f:
    f.readline()
    f.readline()
    f.readline()
    lines = f.readlines()

freq_150 = []
gain_150 = []
phase_150 = []
for element in lines:
    freq_150.append(float(element.split()[0]))
    gain_150.append(float(element.split()[1]))
    phase_150.append(float(element.split()[2]))

amplitude_150 = [np.power(10,g/20) for g in gain_150]
x_150 = [r * np.cos(np.pi / 180 * f) for r,f in zip(amplitude_150, phase_150)]
y_150 = [r * np.sin(np.pi / 180 * f) for r,f in zip(amplitude_150, phase_150)]

#Theoretical
R2 = 150e3
f_150 = np.logspace(np.log10(lower),np.log10(upper),500)
T_150 = xfer(f_150, R2)


### Plotting
figure_0, axis_0 = plt.subplots(3,3, figsize=(10,12))

# 10 kOhm, column 0
axis_0[0,0].semilogx(freq_10,gain_10, label='Experimental')
axis_0[0,0].semilogx(f_10, 20*np.log10(abs(T_10)), 'black', linestyle=':', label='Theoretical')
axis_0[0,0].set_title('$R_2$ = 10 k$\Omega$')
axis_0[0,0].set_ylabel('Gain (dB)')
axis_0[0,0].grid()
axis_0[0,0].legend()

axis_0[1,0].semilogx(freq_10,phase_10, label='Experimental')
axis_0[1,0].semilogx(f_10, (180/np.pi)*np.unwrap(np.angle(T_10)), 'black', linestyle=':', label='Theoretical')
axis_0[1,0].set_ylabel('Phase shift (deg)')
axis_0[1,0].grid()
axis_0[1,0].legend()

axis_0[2,0].plot(x_10, y_10, label='Experimental')
axis_0[2,0].plot(np.real(T_10), np.imag(T_10), 'black', linestyle=':', label='Theoretical')
axis_0[2,0].grid()
axis_0[2,0].legend()

# 33 kOhm, column 1
axis_0[0,1].semilogx(freq_33,gain_33, label='Experimental')
axis_0[0,1].semilogx(f_33, 20*np.log10(abs(T_33)), 'black', linestyle=':', label='Theoretical')
axis_0[0,1].set_title('$R_2$ = 33 k$\Omega$')
axis_0[0,1].grid()
axis_0[0,1].legend()

axis_0[1,1].semilogx(freq_33,phase_33, label='Experimental')
axis_0[1,1].semilogx(f_33, (180/np.pi)*np.unwrap(np.angle(T_33)), 'black', linestyle=':', label='Theoretical')
axis_0[1,1].grid()
axis_0[1,1].legend()

axis_0[2,1].plot(x_33, y_33, label='Experimental')
axis_0[2,1].plot(np.real(T_33), np.imag(T_33), 'black', linestyle=':', label='Theoretical')
axis_0[2,1].grid()
axis_0[2,1].legend()

# 150 kOhm, column 2
axis_0[0,2].semilogx(freq_150,gain_150, label='Experimental')
axis_0[0,2].semilogx(f_150, 20*np.log10(abs(T_150)), 'black', linestyle=':', label='Theoretical')
axis_0[0,2].set_title('$R_2$ = 150 k$\Omega$')
axis_0[0,2].grid()
axis_0[0,2].legend()

axis_0[1,2].semilogx(freq_150,phase_150, label='Experimental')
axis_0[1,2].semilogx(f_150, (180/np.pi)*np.unwrap(np.angle(T_150)), 'black', linestyle=':', label='Theoretical')
axis_0[1,2].grid()
axis_0[1,2].legend()

axis_0[2,2].plot(x_150, y_150, label='Experimental')
axis_0[2,2].plot(np.real(T_150), np.imag(T_150), 'black', linestyle=':', label='Theoretical')
axis_0[2,2].grid()
axis_0[2,2].legend()

#figure_0.suptitle(str('Bode and Nyquist diagrams for $' + trans + '(i \omega)$'))

plt.figure(1)
figure_1, axis_1 = plt.subplots(1,3)
axis_1[0].plot(x_10,y_10, label='Experimental')
axis_1[0].plot(np.real(T_10), np.imag(T_10), 'black', linestyle=':', label='Theoretical')
axis_1[0].set_title('$R_2$ = 10 k$\Omega$')
axis_1[0].legend()

axis_1[1].plot(x_33,y_33, label='Experimental')
axis_1[1].plot(np.real(T_33), np.imag(T_33), 'black', linestyle=':', label='Theoretical')
axis_1[1].set_title('$R_2$ = 33 k$\Omega$')
axis_1[1].legend()

axis_1[2].plot(x_150,y_150, label='Experimental')
axis_1[2].plot(np.real(T_150), np.imag(T_150), 'black', linestyle=':', label='Theoretical')
axis_1[2].set_title('$R_2$ = 150 k$\Omega$')
axis_1[2].legend()

figure_0.savefig('v2.png')
#plt.show()