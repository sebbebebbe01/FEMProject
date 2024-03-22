import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

with open('C:/Users/sebas/.vscode/experimental_physics/frequency_doubling/ir_laser_data.txt') as f:
    f.readline() #Rid of first two info lines
    f.readline()
    lines = f.readlines()

ir_amp = []
ir_watt = []
for x in lines:
    column = 3
    try:
        amp = float(x.split()[0])
        watt = float(x.split()[column])
        if amp < 390:
            break
        ir_amp.append(amp)
        ir_watt.append(watt)
    except:
        break

with open('C:/Users/sebas/.vscode/experimental_physics/frequency_doubling/green_dual_data.txt') as f:
    f.readline() #Rid of first two info lines
    f.readline()
    lines = f.readlines()

dual_amp = []
dual_watt = []
for x in lines:
    column = 4
    try:
        amp = float(x.split()[0])
        watt = float(x.split()[column])
        if amp < 280:
            break
        dual_amp.append(amp)
        dual_watt.append(watt)
    except:
        break

with open('C:/Users/sebas/.vscode/experimental_physics/frequency_doubling/green_filtered_data.txt') as f:
    f.readline() #Rid of first two info lines
    f.readline()
    lines = f.readlines()

fil_amp = []
fil_watt = []
for x in lines:
    column = 5
    try:
        amp = float(x.split()[0])
        watt = float(x.split()[column])
        if amp < 280:
            break
        fil_amp.append(amp)
        fil_watt.append(watt)
    except:
        break

with open('C:/Users/sebas/.vscode/experimental_physics/frequency_doubling/rotation_data.txt') as f:
    f.readline() #Rid of first two info lines
    f.readline()
    lines = f.readlines()

angle = []
rot_power = []
for x in lines:
    if float(x.split()[0])==189:
        break
    angle.append(float(x.split()[0]))
    rot_power.append(float(x.split()[1]))

## Simple plots
fig, axs = plt.subplots(3)
axs[0].plot(ir_amp,ir_watt, 'r', label='IR-laser')
axs[1].plot(fil_amp,fil_watt, 'g:', label='Green laser, filtered')
axs[2].plot(dual_amp,dual_watt, 'g', label='Green laser, unfiltered')
for x in axs:
    x.set(xlabel='Diode current (mA)', ylabel='Power (mW)')
    x.legend()
fig.suptitle('Power vs. current')

plt.figure()
plt.plot(ir_amp,ir_watt, 'r', label='IR laser')
plt.plot(dual_amp,dual_watt, 'g', label='Green laser, unfiltered')
plt.plot(fil_amp,fil_watt, 'g:', label='Green laser, filtered')
plt.xlabel('Diode current (mA)')
plt.ylabel('Output power (mW)')
plt.title('Comparison of power measurements')
plt.legend()

plt.figure()
plt.plot(angle,rot_power, 'g', label='Power')
plt.xlabel('Rotation angle (deg)')
plt.ylabel('Output power (mW)')
plt.title('Power vs. incident angle')
plt.legend()

## Linear regression
plt.figure()
ir_reg = stats.linregress(ir_amp, ir_watt)
label = str(round(ir_reg.slope,4)) + '$x-$' + str(-round(ir_reg.intercept,2)) + ', $R^2$: ' + str(round(ir_reg.rvalue,3))
plt.plot(ir_amp,ir_watt, 'r', label='IR-power')
plt.plot(ir_amp, [ir_reg.intercept + ir_reg.slope*x for x in ir_amp], 'black', linestyle=':', label=label)
plt.legend()
plt.xlabel('Diode current (mA)')
plt.ylabel('Output power (mW)')
plt.title('IR power output and line of best fit')

plt.figure()
dual_reg = stats.linregress(dual_amp, dual_watt)
label = str(round(dual_reg.slope,4)) + '$x-$' + str(-round(dual_reg.intercept,2)) + ', $R^2$: ' + str(round(dual_reg.rvalue,3))
plt.plot(dual_amp,dual_watt, 'g', label='Green laser, unfiltered')
plt.plot(dual_amp, [dual_reg.intercept + dual_reg.slope*x for x in dual_amp], 'black', linestyle=':', label=label)
plt.legend()
plt.xlabel('Diode current (mA)')
plt.ylabel('Output power (mW)')
plt.title('Unfiltered green laser power output and line of best fit')

plt.figure()
fil_reg = stats.linregress(fil_amp, fil_watt)
label = str(round(fil_reg.slope,4)) + '$x-$' + str(-round(fil_reg.intercept,2)) + ', $R^2$: ' + str(round(fil_reg.rvalue,3))
plt.plot(fil_amp,fil_watt, 'g', linestyle='dashed', label='Green laser, filtered')
plt.plot(fil_amp, [fil_reg.intercept + fil_reg.slope*x for x in fil_amp], 'black', linestyle=':', label=label)
plt.legend()
plt.xlabel('Diode current (mA)')
plt.ylabel('Output power (mW)')
plt.title('Filtered green laser power output and line of best fit')

## Unfiltered green vs filtered/0.7
calc = [x - y/0.7 for x,y in zip(dual_watt,fil_watt)]
plt.figure()
plt.plot(dual_amp, calc, 'r--', label='Calculated IR-power')
plt.plot(ir_amp,ir_watt, 'r', label='Measured IR-power')
plt.xlabel('Diode current (mA)')
plt.ylabel('Output power (mW)')
plt.title('Difference between filtered and unfiltered power')
plt.legend()

plt.show()