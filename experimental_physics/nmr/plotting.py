import  numpy  as  np
import matplotlib.pyplot as plt
from scipy import optimize
import fittingfunction

data = 'T1BpWater'
with open(str('C:/Users/sebas/.vscode/experimental_physics/nmr/data/' + data + '.txt')) as f:
    f.readline() #Rid of first three info lines
    f.readline()
    f.readline()
    lines = f.readlines()

delay = []
amplitude = []
for x in lines:
    delay.append(float(x.split()[0])/1000) #In seconds
    amplitude.append(float(x.split()[1]))
delay = np.array(delay)
amplitude = np.array(amplitude)

func = getattr(fittingfunction, data[:4])
[curve_of_best_fit, popt, pcov] = func(delay, amplitude)



uncert = np.sqrt(np.diag(pcov))
print("Std: " + str(uncert))

timeType = str('$T_' + data[1] + '$')
legendNrs = tuple(np.array((popt[0], uncert[0], popt[1], uncert[1])))
if data[2]=='B':
    title = str('Finding ' + timeType + ', using B$_' + data[3] + '$')
else:
    title = 'Finding ' + timeType + ' of ' + data[4:]
plt.plot(delay, amplitude)
plt.plot(delay, curve_of_best_fit, 'r-',
         label=str('fit: $(S_0$=%5.3f$\pm$%5.3f) $\mu$V, ' +timeType+ '=(%5.3f$\pm$%5.3f) s') % legendNrs)
plt.legend()
plt.title(title)
plt.xlabel('t (s)')
plt.ylabel('Attenuation (E/E$_0$)')
plt.grid()
plt.show()