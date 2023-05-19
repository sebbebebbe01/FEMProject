import calfem.core as cfc
import numpy as np

from mesh import *
from matrices import *
import time


# MATERIAL & MESH PARAMTERS/NAMING BOUNDRY PARTS AND ELEMENT TYPES
# COPPER AREA 1, NYLON AREA 2
T_0 = 18 + 273.15
T_inf = 18 + 273.15 #Kelvin
T_max = 143.84823832408833 + 273.15 # From stationary calculations
#T_max = 142.985 + 273.15 #Mac

# Material parameters (only first part for now)
thick = 10#5e-3
ep = [thick]
k1,k2 = 385, 0.26
hfunc = 10**5
a_c = 40
c_1, c_2 = 386, 1500
rho_1, rho_2 = 8930, 1100

# Marker constants
mark_CU = 10
mark_NY = 20
mark_newt = 30
mark_iso = 40
mark_h = 50
# Create a dictionary for different element properties
elprop = {}
elprop[mark_CU] = [k1,c_1,rho_1]
elprop[mark_NY] = [k2,c_2,rho_2]


# Create mesh
mesh = Mesh()
coords, edof, dofs, bdofs, elementmarkers, boundaryElements, ex, ey = mesh.createMesh()

# Arrays with boundary conditions
newtBounds = boundaryElements[mark_newt]
hBounds = boundaryElements[mark_h]

# Create the K-matrix
K = create_Km(coords, edof, dofs, elementmarkers, newtBounds, ex, ey, ep, elprop, mark_newt, a_c)

# Create the f-vectors
f = create_f(coords, dofs, newtBounds, hBounds, thick, a_c, T_inf, hfunc)

# Create C matrix
C = create_C(dofs,edof,ex,ey,elementmarkers,thick,elprop)

# Time interval
t0 = 0
tend = 80
m = 800
dT=(tend-t0)/(m)

M = C+dT*K
T = T_0*np.ones([np.size(dofs),m+1])

t=np.linspace(t0,tend,m+1)

bc = np.array([],'i')
bcVal = np.array([],'i')

a=T_0*np.ones([np.size(dofs),1])
#mesh.showTemp(a,coords,edof)

# Find 90% of T_max
t_step = 0
steps = 0

start = time.time()
while np.max(a)-273.15<0.9*(T_max-273.15): # while max(a)-273.15 < 0.9*(T_max-273.15): För celsius
    rhs = C@a + dT*f
    a,q = cfc.spsolveq(M,rhs,bc,bcVal)
    t_step+=dT
    steps+=1
end = time.time()
print("Time to reach 90% of max temp: ", t_step)
print("Time step: ", dT, " s")
print("Calculation time: ", end-start)
print("Max temp at 90% of max temp: ",np.max(a))
print("Min temp at 90% of max temp: ",np.min(a))

a=T_0*np.ones([np.size(dofs),1])
dT = t_step*0.03 / 5
M = C+dT*K
clim = [291, 300.5] #Colorbar limits, picked manually
for i in range(0,5):
    rhs = C@a + dT*f
    a,q = cfc.spsolveq(M,rhs,bc,bcVal)
    print("Max temp after 3% of 90% of max temp: ",np.max(a))
    print("Min temp after 3% of 90% of max temp: ",np.min(a))
    mesh.showTemp(a,coords,edof,clim)

print("Max temp after 3% of 90% of max temp: ",np.max(a))
print("Min temp after 3% of 90% of max temp: ",np.min(a))
print("Average temp after 3% of 90% of max temp: ", np.average(a))

#mesh.showTemp(a,coords,edof,clim)