import calfem.core as cfc
import numpy as np

from mesh import *
from matrices import *
import time


# Temperatures
T_inf = 18 + 273.15 #Kelvin, temperature outside of the body
#T_max, el_size_factor = 143.78121812271985 + 273.15, 0.2 # Calculated from the stationary problem using el_size_factor = 0.2
#T_max, el_size_factor = 143.84823832408833 + 273.15, 0.02 
T_max, el_size_factor = 143.85509829273457 + 273.15, 0.01 
#T_max, el_size_factor = 143.85819010827964 + 273.15, 0.005

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
mesh = Mesh(el_size_factor=el_size_factor)
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

# Boundary vectors (empty)
bc = np.array([],'i')
bcVal = np.array([],'i')

# Time interval
t0 = 0
tend = 80
m = 800 #Nr of steps
delta_t=(tend-t0)/(m) #Step size

a=T_inf*np.ones([np.size(dofs),1]) # Initial temperature = 18 celsius

## Find 90 % of T_max
t_step = 0
prev_max = T_inf # Used to find a more accurate t_90
M = C+delta_t*K #Help matrix
start = time.time() # Measure computational time
while np.max(a)-273.15 < 0.9*(T_max-273.15):
    rhs = C@a + delta_t*f
    a,q = cfc.spsolveq(M,rhs,bc,bcVal)
    t_step+=delta_t
    new_max = np.max(a)
    delta_temp = new_max - prev_max
    prev_max = new_max
end = time.time()
#Find a more accurate t_90 by linearizing the relation bewtween the last and second to last temps
t_90 = t_step - (np.max(a)-273.15 - 0.9*(T_max-273.15))*delta_t/delta_temp
print("Time to reach 90% of max temp: ", t_90)
print("Time step: ", delta_t, " s")
print("Calculation time: ", end-start)
#print("Max temp at 90% of max temp: ",np.max(a))
#print("Min temp at 90% of max temp: ",np.min(a))

## Show 5 evenely spread snapshots of the first 3 %  of the time it took to reach 90 % of T_max
a=T_inf*np.ones([np.size(dofs),1]) #Reset temperature vector
delta_t = t_90*0.03 / 5 #New step size
M = C+delta_t*K #Help matrix
clim = [291, 300.5] #Colorbar limits, picked manually
for i in range(0,5):
    rhs = C@a + delta_t*f
    a,q = cfc.spsolveq(M,rhs,bc,bcVal)
    mesh.showTemp(a,coords,edof,clim)

print("Max temp after 3% of 90% of max temp: ",np.max(a))
print("Min temp after 3% of 90% of max temp: ",np.min(a))
print("Average temp after 3% of 90% of max temp: ", np.average(a))
