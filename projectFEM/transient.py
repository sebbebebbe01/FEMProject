import calfem.core as cfc
import numpy as np

from mesh import *
from matrices import *


# MATERIAL & MESH PARAMTERS/NAMING BOUNDRY PARTS AND ELEMENT TYPES
# COPPER AREA 1, NYLON AREA 2
T_0 = 18 + 273.15
T_inf = 18 + 273.15 #Kelvin

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
K = create_Km(coords, edof, dofs, elementmarkers, boundaryElements, ex, ey, ep, elprop, mark_newt, a_c)

# Create the f-vectors
f = create_f(coords, dofs, newtBounds, hBounds, thick, a_c, T_inf, hfunc)

# Create C matrix
C = create_C(dofs,edof,ex,ey,elementmarkers,thick,elprop)

# Time interval
t0 = 0
tend = 1000
m = 1000
dT=(tend-t0)/(m)

M = C+dT*K
T = T_0*np.ones([np.size(dofs),m+1])

t=np.linspace(t0,tend,m+1)

bc = np.array([],'i')
bcVal = np.array([],'i')

a=T_0*np.ones([np.size(dofs),1])
#mesh.showTemp(a,coords,edof)

for i in range(0,m):
    rhs = C@a + dT*f
    #print(rhs.shape)
    a,q = cfc.solveq(M,rhs,bc,bcVal)
    #T[:,[i+1]] = a # T[:,i+1] = a.ravel() funkar också
    mesh.showTemp(a,coords,edof)

#mesh.animate(T, coords, edof)