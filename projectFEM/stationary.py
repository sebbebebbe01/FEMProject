import calfem.core as cfc
import numpy as np

from mesh import *
from matrices import *

# MATERIAL & MESH PARAMTERS/NAMING BOUNDRY PARTS AND ELEMENT TYPES
# COPPER AREA 1, NYLON AREA 2
T_inf = 18 + 273.15 #Kelvin

# Material parameters (only first part for now)
thick = 5e-3
ep = [thick]
k1,k2 = 385, 0.26
hfunc = 10**5
a_c = 40

# Marker constants
mark_CU = 10
mark_NY = 20
mark_newt = 30
mark_iso = 40
mark_h = 50
# Create a dictionary for different element properties
elprop = {}
elprop[mark_CU] = [k1] 
elprop[mark_NY] = [k2]

# Create mesh
mesh = Mesh()
coords, edof, dofs, bdofs, elementmarkers, boundaryElements, ex, ey = mesh.createMesh()
#mesh.showSkeleton()

# Arrays with boundary conditions
newtBounds = boundaryElements[mark_newt]
hBounds = boundaryElements[mark_h]

# Create the K-matrix
K = create_Km(coords, edof, dofs, elementmarkers, newtBounds, ex, ey, ep, elprop, mark_newt, a_c)

# Create the f-vectors
f = create_f(coords, dofs, newtBounds, hBounds, thick, a_c, T_inf, hfunc)

# Solve the system
bc = np.array([],'i')
bcVal = np.array([],'i')
T,q = cfc.spsolveq(K,f,bc,bcVal)

print("Max temperature: ",np.max(T)-273.15)
print("Min temperature: ",np.min(T)-273.15)

## VISUALIZE 
#mesh.showTemp(T,coords,edof)
