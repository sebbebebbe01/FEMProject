import calfem.core as cfc
import numpy as np

from mesh import *

## Material data from table 1
E_cu = 128e9
E_ny = 3e9
v_cu = 0.36
v_ny = 0.39
alpha_cu = 17.6e-6
alpha_ny = 80e-6
rho_cu = 8930
rho_ny = 1100
c_cu = 386
c_ny = 1500

# MATERIAL & MESH PARAMTERS/NAMING BOUNDRY PARTS AND ELEMENT TYPES
# COPPER AREA 1, NYLON AREA 2

# Dimensions (see sketch)
L = 0.005 
a = 0.1 * L
b = 0.1 * L
c = 0.3 * L
d = 0.05 * L
h = 0.15 * L
t = 0.05 * L

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

# Create the K-matrix
nDofs = np.size(dofs)
K_stiff = np.zeros((nDofs,nDofs))
D_help = np.eye(2)
for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
    Ke = cfc.flw2te(elx, ely, ep, elprop[elMarker] * D_help)
    cfc.assem(eltopo, K_stiff, Ke)

K_c = np.zeros((nDofs,nDofs))
newtBounds = boundaryElements[mark_newt]
for i in newtBounds:
    segment = np.array(i['node-number-list'])
    dx = coords[segment[0]-1][0] - coords[segment[1]-1][0]
    dy = coords[segment[0]-1][1] - coords[segment[1]-1][1]
    dist = np.sqrt(dx*dx+dy*dy)
    fac = thick*a_c*dist/6
    Ke = np.mat([[2*fac, fac], [fac, 2*fac]])
    cfc.assem(segment,K_c,Ke)

K = K_stiff + K_c

bc = np.array([],'i')
bcVal = np.array([],'i')

# Experimental f_c
fc = np.zeros([nDofs,1])
for element in newtBounds:
    nodes = element['node-number-list']
    dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1])
    fc[nodes[0]-1] += 1/2*dist*thick*a_c*T_inf
    fc[nodes[1]-1] += 1/2*dist*thick*a_c*T_inf
#print(fc)

# Experimental f_b
hBounds = boundaryElements[mark_h]
fb = np.zeros([nDofs,1])
for element in hBounds:
    nodes = element['node-number-list']
    dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1])
    fb[nodes[0]-1] += 0.5*dist*thick*hfunc
    fb[nodes[1]-1] += 0.5*dist*thick*hfunc
fb *= -1
#print(fb)

f = fb + fc
#print(f)

# Solve the system
T,q = cfc.solveq(K,f,bc,bcVal)

print("Max temperature: ",np.max(T))
print("Min temperature: ",np.min(T))

## VISUALIZE 
mesh.showTemp(T,coords,edof)
