import calfem.core as cfc
import numpy as np
import calfem.utils as cfu

from mesh import *
from matrices import *

## Run stationary file, temperature = T
from stationary import *

## Elasticity constants
E_cu = 128e9
E_ny = 3e9
v_cu = 0.36
v_ny = 0.39
alpha_cu = 17.6e-6
alpha_ny = 80e-6

# Väsentliga randvillkor, se mesh.py
mark_iso_stuck = 60
mark_iso_top = 70
mark_iso_right = 80


## Recreate edof to 2 dofs per node
mesh.dofs_per_node = 2
newEdof = np.zeros((edof.shape[0],6), dtype=int)
for i in range(0,edof.shape[0]):
    eltopo = edof[i]
    for j in range(0,3):
        newEdof[i,2*j]   = eltopo[j]*2 - 1
        newEdof[i,2*j+1] = eltopo[j]*2
edof = newEdof

## Add all dofs
newDofs = np.zeros((dofs.size,2), int)
for i in range(0,dofs.size):
    newDofs[i][0] = dofs[i]*2 - 1
    newDofs[i][1] = dofs[i]*2
dofs = newDofs

## Recreate bdofs
for key in bdofs:
    bdof_list = bdofs[key]
    newBdof_list = []
    for dof in bdof_list:
        newBdof_list.append(2*dof-1)
        newBdof_list.append(2*dof)
    bdofs[key] = newBdof_list

ptype=2
D_cu = cfc.hooke(ptype, E_cu, v_cu)[np.ix_([0,1,3],[0,1,3])]
D_ny = cfc.hooke(ptype, E_ny, v_ny)[np.ix_([0,1,3],[0,1,3])]
ep = [ptype, thick]
D = {}
D[mark_CU] = D_cu
D[mark_NY] = D_ny
E = {}
E[mark_CU] = E_cu
E[mark_NY] = E_ny
alpha = {}
alpha[mark_CU] = alpha_cu
alpha[mark_NY] = alpha_ny
v = {}
v[mark_CU] = v_cu
v[mark_NY] = v_ny

# Create a dictionary for different element properties
elprop = {}
elprop[mark_CU] = D_cu
elprop[mark_NY] = D_ny

# Create K and f
K = create_K_elastic(edof,dofs,ex,ey,elementmarkers,ep,D)
f, dT_vec = create_f_elastic(edof,dofs,ex,ey,elementmarkers,ep,D,T,alpha,v)

# Boundary conditions
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_iso_stuck,0)
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_h,0)
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_iso_top,0,2)
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_iso_right,0,1)

#Solvit
a,q = cfc.spsolveq(K, f, bc, bcVal)

# calculate stresses and strains 
ed = cfc.extractEldisp(edof, a)
von_Mises = calc_von_Mises(ed,ex,ey,ep,D,elementmarkers,dofs.shape[0],edof,alpha,E,dT_vec,v)

print("Max nodal von Mises: ", np.max(von_Mises))
print("Min nodal von Mises: ", np.min(von_Mises))

# Visualize
mesh.vis_von_Mises(von_Mises, coords, newEdof)
mesh.elasticVis(von_Mises, a, coords, newEdof)
