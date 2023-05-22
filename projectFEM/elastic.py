import calfem.core as cfc
import numpy as np
import calfem.utils as cfu

from mesh import *
from matrices import *

## Run stationary file, temperature vector = T
from stationary import *

## Elasticity constants
E_cu = 128e9
E_ny = 3e9
nu_cu = 0.36
nu_ny = 0.39
alpha_cu = 17.6e-6
alpha_ny = 80e-6

ptype=2 # Plane strain
ep = [ptype, t_d]

## Create dictionaries for element properties
D_cu = cfc.hooke(ptype, E_cu, nu_cu)[np.ix_([0,1,3],[0,1,3])]
D_ny = cfc.hooke(ptype, E_ny, nu_ny)[np.ix_([0,1,3],[0,1,3])]
D = {}
D[mark_CU] = D_cu
D[mark_NY] = D_ny
E = {}
E[mark_CU] = E_cu
E[mark_NY] = E_ny
alpha = {}
alpha[mark_CU] = alpha_cu
alpha[mark_NY] = alpha_ny
nu = {}
nu[mark_CU] = nu_cu
nu[mark_NY] = nu_ny


# Väsentliga randvillkor, se mesh.py
mark_iso_stuck = 60
mark_top_mirror = 70
mark_right_mirror = 80


## Recreate edof to 2 dofs per node, each node will get it's old dof replaced by 2 new ones
mesh.dofs_per_node = 2
newEdof = np.zeros((edof.shape[0],6), dtype=int)
for i in range(0,edof.shape[0]):
    eltopo = edof[i]
    for j in range(0,3):
        newEdof[i,2*j]   = eltopo[j]*2 - 1
        newEdof[i,2*j+1] = eltopo[j]*2

#Possibly quicker, but calfem.vis.draw_nodal_values() complains about something. It is identical to the one produced above to both size, type, and values
#newEdof = np.array((2*edof[:,0]-1, 2*edof[:,0], 2*edof[:,1] - 1, 2*edof[:,1], 2*edof[:,2] - 1, 2*edof[:,2])).T

## Add all dofs to reflect the 2 dofs per node
newDofs = np.zeros((dofs.size,2), int)
for i in range(0,dofs.size):
    newDofs[i][0] = dofs[i]*2 - 1
    newDofs[i][1] = dofs[i]*2
dofs = newDofs

## Recreate bdofs, mark all new dofs with the correct boundary condition
for key in bdofs:
    bdof_list = bdofs[key]
    newBdof_list = []
    for dof in bdof_list:
        newBdof_list.append(2*dof-1)
        newBdof_list.append(2*dof)
    bdofs[key] = newBdof_list

# Create K and f
K = create_K_elastic(newEdof,dofs,ex,ey,elementmarkers,ep,D)
f, dT_vec = create_f_elastic(newEdof,dofs,ex,ey,elementmarkers,ep,D,T,alpha,nu)

# Boundary conditions
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_iso_stuck,0) # Fastened boundaries, displacement=0
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_h,0) # Boundaries marked h (constant heat flow) are also fastened
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_top_mirror,0,2) # Mirrored at the top, y-displacement=0
bc,bcVal = cfu.applybc(bdofs,bc,bcVal,mark_right_mirror,0,1) # Mirrored to the right, x-displacement=0

#Solvit
a,q = cfc.spsolveq(K, f, bc, bcVal)

# calculate stresses and strains 
ed = cfc.extractEldisp(newEdof, a)
von_Mises = calc_von_Mises(ed,ex,ey,ep,D,elementmarkers,dofs.shape[0],edof,alpha,E,dT_vec,nu) # edof for the thermal, not elastic, problem

print("Max nodal von Mises: ", np.max(von_Mises))
print("Min nodal von Mises: ", np.min(von_Mises))

# Visualize
mesh.vis_von_Mises(von_Mises, coords, newEdof)
mesh.vis_disp_and_stress(von_Mises, a, coords, newEdof)
