import numpy as np
import calfem.core as cfc

from plantml import *

def create_Km(coords, edof, dofs, elementmarkers, newtBounds, ex, ey, ep, elprop, mark_newt, a_c):
    nDofs = np.size(dofs)
    K_stiff = np.zeros((nDofs,nDofs))
    D_help = np.eye(2)
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ke = cfc.flw2te(elx, ely, ep, elprop[elMarker][0] * D_help)
        cfc.assem(eltopo, K_stiff, Ke)

    K_c = np.zeros((nDofs,nDofs))
    for i in newtBounds:
        segment = np.array(i['node-number-list'])
        dx = coords[segment[0]-1][0] - coords[segment[1]-1][0]
        dy = coords[segment[0]-1][1] - coords[segment[1]-1][1]
        dist = np.sqrt(dx*dx+dy*dy)
        fac = ep[0]*a_c*dist/6
        Ke = np.mat([[2*fac, fac], [fac, 2*fac]])
        cfc.assem(segment,K_c,Ke)

    K = K_stiff + K_c
    return K

def create_f(coords, dofs, newtBounds, hBounds, thick, a_c, T_inf, hfunc):
    # Experimental f_c
    fc = np.zeros([np.size(dofs),1])
    for element in newtBounds:
        nodes = element['node-number-list']
        dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1])
        fc[nodes[0]-1] += 1/2*dist*thick*a_c*T_inf
        fc[nodes[1]-1] += 1/2*dist*thick*a_c*T_inf
    #print(fc)

    # Experimental f_b
    fb = np.zeros([np.size(dofs),1])
    for element in hBounds:
        nodes = element['node-number-list']
        dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1])
        fb[nodes[0]-1] += 0.5*dist*thick*hfunc
        fb[nodes[1]-1] += 0.5*dist*thick*hfunc
    fb *= -1
    #print(fb)

    f = fb + fc
    #print(f)
    return f

def create_C(dofs, edof, ex, ey, elementmarkers, thick, elprop):
    C = np.zeros((np.size(dofs),np.size(dofs)))
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ce = plantml(elx,ely,thick*elprop[elMarker][1]*elprop[elMarker][2])
        cfc.assem(eltopo, C, Ce)
    return C