import numpy as np
import calfem.core as cfc

from plantml import *
from scipy.sparse import lil_matrix

## K-matrix
def create_Km(coords, edof, dofs, elementmarkers, newtBounds, ex, ey, ep, elprop, mark_newt, a_c):
    nDofs = np.size(dofs)
    K_stiff = lil_matrix((nDofs,nDofs))
    D_help = np.eye(2)
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ke = cfc.flw2te(elx, ely, ep, elprop[elMarker][0] * D_help)
        K_stiff[np.ix_(eltopo-1,eltopo-1)] += Ke #Same results as cfc.assem()

    K_c = lil_matrix((nDofs,nDofs))
    for i in newtBounds:
        segment = np.array(i['node-number-list'])
        dist = np.linalg.norm(coords[segment[0]-1]-coords[segment[1]-1]) #Distance between nodes
        fac = ep[0]*a_c*dist/6
        Ke = np.mat([[2*fac, fac], [fac, 2*fac]])
        K_c[np.ix_(segment-1,segment-1)] += Ke #Same results as cfc.assem()
    
    K = K_stiff + K_c
    return K


def create_f(coords, dofs, newtBounds, hBounds, thick, a_c, T_inf, hfunc):
    # Create f_c
    fc = np.zeros([np.size(dofs),1])
    for element in newtBounds:
        nodes = element['node-number-list']
        dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1]) #Distance between nodes
        fc[nodes[0]-1] += 0.5*dist*thick*a_c*T_inf
        fc[nodes[1]-1] += 0.5*dist*thick*a_c*T_inf

    # Create f_b
    fb = np.zeros([np.size(dofs),1])
    for element in hBounds:
        nodes = element['node-number-list']
        dist = np.linalg.norm(coords[nodes[0]-1]-coords[nodes[1]-1]) #Distance between nodes
        fb[nodes[0]-1] += 0.5*dist*thick*hfunc
        fb[nodes[1]-1] += 0.5*dist*thick*hfunc

    f = fb + fc
    return f


def create_C(dofs, edof, ex, ey, elementmarkers, thick, elprop):
    C = lil_matrix((np.size(dofs),np.size(dofs)))
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ce = plantml(elx,ely,thick*elprop[elMarker][1]*elprop[elMarker][2])
        C[np.ix_(eltopo-1,eltopo-1)] += Ce #Same results as cfc.assem()
    return C

def create_K_elastic(edof, dofs, ex, ey, elementmarkers, ep, D):
    K = lil_matrix((np.size(dofs),np.size(dofs)))
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ke = cfc.plante(elx,ely,ep,D[elMarker])
        K[np.ix_(eltopo-1,eltopo-1)] += Ke #Same results as cfc.assem()
    return K

def create_f_elastic(edof, dofs, ex, ey, elementmarkers, ep, D, T, alpha,v):
    f = np.zeros([np.size(dofs),1])
    dT_vec = np.zeros([np.size(edof),1])
    j = 0
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        temps = []
        for i in range(1,6,2):
            temps.append(T[int(eltopo[i]/2-1)])
        dTemp = np.average(temps) - 18 - 273.15
        dT_vec[j] = dTemp
        j+=1

        vector = np.array([1,1,0])
        es = (1+v[elMarker])*dTemp*alpha[elMarker]*D[elMarker]@vector
        fe = cfc.plantf(elx,ely,ep,es)

        for i, f_temp in zip(eltopo, fe):
            f[i-1]+=f_temp
            
    return f, dT_vec # Returns force vector and (constant) temperature on element

def calc_von_Mises(ed,ex,ey,ep,D,elementmarkers,nrNodes,edof,alpha,E,dT_vec,v):
    vM_el = []
    ind = [0,1]
    j=0
    for (eled,elx,ely,elMarker) in zip (ed,ex,ey,elementmarkers):
        es, et = cfc.plants(elx,ely,ep,D[elMarker],eled)
        es[:,ind] -= alpha[elMarker]*E[elMarker]*dT_vec[j] / (1-  2*v[elMarker])
        sigma_zz = v[elMarker] * (es[0,0] + es[0,1]) - alpha[elMarker]*E[elMarker]*dT_vec[j]
        j+=1
        vM_el.append(np.sqrt( pow(es[0,0],2) + pow(es[0,1],2) + pow(sigma_zz,2) - es[0,0]*es[0,1] - es[0,1]*sigma_zz - sigma_zz*es[0,0] + 3*pow(es[0,2],2)))
    
    print(np.max(vM_el))

    trackRecord = np.zeros((nrNodes,2))
    
    for i in range(0,len(vM_el)):
        eltopo = edof[i] ## Notera gammal edof
        for j in range(1,6,2): 
            node = int(eltopo[j]/2)
            trackRecord[node-1,0] += vM_el[i]
            trackRecord[node-1,1] += 1
    
    vM_node = trackRecord[:,0] / trackRecord[:,1]
    return vM_node