import numpy as np
import calfem.core as cfc

from plantml import *
from scipy.sparse import lil_matrix

def create_Km(coords, edof, dofs, elementmarkers, newtBounds, ex, ey, ep, elprop, a_c):
    """
    Creates the K-matrix for the temperature distribution using K_stiff
    for all elements and K_c for all boundary elements with Newton conditions.
    Returns:
        K
    """
    nDofs = np.size(dofs)
    K_stiff = lil_matrix((nDofs, nDofs))
    D_help = np.eye(2)
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ke = cfc.flw2te(elx, ely, ep, elprop[elMarker][0] * D_help)
        K_stiff[np.ix_(eltopo - 1, eltopo - 1)] += Ke  # Same results as cfc.assem()

    K_c = lil_matrix((nDofs, nDofs))
    for element in newtBounds:
        nodes = np.array(
            element["node-number-list"]
        )  # Extract nodes of line-segment and converts it to an array for np.ix_ to work
        dist = np.linalg.norm(
            coords[nodes[0] - 1] - coords[nodes[1] - 1]
        )  # Distance between nodes
        help_factor = ep[0] * a_c * dist / 6  # Manually derived, see report
        Ke = help_factor * np.mat([[2, 1], [1, 2]])
        K_c[np.ix_(nodes - 1, nodes - 1)] += Ke  # Same results as cfc.assem()

    K = K_stiff + K_c
    return K


def create_f(coords, dofs, newtBounds, hBounds, t_d, a_c, T_inf, hfunc):
    """
    Creates the boundary force vector for the temperature distribution using f_c
    for all boundary elements with Newton conditions and f_b for all elements with
    a given heat insertion (h).
    Returns:
        f
    """
    # Create f_c, see report for motivation
    f_c = np.zeros([np.size(dofs), 1])
    for element in newtBounds:
        nodes = element["node-number-list"]  # Extract nodes of line-segment
        dist = np.linalg.norm(
            coords[nodes[0] - 1] - coords[nodes[1] - 1]
        )  # Distance between nodes
        f_c[nodes[0] - 1] += 0.5 * dist * t_d * a_c * T_inf
        f_c[nodes[1] - 1] += 0.5 * dist * t_d * a_c * T_inf

    # Create f_b, see report for motivation
    f_b = np.zeros([np.size(dofs), 1])
    for element in hBounds:
        nodes = element["node-number-list"]
        dist = np.linalg.norm(
            coords[nodes[0] - 1] - coords[nodes[1] - 1]
        )  # Distance between nodes
        f_b[nodes[0] - 1] += 0.5 * dist * t_d * hfunc
        f_b[nodes[1] - 1] += 0.5 * dist * t_d * hfunc

    f = f_b + f_c
    return f


def create_C(dofs, edof, ex, ey, elementmarkers, t_d, elprop):
    """
    Creates the transient temperature matrix C.
    Returns:
        C
    """
    C = lil_matrix((np.size(dofs), np.size(dofs)))
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ce = plantml(
            elx, ely, t_d * elprop[elMarker][1] * elprop[elMarker][2]
        )  # plantml from canvas, not calfem.core
        C[np.ix_(eltopo - 1, eltopo - 1)] += Ce  # Same results as cfc.assem()
    return C

def create_K_elastic(edof, dofs, ex, ey, elementmarkers, ep, D):
    """
    Creates the K-matrix for the elastic problem.
    Returns:
        K
    """
    K = lil_matrix((np.size(dofs), np.size(dofs)))
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        Ke = cfc.plante(elx, ely, ep, D[elMarker])
        K[np.ix_(eltopo - 1, eltopo - 1)] += Ke  # Same results as cfc.assem()
    return K


def create_f_elastic(edof, dofs, ex, ey, elementmarkers, ep, D, T, alpha, nu):
    """
    Creates the element force vector f for the elastic problem as well as the elemental
    temperatures.
    Returns:
        f, dT_vec
    """
    f = np.zeros([np.size(dofs), 1])
    dT_vec = np.zeros([np.size(edof), 1])
    j = 0
    for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
        # Calculating delta T for given element by taking the mean of corresponding nodal values.
        temps = []
        for i in range(1, 6, 2):
            temps.append(T[int(eltopo[i] / 2 - 1)])
        dTemp = (
            np.average(temps) - 18 - 273.15
        )  # T_inf = 18 C, we want temperature difference so we substract innitial temp (in Kelvin)
        dT_vec[j] = dTemp  # ...and saving values of delta T
        j += 1

        # Compute internal element force vector
        vector = np.array([1, 1, 0])
        es = (1 + nu[elMarker]) * dTemp * alpha[elMarker] * D[elMarker] @ vector
        fe = cfc.plantf(elx, ely, ep, es)

        # Add contribution to total force vector
        for i, f_temp in zip(eltopo, fe):
            f[i - 1] += f_temp

    return f, dT_vec  # Returns force vector and temperature on element


def calc_von_Mises(
    ed, ex, ey, ep, D, elementmarkers, nrNodes, edof, alpha, E, dT_vec, nu
):
    """
    Creates the nodal von Mises-stress vector by first creating the elemental stresses
    from core.plants() and then averaging each nodes adjacent elemental values onto the node.
    Returns:
        vM_node
    """
    # Start by calculating elemental von Mises stresses
    vM_el = []
    ind = [
        0,
        1,
    ]  # We just want to substract thermal stresses from total normal stresses
    j = 0
    max_stress={}
    max_stress[10]=0
    max_stress[20]=0
    for eled, elx, ely, elMarker in zip(ed, ex, ey, elementmarkers):
        es, et = cfc.plants(elx, ely, ep, D[elMarker], eled)
        es[:, ind] -= (
            alpha[elMarker] * E[elMarker] * dT_vec[j] / (1 - 2 * nu[elMarker])
        )  # See report
        sigma_zz = (
            nu[elMarker] * (es[0, 0] + es[0, 1])
            - alpha[elMarker] * E[elMarker] * dT_vec[j]
        )  # Manual calculation of sigma_zz (not given by plants)
        j += 1
        vM = np.sqrt(
                pow(es[0, 0], 2)
                + pow(es[0, 1], 2)
                + pow(sigma_zz, 2)
                - es[0, 0] * es[0, 1]
                - es[0, 1] * sigma_zz
                - sigma_zz * es[0, 0]
                + 3 * pow(es[0, 2], 2)
        ) # See report
        vM_el.append(vM)  
        if vM > max_stress[elMarker]:
            max_stress[elMarker] = vM

    print("Max elemental stress, copper: ",max_stress[10][0])
    print("Max elemental stress,  nylon: ",max_stress[20][0])

    # Add neighboring elemental values of a node to a help matrix and divide by the number of adjacent elements to get nodal stesses
    nodal_vM_help = np.zeros((nrNodes, 2))
    for i in range(0, len(vM_el)):
        eltopo = edof[i]  # Note this is the edof for the thermal, not elastic, problem
        for j in range(3):
            node = eltopo[j]
            nodal_vM_help[node - 1, 0] += vM_el[i]
            nodal_vM_help[node - 1, 1] += 1

    vM_node = nodal_vM_help[:, 0] / nodal_vM_help[:, 1]
    return vM_node