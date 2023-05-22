import numpy as np
import calfem.core as cfc

import calfem.geometry as cfg
import calfem.mesh as cfm
from calfem.qt5 import *
import calfem.vis as cfv

class Mesh:
    def __init__(self, el_type=2, dofs_per_node=1, el_size_factor=0.02):
        self.el_type, self.dofs_per_node, self.el_size_factor = (
            el_type,
            dofs_per_node,
            el_size_factor,
        )

    def createMesh(self):
        """
        Creates the mesh with the dimensions given in the problem description

        Returns:
            coords, edof, dofs, bdofs, elementmarkers, boundaryElements, ex, ey
        """

        # Dimensions (see sketch)
        L = 0.005
        a = 0.1 * L
        b = 0.1 * L
        c = 0.3 * L
        d = 0.05 * L
        h = 0.15 * L
        t = 0.05 * L

        # Marker constants
        mark_CU = 10 #Copper
        mark_NY = 20 #Nylon
        mark_newt = 30 #Boundaries following Newton's law of cooling
        mark_iso = 40 #Homogeneous Dirichlet conditions (isolated), these boundaries are also fastened (no displacement)
        mark_h = 50 #Constant heat flow boundary condition
        mark_iso_stuck = 60 #Isolated and fastened
        mark_top_mirror = 70 #Homogeneous Dirichlet conditions in y-direction (due to mirror-symmetry)
        mark_right_mirror = 80 #Homogeneous Dirichlet conditions in x-direction (due to mirror-symmetry)

        el_type, dofs_per_node, el_size_factor = (
            self.el_type,
            self.dofs_per_node,
            self.el_size_factor,
        )

        ## GEOMETRY PART AND VISUALISATION

        # Create geometry
        g = cfg.Geometry()
        self.g = g

        # Add points
        points = [
            [-L, -b - a],
            [-L, -0.5 * L],
            [-L + c + d, -0.5 * L],
            [-L + c + d, -b - a],
            [-L + a + t, -b - a],
            [-L + a + t, -b - a - h],
            [-L + a, -b - a - h],
            [-L + a, -b - a],
            [-L, -b],
            [-L, 0],
            [-L + a, 0],
            [-L + a, -b],
            [-L + a + c, -b],
            [-L + a + c + d, -b - d],
            [-L + a + c + d, -0.5 * L + d], 
            [-2 * d, -0.2 * L],
            [0, -0.2 * L],
            [0, -0.2 * L - d],
            [-2 * d, -0.2 * L - d],
            [-L + a + c + d, -0.5 * L],
        ]

        for x, y in points:
            g.point([x, y])

        # Add lines
        g.spline([0, 1], marker=mark_iso_stuck)
        g.spline([1, 2], marker=mark_iso)
        g.spline([2, 3])
        g.spline([3, 4])
        g.spline([4, 5])
        g.spline([5, 6])
        g.spline([6, 7])
        g.spline([7, 0])
        g.spline([0, 8], marker=mark_iso_stuck)
        g.spline([8, 9], marker=mark_h)
        g.spline([9, 10], marker=mark_top_mirror)  # because of symmetry
        g.spline([10, 11], marker=mark_newt)
        g.spline([11, 12], marker=mark_newt)
        g.spline([12, 13], marker=mark_newt)
        g.spline([13, 14], marker=mark_newt)
        g.spline([14, 15], marker=mark_newt)
        g.spline([15, 16], marker=mark_newt)
        g.spline([16, 17], marker=mark_right_mirror)  # because of symmetry
        g.spline([17, 18], marker=mark_newt)
        g.spline([18, 19], marker=mark_newt)
        g.spline([19, 2], marker=mark_iso)

        # Add surfaces
        g.surface(
            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            marker=mark_CU,
        )
        g.surface([0, 1, 2, 3, 4, 5, 6, 7], marker=mark_NY)

        # Create mesh
        mesh = cfm.GmshMeshGenerator(g)
        mesh.el_size_factor = el_size_factor
        mesh.el_type = el_type
        mesh.dofs_per_node = dofs_per_node
        mesh.return_boundary_elements = True
        coords, edof, dofs, bdofs, elementmarkers, boundaryElements = mesh.create()
        ex, ey = cfc.coordxtr(edof, coords, dofs)
        self.coords, self.edof = coords, edof

        return coords, edof, dofs, bdofs, elementmarkers, boundaryElements, ex, ey

    def showSkeleton(self):
        """
        Draws the geometry of the object and its mesh
        """
        cfv.figure()
        cfv.drawGeometry(self.g, title="Geometry")

        cfv.figure()
        cfv.draw_mesh(
            coords=self.coords,
            edof=self.edof,
            dofs_per_node=self.dofs_per_node,
            el_type=self.el_type,
            filled=True,
            title="Mesh (units in mm)",
        )
        cfv.showAndWait()

    def showTemp(self, T, coords, edof, clim = None):
        """
        Mirrors and draws the temperature distribution of the object with (optionally) color bar limits clim = [a,b]
        """
        el_type = self.el_type  # Type of mesh
        dofs_per_node = self.dofs_per_node  # Factor that changes element sizes

        cfv.figure()
        cfv.draw_nodal_values(
            T,
            coords,
            edof,
            dofs_per_node,
            el_type,
            draw_elements=False,
            title="Temperature",
            clim=clim
        )
        cfv.draw_nodal_values(
            T,
            np.transpose(np.array((coords[:, 0], -coords[:, 1]))),
            edof,
            dofs_per_node,
            el_type,
            draw_elements=False,
            title="Temperature",
            clim=clim
        )
        cfv.draw_nodal_values(
            T,
            -coords,
            edof,
            dofs_per_node,
            el_type,
            draw_elements=False,
            title="Temperature",
            clim=clim
        )
        cfv.draw_nodal_values(
            T,
            np.transpose(np.array((-coords[:, 0], coords[:, 1]))),
            edof,
            dofs_per_node,
            el_type,
            draw_elements=False,
            title="Temperature",
            clim=clim
        )
        cfv.colorBar()
        cfv.showAndWait()

    def vis_von_Mises(self, von_Mises, coords, edof):
        """
        Mirrors and draws the nodal von Mises strains
        """
        cfv.figure()
        cfv.draw_nodal_values(
            von_Mises,
            coords,
            edof,
            self.dofs_per_node,
            self.el_type,
            title="von Mises",
            draw_elements=False
        )
        cfv.draw_nodal_values(
            von_Mises,
            np.array((coords[:, 0], -coords[:, 1])).T,
            edof,
            self.dofs_per_node,
            self.el_type,
            title="von Mises",
            draw_elements=False
        )
        cfv.draw_nodal_values(
            von_Mises,
            -coords,
            edof,
            self.dofs_per_node,
            self.el_type,
            title="von Mises",
            draw_elements=False
        )
        cfv.draw_nodal_values(
            von_Mises,
            np.transpose(np.array((-coords[:, 0], coords[:, 1]))),
            edof,
            self.dofs_per_node,
            self.el_type,
            title="von Mises",
            draw_elements=False
        )
        cfv.show_and_wait()

    def vis_disp_and_stress(self, von_Mises, a, coords, edof):
        """
        Mirrors and draws the displacements, a, of the nodes with nodal values (von Mises stresses), von_Mises.
        """

        ## Puts the displacement values in a in a (nNodes, 2) array instead of (2*nNodes,1) in order to mirror
        b = np.zeros((int(a.size / 2), 2))
        i = 0
        j = 0
        for x in a:
            b[int(i / 2), j] = x
            i += 1
            j = (j + 1) % 2

        el_type = self.el_type  # Type of mesh
        dofs_per_node = self.dofs_per_node  # Factor that changes element sizes
        magnfac=1.0
        showOldMesh = True
        title = 'Displacements (magnfac = '+ str(magnfac) + ') and stresses'
        von_Mises = von_Mises.tolist()
        cfv.figure()
        cfv.draw_displacements(
            b,
            coords,
            edof,
            dofs_per_node,
            el_type,
            von_Mises,
            draw_undisplaced_mesh=showOldMesh,
            title=title,
            magnfac=magnfac,
        )
        cfv.draw_displacements(
            np.array((b[:, 0], -b[:, 1])).T,
            np.array((coords[:, 0], -coords[:, 1])).T,
            edof,
            dofs_per_node,
            el_type,
            von_Mises,
            draw_undisplaced_mesh=showOldMesh,
            title=title,
            magnfac=magnfac,
        )
        cfv.draw_displacements(
            -b,
            -coords,
            edof,
            dofs_per_node,
            el_type,
            von_Mises,
            draw_undisplaced_mesh=showOldMesh,
            title=title,
            magnfac=magnfac,
        )
        cfv.draw_displacements(
            np.array((-b[:, 0], b[:, 1])).T,
            np.transpose(np.array((-coords[:, 0], coords[:, 1]))),
            edof,
            dofs_per_node,
            el_type,
            von_Mises,
            draw_undisplaced_mesh=showOldMesh,
            title=title,
            magnfac=magnfac,
        )
        cfv.show_and_wait()
