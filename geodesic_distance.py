import numpy as np
import scipy
import trimesh

import laplace

class GeodesicDistance(object):
    def __init__(self, mesh : trimesh.Trimesh):
        self.t = 100 * np.mean(mesh.edges_unique_length)**2
        self.A = laplace.build_mass_matrix(mesh)
        self.Lc = laplace.build_cotangens_matrix(mesh)
        self.M = self.A - self.t*self.Lc
        

        self.lu_M = scipy.sparse.linalg.splu(self.M)
        self.lu_Lc = scipy.sparse.linalg.splu(self.Lc)

        self.num_vertices = len(mesh.vertices)

    def get_distances(self, idx):
        K = np.zeros((self.num_vertices, ), dtype=self.A.dtype)
        K[idx] = 1.0

        L = self.lu_M.solve(K)