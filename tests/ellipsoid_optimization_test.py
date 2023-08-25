"""
aclf_test.py
Description:
    Tests some of the methods defined in aclf.py
"""

import unittest

import time

import polytope as pc
import gurobipy as gp
import pypoman as ppm
import numpy as np

import sys
sys.path.append('../')
from src.aclf import find_ellipsoid_for_intermediate_set

from pydrake.all import HPolyhedron, MaximumVolumeInscribedEllipsoid


class FindingEllipsoids2(unittest.TestCase):

    """
    test_optimization1
    Description:
        Tests the optimization with a value that should be feasible.
    """
    def test_optimization1(self):
        # Constants
        X_interm = HPolyhedron(
            np.vstack((np.eye(2),-np.eye(2))),
            np.array([1.4,1.3,1.2,1.1])
        )
        Theta = HPolyhedron( np.array([[1.0],[-1.0]]) , np.array([0.8,-0.5]) )

        print("X_interm =")
        print(X_interm.A())
        print(X_interm.b())

        ellipsoid_out = MaximumVolumeInscribedEllipsoid(X_interm)

        return

        Theta_vertices = ppm.compute_polytope_vertices(Theta.A,Theta.b)
        theta_dim = len(Theta_vertices[0])

        dim = x_dim + theta_dim

        # Algorithm

        m = gp.Model()

        # Create Optimizaiton Variables
        P = m.addMVar((dim,dim))
        c = m.addMVar((dim,1))
        Q = m.addMVar((dim,dim))

        z_vertices = []
        for vertex_idx in range(1,len(Theta_vertices)):
            z_vertices.append( m.addMVar((dim,1),lb=-gp.GRB.INFINITY) )

        # Add Some Constraints
        # 1. Containment BY the polytope X_interm
        num_rows = len(X_interm.A)
        for row_ind in range(num_rows):
            A_i = X_interm.A[row_ind]
            b_i = X_interm.b[row_ind]
            A_i_tilde = np.hstack(
                (np.reshape(A_i,(1,x_dim)),np.zeros((1,theta_dim)))
            )
            m.addConstr(
                A_i_tilde @ Q @ A_i_tilde.T <= b_i - A_i_tilde @ c
            )
        
        # Try To Solve Model
        m.write("logs/find_ellipsoid-test2.rlp")
        start = time.time()
        m.optimize()
        end = time.time()
        print('solving it takes %.3f s'%(end - start))

        self.assertTrue(end-start > 0)
        for row_index in range(dim):
            for col_index in range(dim):
                self.assertTrue(P.X[row_index,col_index] >= 0)

if __name__ == '__main__':
    unittest.main()