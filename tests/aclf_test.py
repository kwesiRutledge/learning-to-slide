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


class FindingEllipsoids(unittest.TestCase):

    def test_feasibility1(self):
        # Constants
        P_interm = pc.Polytope(
            np.vstack((np.eye(3),-np.eye(3))),
            np.array([1.0,2.0,3.0,1.2,2.2,-0.2])
        )
        Theta = pc.Polytope( np.array([[1.0],[-1.0]]) , np.array([2.0,-0.5]) )

        dim = P_interm.dim

        Theta_vertices = ppm.compute_polytope_vertices(Theta.A,Theta.b)

        # Algorithm

        m = gp.Model()

        # Create Optimizaiton Variables
        P = m.addMVar((dim,dim))
        c = m.addMVar((dim,1))

        # Add Some Constraints

        # Try To Solve Model
        m.write("logs/find_ellipsoid.rlp")
        start = time.time()
        m.optimize()
        end = time.time()
        print('solving it takes %.3f s'%(end - start))

        self.assertTrue(end-start > 0)
        for row_index in range(dim):
            for col_index in range(dim):
                self.assertTrue(P.X[row_index,col_index] >= 0)

    def test_feasibility2(self):
        # Constants
        P_interm = pc.Polytope(
            np.vstack((np.eye(3),-np.eye(3))),
            np.array([1.0,2.0,3.0,1.2,2.2,-0.2])
        )
        Theta = pc.Polytope( np.array([[1.0],[-1.0]]) , np.array([2.0,-0.5]) )

        dim = P_interm.dim

        Theta_vertices = ppm.compute_polytope_vertices(Theta.A,Theta.b)

        # Algorithm

        m = gp.Model()

        # Create Optimizaiton Variables
        P = m.addMVar((dim,dim))
        c = m.addMVar((dim,1))

        z_vertices = []
        for vertex_idx in range(1,len(Theta_vertices)):
            z_vertices.append( m.addMVar((dim,1),lb=-gp.GRB.INFINITY) )

        # Add Some Constraints

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

    """
    test_optimization1
    Description:
        Tests the optimization with a value that should be feasible.
    """
    def test_optimization1(self):
        # Constants
        X_interm = pc.Polytope(
            np.vstack((np.eye(2),-np.eye(2))),
            np.array([1.4,1.3,1.2,1.1])
        )
        Theta = pc.Polytope( np.array([[1.0],[-1.0]]) , np.array([0.8,-0.5]) )

        x_dim = X_interm.dim

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


    """
    test_optimization2
    Description:
        Tests the optimization with the nonconvex constraints.
    """
    def test_optimization2(self):
        # Constants
        X_interm = pc.Polytope(
            np.vstack((np.eye(2),-np.eye(2))),
            np.array([1.4,1.3,1.2,1.1])
        )
        Theta = pc.Polytope( np.array([[1.0],[-1.0]]) , np.array([0.8,-0.5]) )

        x_dim = X_interm.dim

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
        for vertex_idx in range(len(Theta_vertices)):
            z_vertices.append( m.addMVar((dim,1),lb=-gp.GRB.INFINITY) )

        # Add Some Constraints
        # 1. Containment BY the polytope X_interm
        num_rows = len(X_interm.A)
        for row_ind in range(num_rows):
            A_i = X_interm.A[row_ind]
            b_i = X_interm.b[row_ind]
            A_i_tilde = np.hstack(
                ( np.reshape(A_i,(1,x_dim)),np.zeros((1,theta_dim)) )
            )
            # m.addConstr(
            #     A_i_tilde @ Q @ A_i_tilde.T <= b_i - A_i_tilde @ c
            # )
            m.addConstr(
                A_i_tilde @ P @ P.T @ A_i_tilde.T <= b_i - A_i_tilde @ c
            )
        
        # 2. Containment OF the polytope Theta
        R_Theta = np.hstack((np.zeros((theta_dim,x_dim)),np.eye(theta_dim)))
        for vertex_idx in range(len(Theta_vertices)):
            z_i = z_vertices[vertex_idx]
            v_Theta_i = Theta_vertices[vertex_idx]
            m.addConstr(
                v_Theta_i == R_Theta @ P @ z_i + R_Theta @ c
            )

        # Try To Solve Model
        m.setParam(gp.GRB.Param.NonConvex, 2)
        m.write("logs/find_ellipsoid-test3-nonconvex.rlp")
        start = time.time()
        m.optimize()
        end = time.time()
        print('solving it takes %.3f s'%(end - start))

        self.assertTrue(end-start > 0)
        for row_index in range(dim):
            for col_index in range(dim):
                self.assertTrue(P.X[row_index,col_index] >= 0)

    """
    test_optimization3
    Description:
        Tests the optimization with the nonconvex constraints
        when ran in the function.
    """
    def test_optimization3(self):
        # Constants
        X_interm = pc.Polytope(
            np.vstack((np.eye(2),-np.eye(2))),
            np.array([1.4,1.3,1.2,1.1])
        )
        Theta = pc.Polytope( np.array([[1.0],[-1.0]]) , np.array([0.8,-0.5]) )

        x_dim = X_interm.dim

        # Algorithm

        timing_a, status_a, P_a, c_a = find_ellipsoid_for_intermediate_set(X_interm, Theta)

        print(timing_a)
        print("status_a = ",status_a)
        print(P_a)
        print(c_a)

        self.assertTrue(timing_a > 0)
        dim = X_interm.dim + Theta.dim
        for row_index in range(dim):
            for col_index in range(dim):
                self.assertTrue(P_a[row_index,col_index] >= 0)

if __name__ == '__main__':
    unittest.main()