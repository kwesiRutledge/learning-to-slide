"""
aclf.py
Description:
    This file contains functions and constants related to our new adaptive clf.
"""

import time, copy

import polytope as pc
import gurobipy as gp
import pypoman as ppm

import numpy as np

"""
MVar_determinant
Description:
    Computes the determinant of the gurobipy matrix variable mvar_in.
"""
def MVar_determinant(A:list[list[gp.Var]]) -> gp.Var:
    # Constants
    dim = len(A)
    indices = list(range(dim))

    # Input Checking
    if dim > 3:
        raise( Exception("This determinant function is only defined for a 3-dimensional variable") )

    # Base Case of Recursive Algorithm
    if dim == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val

    # Recursive calls (define submatrix for focus column and call this function)
    total = 0
    for fc in indices:
        As = A.copy()
        As = As[1:]
        height = len(As)

        for i in range(height):
            As[i] = As[i][0:fc] + As[i][fc+1:]

        print(As)
        sign = (-1) ** (fc % 2)
        sub_det = MVar_determinant(As)
        total += sign * A[0][fc]*sub_det

    return total


"""
find_ellipsoid_for_intermediate_set
Description:

Inputs:
    Theta = Polytope that defines the set of feasible parameters for the system.
"""
def find_ellipsoid_for_intermediate_set( X_interm:pc.Polytope, Theta:pc.Polytope ):
    # Constants
    x_dim = X_interm.dim
    theta_dim = Theta.dim
    dim = x_dim + theta_dim

    Theta_vertices = ppm.compute_polytope_vertices(Theta.A,Theta.b)

    # Algorithm

    m = gp.Model()

    # Create Optimizaiton Variables
    P = m.addMVar((dim,dim))
    c = m.addMVar((dim,1))

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

    # 3. Create Log Det expression for objective
    det = m.addVar(lb=-gp.GRB.INFINITY,vtype=gp.GRB.CONTINUOUS,name="determinant")
    logDet = m.addVar(lb=-gp.GRB.INFINITY,vtype=gp.GRB.CONTINUOUS,name="logDet")

    print(P)

    if dim != 3:
        raise( Exception("This function is currently designed to compute the three-dimensional determinant, but your problem has dimension " + str(dim) + ".") )

    m.addConstr(det == MVar_determinant(P.tolist()))
    m.addGenConstrLog(det,logDet,name="logDetConstr")

    # Create Objective
    m.setObjective(logDet,gp.GRB.MAXIMIZE)

    # Try To Solve Model
    m.setParam(gp.GRB.Param.NonConvex, 2)
    m.write("find_ellipsoid.rlp")
    start = time.time()
    m.optimize()
    end = time.time()
    print('solving it takes %.3f s'%(end - start))

    # Return Values if status is good
    P_opt_val = np.zeros((dim,dim))
    c_opt_val = np.zeros((dim,1))
    if (m.Status == 2): #2 means the problem is solved and optimal.
        P_opt_val = P.X
        c_opt_val = c.X

    return end-start, m.Status, P_opt_val, c_opt_val



