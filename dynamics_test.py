import numpy as np


# Define New Class fo Dynamical System
"""
SlidingBlock
Description:
    This describes the dynamics of the sliding block system which describes the motion of a block sliding on the ground
    when it is being pushed (i.e. force being applied by a manipulator) and has friction with the ground.
Notes:
    Dynamics structure

        x-dot(t) = f(x) + F(x) * theta + ( g(x) + G(x) * theta ) * u

    Right now, I assume that the block is only moving in the positive directions (v_x > 0)
"""
class SlidingBlock():
    def __init__(self,mass_in=1.0,mu_kinetic_in=0.25):
        self.m = mass_in # Mass in kg of block
        self.mu_k = mu_kinetic_in # Coefficient of kinetic friction on the block

        self.n_x = 2 # Dimension of the state space (x)
        self.n_u = 2 # Dimensions of the input space (u)
    
    """
    f(x)
    Description:
        Drift term of the dynamics that is INDEPENDENT of unknown parameter.
    """
    def f(self,x):
        if len(x) != self.n_x:
            raise(Exception("State input to f() has dimension " + str(len(x)) + " expected 2." ))

        return np.array([[0,1.0],[0,0]]).dot(x)

    """
    F(x)
    Description:
        Drift term that is multiplied by the unknown parameter before being added to the other terms.
    """
    def F(self,x):
        # Error Handling
        if len(x) != self.n_x:
            raise(Exception("State input to F() should have dimension 2, instead received vector of length" + str(len(x)) + "." ))

        # Create constant
        g = 9.81 # For gravity
        return np.array([0.0,-g])

    """
    g(x)
    Description:
        Constant term that is multiplied by the input.
    """
    def g(self,x):
        # Error Handling
        if len(x) != self.n_x:
            raise(Exception("State input to g() should have dimension 2, instead received vector of length" + str(len(x)) + "." ))

        # Get constants
        m = self.m
        return np.array([[0.0,0.0],[1/m,0.0]])

    """
    G(x)
    Description:
        Term that is multiplied by the parameter to create a term that is multiplied by input.
        This term is added to the dynamics.
    """
    def G(self,x):
        # Error Handling
        if len(x) != self.n_x:
            raise(Exception("State input to g() should have dimension 2, instead received vector of length" + str(len(x)) + "." ))

        # Get constants
        m = self.m
        mu_k = self.mu_k
        return np.array([[0.0,0.0],[0.0,(mu_k/m)]])


    """
    dynamics(x,u)
    Description:
        Term that produces the derivative of the state as a function of the current state and input (the )
    """
    def dynamics(self,x,u):
        # Error Handling
        if len(x) != self.n_x:
            raise(Exception("State input to dynamics() should have dimension 2, instead received vector of length" + str(len(x)) + "." ))
        if len(u) != self.n_u:
            raise(Exception("Input input to dynamics() should have dimension 2, instead received vector of length" + str(len(u)) + "." ))

        # Compute the derivative using the previous components
        eps0 = 1e-4
        mu_k = self.mu_k
        if x[1] > eps0:
            return self.f(x) + self.F(x)*mu_k + ( self.g(x) + self.G(x)*mu_k ).dot(u)
        elif (-eps0 <= x[1]) and (x[1] <= eps0):
            return np.array([0.0,0.0])
        else:
            x_dot = self.f(x) + self.F(x)*mu_k + ( self.g(x) + self.G(x)*mu_k ).dot(u)
            x_dot[1] = 0.0
            return x_dot


# Create a simple example to show how to compute the derivatives of everything.
print("Testing derivative at multiple places where velocity and inputs are zero. (If both are zero, then block should not move).")

tempblock = SlidingBlock()

x0 = np.array([0.5,0.0])
u0 = np.array([0.0,0.0])

print("dynamics(x0,u0) = " + str(tempblock.dynamics(x0,u0)))

x1 = np.array([0.0,0.0])
u1 = np.array([0.0,-0.50])
print("dynamics(x1,u1) = " + str(tempblock.dynamics(x1,u1)))

print(" ")

print("What I'm realizing is that my model is ignoring a lot of frictional characterstics. We are not assuming that it works starting from rest, which is kind of interesting.")
print(" ")

# Testing simulating individual trajectories
from scipy.integrate import odeint
import matplotlib.pyplot as plt

print("Testing odeint")

def f0(y,t,params):
    # Create Constants
    mu_k = params[0]
    u0 = np.zeros((2,))
    
    block1 = SlidingBlock(1.0,mu_k)
    
    return block1.dynamics(y,u0)

# Integrate using odeint
x2 = [0.5,3.5]
mu_k = 0.25
params = [mu_k]

tstop = 10
deltat = 0.1
t = np.arange(0.0,tstop,deltat)

psoln = odeint(f0,x2,t,args=(params,))

#print("odeint's solution:")
#print(psoln)

plt.figure()
plt.plot(psoln)
plt.legend(
    ['x0','x1']
)
plt.title("Trajectory of the sliding block (intial state = x2)")

plt.savefig("images/dynamics_test1.png")