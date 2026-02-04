# 3rd-party
import numpy as np

# local (controlbook)
from . import params as P
from ..common.dynamics_base import DynamicsBase
from . import eom_generated


class HummingbirdDynamics(DynamicsBase):
    def __init__(self, alpha=0.0):
        super().__init__(
            # Initial state conditions
            state0=np.array(
                [P.phi0, P.theta0, P.psi0, P.phidot0, P.thetadot0, P.psidot0]
            ),
            # Input limits [f_l, f_r]
            u_max=np.inf,  # we will saturate later based on input limits
            u_min=0.0,  # no negative forces
            # Time step for integration
            dt=P.ts,
        )
        # See params.py/textbook for details on these parameters

        # Parameter randomization is introduced in Chapter 10, you do not need
        # to include it before then
        # NOTE: Keys must match the parameter names in eom_generated.py functions
        self.eom_params = {
            "m1": self.randomize_parameter(P.m1, alpha),
            "m2": self.randomize_parameter(P.m2, alpha),
            "m3": self.randomize_parameter(P.m3, alpha),
            "J_1x": self.randomize_parameter(P.J1x, alpha),
            "J_1y": self.randomize_parameter(P.J1y, alpha),
            "J_1z": self.randomize_parameter(P.J1z, alpha),
            "J_2x": self.randomize_parameter(P.J2x, alpha),
            "J_2y": self.randomize_parameter(P.J2y, alpha),
            "J_2z": self.randomize_parameter(P.J2z, alpha),
            "J_3x": self.randomize_parameter(P.J3x, alpha),
            "J_3y": self.randomize_parameter(P.J3y, alpha),
            "J_3z": self.randomize_parameter(P.J3z, alpha),
            "ell_1": self.randomize_parameter(P.ell1, alpha),
            "ell_2": self.randomize_parameter(P.ell2, alpha),
            "ell_3x": self.randomize_parameter(P.ell3x, alpha),
            "ell_3y": self.randomize_parameter(P.ell3y, alpha),
            "ell_3z": self.randomize_parameter(P.ell3z, alpha),
            "ell_T": self.randomize_parameter(P.ellT, alpha),
            "d": self.randomize_parameter(P.d, alpha),
            "g": P.g,
        }

        p = self.eom_params
        self.km = P.g * (p["m1"] * p["ell_1"] + p["m2"] * p["ell_2"]) / p["ell_T"]
        self.B = P.beta * np.eye(3)  # friction-based terms

    ############################################################################
    # The following 4 methods are where we use the generated functions from
    # SymPy to calculate the terms in the equations of motion: the mass matrix
    # M, the coriolis vector C, the gravity terms dP_dq, and the generalized
    # force vector tau. Each of these methods are just wrappers to call the
    # generated functions with the appropriate parameters.
    # See the end of h3_generate_EL.py for an example of how to call the
    # generated functions.

    # NOTE: IF YOU DID NOT GENERATE THESE FUNCTIONS, YOU WOULD NEED TO
    # IMPLEMENT THE CALCULATION OF THESE TERMS BY HAND.

    # IF YOU DID GENERATE THEM, JUST CALL THE FUNCTIONS HERE BELOW. This would
    # look like the following:
        # eom_generated.calculate_Y(x, **self.eom_params) OR
        # eom_generated.calculate_Y(x, u, **self.eom_params)
    # SEE eom_generated.py FOR THE EXACT FUNCTION NAMES AND ARGUMENTS.


    def calculate_M(self, x):
        # M = 
        return M

    def calculate_C(self, x):
        # C
        return C

    def calculate_dP_dq(self, x):
        # dP_dq = 
        return dP_dq

    def calculate_tau(self, x, u):
        # tau = 
        return tau

    ############################################################################

    def f(self, x, u):
        """
        Full nonlinear dynamic model `xdot = f(x,u)`.

        Args:
            x (1D numpy array): state vector [q, qdot] or
                [phi, theta, psi, phidot, thetadot, psidot]
            u (1D numpy array): input vector [f_l, f_r] (later will change to PWM values)
        Returns:
            xdot (1D numpy array): time derivative of state vector
                [phidot, thetadot, psidot, phiddot, thetaddot, psiddot]
        """
        # Ensure inputs are flattened for consistent indexing
        x = x.flatten()
        u = u.flatten()
        
        # Re-label terms for readability
        qdot = x[3:6]

        #TODO: Find qddot using the equations of motion, then formulate and return xdot
        # qddot =

        xdot = np.concatenate((qdot, qddot))
        return xdot

    def h(self):
        """
        Output function `y = h(x)`, where x is self.state and does not need to
        be passed in. Here we assume we can measure all the generalized
        coordinates, i.e., position states.

        Returns:
            y (1D numpy array): output vector [phi, theta, psi]
        """
        #TODO: return the measured outputs based on self.state - these would be the first three states
        # y = 
        return y

    def update(self, pwm):
        """
        Integrates the system dynamics forward one time step using RK4.
        This version accepts PWM inputs and converts them to forces.

        Args:
            pwm (NDArray[np.float64]): input vector [pwm_l, pwm_r].
        Returns:
            y (NDArray[np.float64]): measured output vector [phi, theta, psi].
        """
        # saturate pwm
        pwm = np.clip(pwm, 0.0, 1.0)

        # convert pwm to motor forces [f_l, f_r]
        u = pwm * self.km

        # call update from parent class (DynamicsBase)
        y = super().update(u)
        return y
