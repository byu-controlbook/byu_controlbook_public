# %%
# local (controlbook)
from case_studies.common import sym_utils as su

# local imports (from this folder)
from case_studies.B_pendulum.generate_KE import *

# This makes it so printing from su only happens when running this file directly
su.enable_printing(__name__ == "__main__")

# %%[markdown]
# The code imported from above shows how we defined q, q_dot, and necessary system parameters. Then we used position, velocity, and angular velocity to calculate kinetic energy.

# %%
# defining potential energy
P = (
    m1 * g * ell / 2.0 * (cos(theta) - 1)
)  # this is "mgh", where "h" is a function of generalized coordinate "q"

# can also do the following to get the same answer
#   g_vec = Matrix([[0], [g], [0]])  # defining gravity in the direction that increases potential energy
#   P = m1*g_vec.T@p1
#   P = P[0,0]


# calculate the lagrangian, using simplify intermittently can help the equations to be
# simpler, there are also options for factoring and grouping if you look at the sympy
# documentation.
L = simplify(K - P)

# %%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like friction and tau)
EL_case_studyB = simplify(diff(diff(L, qdot), t) - diff(L, q))
# EL_case_studyB = simplify(L.diff(q_dot).diff(t) - L.diff(q))

display(Math(vlatex(EL_case_studyB)))


# %%
############################################################
### Including friction and generalized forces, then solving for highest order derivatives
############################################################

# these are just convenience variables
zd = z.diff(t)
zdd = zd.diff(t)
thetad = theta.diff(t)
thetadd = thetad.diff(t)

# defining symbols for external force and friction
F, b = symbols("F, b")

# defining the right-hand side of the equation and combining it with E-L part
RHS = Matrix([[F - b * zd], [0]])
full_eom = EL_case_studyB - RHS

# finding and assigning zdd and thetadd
# if our eom were more complicated, we could rearrange, solve for the mass matrix, and invert it to move it to the other side and find qdd and thetadd
result = simplify(sp.solve(full_eom, (zdd, thetadd)))

# TODO - add an example of finding the same thing, but not using sp.solve


# result is a Python dictionary, we get to the entries we are interested in
# by using the name of the variable that we were solving for
zdd_eom = result[zdd]  # EOM for zdd, as a function of states and inputs
thetadd_eom = result[thetadd]  # EOM for thetadd, as a function of states and inputs

display(Math(vlatex(zdd_eom)))
display(Math(vlatex(thetadd_eom)))


# %% [markdown]
# OK, now we can get the state variable form of the equations of motion.

# %%
import case_studies.B_pendulum.params as P
import numpy as np

# defining fixed parameters that are not states or inputs (like g, ell, m1, m2, b)
# can be done like follows:
# params = [(m1, P.m1), (m2, P.m2), (ell, P.ell), (g, P.g), (b, P.b)]

# but in this example, I want to keep the masses, length, and damping as variables so
# that I can simulate uncertainty in those parameters in real life.
params = [(g, P.g)]


# substituting parameters into the equations of motion
zdd_eom = zdd_eom.subs(params)
thetadd_eom = thetadd_eom.subs(params)

# now defining the state variables that will be passed into f(x,u)
state = sp.Matrix([z, theta, zd, thetad])
ctrl_input = sp.Matrix([F])

# defining the function that will be called to get the derivatives of the states
state_dot = sp.Matrix([zd, thetad, zdd_eom, thetadd_eom])



# %% [markdown]
# The next step is to save this function "f" so that we can use it with a numerical integrator, like
# scipy.integrate.ivp.solve_ivp or the rk4 functions in the case studies. To save this function, we do the
# following:

# %%
if __name__ == "__main__":
    from case_studies import B_pendulum

    # make sure printing only happens when running this file directly
    su.enable_printing(__name__ == "__main__")

    su.write_eom_to_file(state, ctrl_input, [m1, m2, ell, b], B_pendulum, eom=state_dot)

    import numpy as np
    from case_studies.B_pendulum import eom_generated
    import importlib

    importlib.reload(eom_generated)  # reload in case it was just generated/modified
    P = B_pendulum.params

    param_vals = {
        "m1": P.m1,
        "m2": P.m2,
        "ell": P.ell,
        "b": P.b,
    }

    x_test = np.array([0.0, 0.0, 0.0, 0.0])
    u_test = np.array([1.0])

    x_dot_test = eom_generated.calculate_eom(x_test, u_test, **param_vals)
    print("\nx_dot_test from generated function = ", x_dot_test)
    # should match what was printed earlier when we called eom directly
