# %%
from case_studies.A_arm.generate_state_variable_form import *

# %%
############################################################
### Defining vectors for x_dot, x, and u, then taking partial derivatives
############################################################

# defining derivative of states, states, and inputs symbolically
### for this first one, keep in mind that thetadd_eom is actually a full
### row of equations, while thetadd is just the symbolic variable itself.
state_variable_form = Matrix([[thetad], [thetadd_eom]])
states = Matrix([[theta], [thetad]])
inputs = Matrix([[tau]])


# finding the jacobian with respect to states (A) and inputs (B)
A = state_variable_form.jacobian(states)
B = state_variable_form.jacobian(inputs)

# sub in values for equilibrium points (x_e, u_e) or (x_0, u_0)
A_lin = simplify(
    A.subs([(thetad, 0.0), (theta, 0.0), (tau, m * g * ell / 2.0 * cos(0))])
)
B_lin = simplify(
    B.subs([(thetad, 0.0), (theta, 0.0), (tau, m * g * ell / 2.0 * cos(0))])
)

display(Math(vlatex(A_lin)))
display(Math(vlatex(B_lin)))


# %%
# TODO - this form of substitution is OK because they are all 0, however, substituting
# 0 (or any constant) for theta, will also make thetad zero. This is a known issue - (see "substitution"
# at https://docs.sympy.org/latest/modules/physics/mechanics/advanced.html). Future
# examples should substitute for dynamic symbols after the EOM are derived, and before
# substituting in the equilibrium points. It may be OK because velocity variables like
# thetad and zd should be zero at an equilibrium point, and would be for any constant value
# of the generalized coordinates in the state vector, but it is something to be aware of.

# see also use of "msubs" at link above for speed up in substitution in large expressions.
