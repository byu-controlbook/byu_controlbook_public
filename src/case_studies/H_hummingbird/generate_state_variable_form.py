# %%
################################################################################
# This file is meant to be run interactively with VSCode's Jupyter extension.
# It works as a regular Python script as well, but the printed results will not
# display as nicely.
################################################################################

# %% [markdown]
# # Lab H.3: Equations of Motion
# ### Load H.2

# %%
# It is generally not recommended to import * (everything) from any module or
# package, but in this case we are essentially extending the previous lab file...
# they would normally be in the same file, but we are breaking them up to have
# one file per lab assignment.
from case_studies.H_hummingbird.generate_KE import *
from case_studies.common import sym_utils as su

# This makes it so printing from su only happens when running this file directly
su.enable_printing(__name__ == "__main__")

# %%[markdown]
# # Part 1: Calculate Kinetic Energy, Potential Energy

# %%
# TODO calculate the kinetic energy using the definition with the mass matrix "M"
# and "qdot"
# K =

# then make sure to grab just the scalar part of the result
K = K[0, 0]


# TODO: define symbols needed for potential energy and the RHS 
# (e.g. friction, force_left, force_right, etc.) of the equations 
# of motion.


# TODO now calculate the potential energy "P"
# P = 

su.printeq("P", P)



# %% [markdown]
# # Part 2: Calculate C, dP/dq, and tau
# ### Generalized Forces (tau):
# %%
# TODO: now calculate and define tau
# tau = sp.Matrix([])


# Group terms together for readability (this is to help with checking the result
# but is not strictly necessary)
tau[-1] = sp.collect(tau[-1], [sp.sin(theta), sp.sin(phi) * sp.cos(theta)])
tau[-1] = sp.collect(tau[-1], [f_l + f_r, f_l - f_r])

su.printeq("tau", tau)


# %% [markdown]
# ### Coriolis Matrix (C):
# The lab manual does not show what C is directly, but rather the intermediate
# terms of Mdot, dM/dphi, dM/dtheta, and dM/dpsi.
# %%
# TODO: calculate Mdot and verify with lab manual
# Mdot = 


#%%[markdown]
# The following symbolic manipulations are just to make Mdot match the lab manual
# more closely. They are not strictly necessary.

# %%
# Print Mdot with some terms substituted to make it easier to read
Mdot22 = Mdot[1, 1]
Mdot23 = Mdot[1, 2]
Mdot33 = Mdot[2, 2]
print_subs = {
    Mdot22: sp.Symbol("Mdot_22"),
    Mdot33: sp.Symbol("Mdot_33"),
    Mdot23: sp.Symbol("Mdot_23"),
}
su.printeq("Mdot", Mdot.subs(print_subs))

# Note: using Python 3.13.7 and SymPy 1.14.0 with these simplifications matched
# the lab manual results pretty well, but other versions may look different

# Wrap Mdot inside of sp.Matrix call so that individual elements can be modified
Mdot = sp.Matrix(Mdot)

Mdot22 = sp.collect(Mdot22, 2 * sp.sin(phi) * sp.cos(phi) * qdot[0])
su.printeq("Mdot_22", Mdot22)
Mdot[1, 1] = Mdot22  # replace with simplified version

Mdot23 = sp.trigsimp(sp.factor(Mdot23))
su.printeq("Mdot_23", Mdot23)
Mdot[1, 2] = Mdot23  # replace with simplified version
Mdot[2, 1] = Mdot23  # Mdot is symmetric

Mdot33 = sp.factor(Mdot33)
Mdot33 = sp.collect(Mdot33, [qdot[0], qdot[1] * sp.sin(theta)])
Mdot33 = sp.expand_trig(sp.trigsimp(sp.collect(Mdot33, Mdot33.free_symbols)))
su.printeq("Mdot_33", Mdot33)
Mdot[2, 2] = Mdot33  # replace with simplified version



# %%
# TODO: calculate the partial derivatives of M with respect to each generalized
# coordinate and verify with lab manual (if they do not match perfectly, you can
# try using sp.trigsimp(), sp.simplify(), sp.collect(), etc., but you may want
# to just move on and use numerical comparisons in the testDynamics.py file to
# verify your matrices are correct)

# Note: using Python 3.13.7 and SymPy 1.14.0 with these simplifications matched
# the lab manual results pretty well, but other versions may look different

# %%

#dM_dphi = 
#dM_dphi = sp.simplify(dM_dphi)

su.printeq("dM/dϕ", dM_dphi)

# %%
# dM_dtheta =

# substitute N33 just for printing to match the lab manual
N33 = dM_dtheta[2, 2]
print_subs = {N33: sp.Symbol("N_33")}
su.printeq("dM/dθ", dM_dtheta.subs(print_subs))

N33 = sp.trigsimp(sp.collect(N33, N33.free_symbols))
N33 = sp.trigsimp(sp.collect(N33, N33.free_symbols))  # repeat simplifies further
su.printeq("N_33", N33)

# wrap dM_dtheta inside of sp.Matrix call so that individual elements can be modified
dM_dtheta = sp.Matrix(dM_dtheta)
dM_dtheta[2, 2] = N33

# %%
# dM_dpsi =

su.printeq("dM/dψ", dM_dpsi)


# %% [markdown]
# ### Coriolis Vector (C):
# %%
# TODO: calculate C
# intermediate terms may be helpful, but calculate C -> 
# C = 


# can print C directly, but it is very long
# su.printeq("C", C)

# %% [markdown]
# ### Partial Derivative of Potential Energy (dP/dq):
# %%
# dP_dq = 

su.printeq("dP/dq", dP_dq)

# %% [markdown]
# # Step 3: Create functions for M, C, dP/dq, and tau
# %%
# TODO: now you can either run the code below to generate functions to calculate
# M, C, tau, and dP_dq (where you can import the generated functions into your 
# dynamics class), or you can hand-code the results you have found above into the 
# hummingbirdDynamics.py file directly.

# Define state and input vectors for function generation
phidot, thetadot, psidot = qdot
state = sp.Matrix([phi, theta, psi, phidot, thetadot, psidot])
u = sp.Matrix([f_l, f_r])

# Create a list of all parameters needed for the functions
# TODO: you may need to modify this to match the variable names you created above
sys_params = list(
    [m1, m2, m3, J1x, J1y, J1z, J2x, J2y, J2z, J3x, J3y, J3z]
    + [ell_1, ell_2, ell_3x, ell_3y, ell_3z, ell_T, d, g]
)

# Everything below this point will only run if this file is executed directly,
# not if it is imported as a module.
if __name__ == "__main__":
    from case_studies import H_hummingbird

    su.write_eom_to_file(
        state, u, sys_params, H_hummingbird, M=M, C=C, dP_dq=dP_dq, tau=tau
    )

    # %% [markdown]
    # # Part 4: (Optional) Test your generated EOM functions
    # %%
    # You can use your generated EOM functions as follows:
    import numpy as np
    from case_studies.H_hummingbird import eom_generated
    import importlib

    importlib.reload(eom_generated)  # reload in case it was just generated/modified

    P = H_hummingbird.params

    param_vals = {
        "m_1": P.m1,
        "m_2": P.m2,
        "m_3": P.m3,
        "J_1x": P.J1x,
        "J_1y": P.J1y,
        "J_1z": P.J1z,
        "J_2x": P.J2x,
        "J_2y": P.J2y,
        "J_2z": P.J2z,
        "J_3x": P.J3x,
        "J_3y": P.J3y,
        "J_3z": P.J3z,
        "ell_1": P.ell1,
        "ell_2": P.ell2,
        "ell_3x": P.ell3x,
        "ell_3y": P.ell3y,
        "ell_3z": P.ell3z,
        "ell_T": P.ellT,
        "d": P.d,
        "g": P.g,
    }

    x_test = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    u_test = np.array([0.24468 / 2, 0.22468 / 2])

    M_val = eom_generated.calculate_M(x_test, **param_vals)
    print("\nM_val:")
    print(M_val, "\n")  # should be np.diag([0.0002, 0.0126, 0.0127])

    C_val = eom_generated.calculate_C(x_test, **param_vals)
    print("C_val:")
    print(C_val, "\n")  # should be all zeros

    dP_dq_val = eom_generated.calculate_dP_dq(x_test, **param_vals)
    print("dP_dq_val:")
    print(dP_dq_val, "\n")  # should be [[0], [0.0833], [0]]

    tau_val = eom_generated.calculate_tau(x_test, u_test, **param_vals)
    print("tau_val:")
    print(tau_val, "\n")  # should be [[0.0012], [0.0833], [0]]



    