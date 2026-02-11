# %%
################################################################################
# This file is meant to be run interactively with VSCode's Jupyter extension.
# It works as a regular Python script as well, but the printed results will not
# display as nicely.
################################################################################

# %% [markdown]
# # Lab H.4: Linearize the Dynamics
# ### Load H.3 and H.2
# %%
# It is generally not recommended to import * (everything) from any module or
# package, but in this case we are essentially extending the previous lab file...
# they would normally be in the same file, but we are breaking them up to have
# one file per lab assignment.
from case_studies.H_hummingbird.generate_state_variable_form import *

# This makes it so printing from su only happens when running this file directly
su.enable_printing(__name__ == "__main__")

# %% [markdown]
# # Part 1: Form EOM (LHS and RHS) using F and $\tau$ as control inputs
# #### Convert control inputs from $(f_l, f_r)$ to $(F, \tau)$
# %%
# Previously, we defined generalized forces "tau" in terms of f_l and f_r
# because those are closer to the actual physical inputs (PWM → forces).
# However, for SISO or xfer function-based control design, we split the 
# dynamics into:
#   - Longitudinal (pitch): controlled by total thrust F = f_l + f_r
#   - Lateral (roll/yaw): controlled by body torque τ = d(f_l - f_r)
#
# We'll algebraically solve for f_l and f_r in terms of F and τ using sp.solve
F, T = sp.symbols("F tau")  # T for torque (tau already used for generalized forces)

# %%
# Set up the system of equations: F = f_l + f_r  and  τ = d(f_l - f_r)
# TODO: Set up and solve the system of equations for f_l and f_r
# HINT: Create two equations using sp.Eq(): F = f_l + f_r  and  τ = d(f_l - f_r)
# HINT: Use sp.solve([eq1, eq2], [f_l, f_r]) to solve the system
# HINT: Extract solutions using solution[f_l] and solution[f_r]

# fl_as_FT = None  # TODO: replace with solution for f_l
# fr_as_FT = None  # TODO: replace with solution for f_r

su.printeq("f_l", fl_as_FT)
su.printeq("f_r", fr_as_FT)

tau_FT = tau.subs({f_l: fl_as_FT, f_r: fr_as_FT})
tau_FT = sp.simplify(tau_FT)

# Verify that f_l and f_r are completely eliminated
assert not tau_FT.has(f_l) and not tau_FT.has(f_r), "f_l or f_r still in tau_FT!"
su.printeq("tau_FT", tau_FT)

# %% [markdown]
# #### Form full EOM using new control inputs
# Now that tau is in terms of F and τ, form: $M \ddot{q} + C + \frac{\partial P}{\partial q} = \tau - B \dot{q}$
# %%
# TODO: define qddot, then form LHS and RHS of EOM using tau_FT, and finally create the full EOM
# HINT: qddot can be found by differentiating qdot with respect to time
# HINT: RHS = tau_FT - B_mat @ qdot, where B_mat is the friction matrix (B_mat = sp.eye(3, 3) * beta)

# eom = None  # TODO: replace with full EOM
su.printsym(eom)

# %% [markdown]
# # Part 2: Find equilibrium values
# We need equilibrium $(q_e, F_e, \tau_e)$ such that $\ddot{q} = \dot{q} = 0$
# %% [markdown]
# #### Set all derivatives to zero and replace state/input variables with equilibrium symbols
# %%
# Create equilibrium symbols for states and inputs
phi_eq, theta_eq, psi_eq = sp.symbols("phi_e theta_e psi_e")
F_eq, T_eq = sp.symbols("F_e tau_e")

# Build substitution dictionary: replace angles with equilibrium symbols
# TODO build a dictionary called eq_subs that maps phi, theta, psi to their equilibrium symbols
# HINT: eq_subs = {phi: phi_eq, theta: theta_eq,...}

# TODO add first and second derivatives (all zero at equilibrium) to eq_subs
# HINT: eq_subs[phidot] = 0, eq_subs[phiddot] = 0, etc.

# TODO add control inputs F_e and T_e (replace with equilibrium symbols) to eq_subs
# HINT: eq_subs[F] = F_eq, etc. 

# TODO Apply all substitutions to the EOM
# eom_eq = None  # TODO: replace with EOM after substitutions
su.printsym(eom_eq)

# %% [markdown]
# #### Substitute $\tau_e=0$
# From inspection the 1st row (phi equation), we can see that τ_e must be 0 for equilibrium.
# This makes physical sense: any net body torque would cause rotational motion.
# %%
T_eq_val = 0
eom_eq = eom_eq.subs({T_eq: T_eq_val})
su.printsym(eom_eq)

# %% [markdown]
# #### Substitute $\phi_e=0$
# Looking at the psi equation, one of three must be true: F_e=0, φ_e=0, or θ_e=0
# We need F_e ≠ 0 to balance gravity, so we must choose between φ_e and θ_e.
# Physical reasoning: If φ ≠ 0, the thrust vector has a lateral component causing
# sideways motion. However, a pitched orientation (θ ≠ 0) can be stable.
# Therefore, we set φ_e = 0.

# %%
phi_eq_val = 0
eom_eq = eom_eq.subs({phi_eq: phi_eq_val})
su.printsym(eom_eq)

# %% [markdown]
# #### Substitute $\theta_e=0$
# Technically, θ_e can be any value, but we need a specific equilibrium point.
# We choose θ_e = 0 because it's a natural operating point (level hover).

# %%
theta_eq_val = 0
eom_eq = eom_eq.subs({theta_eq: theta_eq_val})
su.printsym(eom_eq)

# %% [markdown]
# #### Solve for $F_e$
# %%
# TODO Now solve for F_e from the remaining equation (theta equation)
# Hint: Use sp.solve(eom_eq, F_eq) to solve for F_e
F_eq_expr = None  # TODO: replace with solution for F_e

F_eq_expr = sp.factor(F_eq_expr)  # factors out g
su.printeq("F_e", F_eq_expr)

# ψ_e does not appear in the equations, so we choose ψ_e = 0 for simplicity
psi_eq_val = 0

# %% [markdown]
# # Part 3: Longitudinal Dynamics (Pitch Control)
# #### Step 1: Extract longitudinal EOM
# %%
# Use the middle row (theta equation) from the full EOM for longitudinal dynamics
# TODO: extract the middle row (theta equation) from the full EOM and assign it to eom_lon
eom_lon = None  # TODO: replace with theta equation from full EOM
su.printsym(sp.expand_trig(eom_lon))

# %% [markdown]
# #### Step 2: Decouple from lateral dynamics
# %%
# Assume lateral motion is at equilibrium (φ=ψ=0, their derivatives=0)
# TODO perform necessary substitutions to decouple the longitudinal EOM from lateral states and inputs
# HINT: Substitute φ=0, ψ=0, and their derivatives=0 into eom_lon
# HINT: Also substitute F=F_e (equilibrium thrust) to focus on deviations from equilibrium

eom_lon = None  # TODO: replace with eom_lon after substitutions
su.printsym(eom_lon)

# %% [markdown]
# #### Step 3: Feedback linearization
# %%
# Cancel nonlinear terms by choosing F = F_fl + F_ctrl where F_fl cancels gravity.
# The linear terms (with θ̇ and θ̈) remain unchanged.
# TODO perform feedback linearization to cancel nonlinear terms and get a linearized 
# EOM in terms of F_ctrl

# HINT: Set θ̇=0 and θ̈=0 to isolate the nonlinear (gravity) terms, then solve for F_fl 
# that cancels those terms. Finally, substitute F = F_fl + F_ctrl to get the linearized EOM.

F_fl = None  # TODO: replace with solution for F_fl
su.printeq("F_fl", F_fl)

# Apply feedback linearization: F = F_fl + F_ctrl
F_ctrl = sp.symbols("F_ctrl")
feedback_lin_subs = {F: F_fl + F_ctrl}
linearized_eom_lon = eom_lon.subs(feedback_lin_subs)
linearized_eom_lon = sp.simplify(linearized_eom_lon)

# Solve for θ̈ in terms of F_ctrl
# TODO solve for θ̈ in terms of F_ctrl and rearrange to get the linearized EOM
# HINT: Use sp.solve get θ̈

linearized_eom_lon = None  # TODO: replace with linearized EOM in terms of θ̈ and F_ctrl
su.printsym(linearized_eom_lon)

# %% [markdown]
# #### Step 4: Simplify by neglecting friction
# %%
# TODO Set β=0 to ignore rotational damping (typically small for drones)

thetaddot_eom = None  # TODO: replace with simplified θ̈ EOM after setting β=0
su.printeq(thetaddot, thetaddot_eom)

# %% [markdown]
# # Part 4: Lateral Dynamics (Roll and Yaw Control)
# #### Step 1: Extract lateral EOM
# %%
# Use the first (phi) and third (psi) rows from the full EOM
eom_lat = sp.Eq(
    sp.Matrix([eom.lhs[0], eom.lhs[2]]),
    sp.Matrix([eom.rhs[0], eom.rhs[2]])
)
su.printsym(eom_lat)

# %% [markdown]
# #### Step 2: Decouple from longitudinal dynamics
# %%
# TODO: Assume longitudinal motion is at equilibrium: θ=0, θ̇=0, F=F_e
# TODO: Also neglect friction (β=0) for simplicity
# HINT: Substitute the longitudinal equilibrium values for θ, θ̇, F, and β into the lateral EOM

eom_lat = None  # TODO: replace with eom_lat after substitutions
su.printsym(eom_lat)

# %% [markdown]
# #### Step 3: Jacobian linearization
# %%
# Linearize using Jacobian: ẋ = f(x,u) → δẋ = A·δx + B·δu
# State: x = [φ, ψ, φ̇, ψ̇]ᵀ,  Input: u = [τ]ᵀ
# TODO: Solve for φ̈ and ψ̈ from eom_lat, then create xdot_lat, x_lat, u_lat, and compute Jacobians 
# HINT: Use sp.solve to get φ̈ and ψ̈ in terms of states and inputs

# A_lat = None  # TODO: replace with Jacobian df/dx 
# B_lat = None  # TODO: replace with Jacobian df/du 

# substitute equilibrium values into the Jacobians
A_lat = sp.cancel(A_lat.subs(lat_eq_vals))  # cancel cleans up signs
B_lat = B_lat.subs(lat_eq_vals)

# Write linearized EOM using tilde notation (deviations from equilibrium)
# We focus on the acceleration equations (last 2 rows)

x_tilde_dot_lat = sp.Matrix(
    [
        su.TildeSymbol("phi", 1),
        su.TildeSymbol("psi", 1),
        su.TildeSymbol("phi", 2),
        su.TildeSymbol("psi", 2),
    ]
)
x_tilde_lat = sp.Matrix(
    [
        su.TildeSymbol("phi"),
        su.TildeSymbol("psi"),
        su.TildeSymbol("phi", 1),
        su.TildeSymbol("psi", 1),
    ]
)
u_tilde_lat = sp.Matrix([su.TildeSymbol("tau")])

linearized_eom_lat = sp.Eq(
    x_tilde_dot_lat[2:, :], A_lat[2:, :] @ x_tilde_lat + B_lat[2:, :] @ u_tilde_lat
)
su.printsym(linearized_eom_lat)

# %%
