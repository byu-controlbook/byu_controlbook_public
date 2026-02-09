# Hummingbird Parameter File
import numpy as np

# Slider Starting Values:
H_number = 4

# km value for simulation is self-consistent and does not need to be identified directly.
km = 0.367

t_flip = 10  # Square wave Duration

phi_ref = 0.0
phi_kp = 0.0
phi_kd = 0.0

theta_ref = 0.0
theta_kp = 0.0
theta_kd = 0.0
theta_ki = 0.0

psi_ref = 0.0
psi_kp = 0.0
psi_kd = 0.0
psi_ki = 0.0

# Update Parameters
Ts_Animation = 0.05  # Refresh rate of the animation
Ts_GUI = 0.1  # Refresh rate of the GUI (Updating encoder values, sending gains)

# Characters for communication
char_pause = "P"  # Pause
char_play = "L"  # Play
char_zero = "Z"  # Zero
char_km = "K"  # Tune Km, Lab H4

char_pitch = "T"  # Tune Pitch (Theta), Lab H7
char_roll = "F"  # Tune Roll (Phi), Lab H8 part 1
char_yaw = "S"  # Tune Yaw (Psi), Lab H8 part 2

char_full = "O"  # Tune Roll, Pitch, Yaw, Lab H9

# Characters to indicate order for sending
char_val1 = "A"
char_val2 = "B"
char_val3 = "C"
char_val4 = "D"
char_val5 = "E"

char_val6 = "G"
char_val7 = "H"
char_val8 = "I"
char_val9 = "J"

char_val10 = "M"
char_val11 = "N"

# Communication Variables
specialByte = 124
startMarker = 125
endMarker = 126

Zero = False
Pause = False
Play = False
reading = False
bytesRec = 0
tempBuffer = ""
GUI_msg_d1 = ""
GUI_msg_d2 = ""

# Variable for Storing GUI Encoder Values
phi = 0.0
theta = 0.0
psi = 0.0

# Variables for H4, Tuning Km
km_min = 0.2
km_max = 0.5
km_step = 0.001

# Variables for H7, Tuning Pitch (Theta)
theta_ref_min = -30
theta_ref_max = 30
theta_ref_step = 1

theta_kp_min = 0
theta_kp_max = 2
theta_kp_step = 0.01

theta_ki_min = 0
theta_ki_max = 1
theta_ki_step = 0.01

theta_kd_min = 0
theta_kd_max = 1
theta_kd_step = 0.01

# Variables for H8 - Part 1
phi_ref_min = -30
phi_ref_max = 30
phi_ref_step = 1

phi_kp_min = 0
phi_kp_max = 0.2
phi_kp_step = 0.0005

phi_kd_min = 0
phi_kd_max = 0.005
phi_kd_step = 0.0005

# Variables for H8 - Part 2
psi_ref_min = -45
psi_ref_max = 45
psi_ref_step = 1

psi_kp_min = 0
psi_kp_max = 1
psi_kp_step = 0.01

psi_ki_min = 0
psi_ki_max = 0.5
psi_ki_step = 0.01

psi_kd_min = 0
psi_kd_max = 0.2
psi_kd_step = 0.01


# Parameters for GUI Generation
def find_grid_location(Lab_Layout, Frame):
    for i, sublist in enumerate(Lab_Layout):
        # Check if the element is in the sublist
        if Frame in sublist:
            # Find the column index of the element in the sublist
            j = sublist.index(Frame)
            # Store both row and column indexes in location
            location = (i, j)
            return location
    return (-1, -1)


controller_km = "Km"
controller_roll = "Roll"
controller_pitch = "Pitch"
controller_yaw = "Yaw"
controller_full = "Full"

controller = controller_km

Frame_Km = "Km"
Frame_Encoder = "Encoder"

Frame_Pitch = "Pitch"
Flip_Pitch = False

Frame_Roll = "Roll"
Flip_Roll = False

Frame_Action = "Action"

Frame_Yaw = "Yaw"
Flip_Yaw = False

Frame_Control = "Control"

Frame_Action_Full = "Full"

Frame_Empty = "Empty"

Tune_km = [[Frame_Km], [Frame_Encoder]]  # Lab H4

Tune_Pitch = [[Frame_Km], [Frame_Pitch], [Frame_Encoder]]  # Lab H7

Tune_Lat = [
    [Frame_Km, Frame_Action],
    [Frame_Roll, Frame_Yaw],
    [Frame_Encoder, Frame_Control],
]  # Lab H8

Tune_Full = [
    [Frame_Empty, Frame_Km, Frame_Action_Full],
    [Frame_Roll, Frame_Pitch, Frame_Yaw],
    [Frame_Empty, Frame_Encoder, Frame_Control],
]  # Lab H10

# Define the dictionary
tuning_parameters = {
    4: Tune_km,
    7: Tune_Pitch,
    8: Tune_Lat,
    10: Tune_Full,
}
