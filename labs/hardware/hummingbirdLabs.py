import sys
import serial

import arduinoCommunication as ardCom
import hummingbirdParam as P
import tkinter as tk
import numpy as np

from hummingbirdAnimation import HummingbirdAnimation
from threading import Thread
from HB_GUI import HB_GUI

# ----------------GUI Options--------------------
P.H_number = 10  # Labs 4, 7, 8, 10
animate = False
P.controller = P.controller_pitch if P.H_number == 7 else P.controller_km

# ----------------Check for Serial Connection--------------------
try:
    # Change the com port
    ser = serial.Serial("/dev/ttyUSB0", 115200)
except Exception as e:
    print()
    print("Hummingbird not connected")
    print("Check COM Port")
    print()
    print(e)
    print()
    sys.exit()

# ----------------Start Arduino Thread--------------------
thread = Thread(target=lambda: ardCom.read_from_arduino(ser))
thread.daemon = True  # This thread dies when the main thread (program) exits
thread.start()

# ----------------Initialize the Root Event Handler--------------------
root = tk.Tk()
root.title("Controller Interface")
GUI_Frames = P.tuning_parameters[P.H_number]
GUILayout = HB_GUI(root, GUI_Frames)

# ----------------Initialize Animation--------------------
if animate:
    animation = HummingbirdAnimation()

    def update_animation():
        global t
        t = t + P.Ts_Animation
        t = np.mod(t, 29)
        animation.update(t, np.array([[P.phi], [P.theta], [P.psi]]))
        root.after(int(P.Ts_Animation * 1000), update_animation)

    update_animation()


# ----------------GUI Communication--------------------
def update_GUI():
    GUILayout.GUI_Values(P.phi, P.theta, P.psi)
    ardCom.send_gains(ser)
    root.after(int(P.Ts_GUI * 1000), update_GUI)


update_GUI()


def flip_phi():
    GUILayout.Flip_Phi()
    root.after(P.t_flip * 1000, flip_phi)


def flip_theta():
    GUILayout.Flip_Theta()
    root.after(P.t_flip * 1000, flip_theta)


def flip_psi():
    GUILayout.Flip_Psi()
    root.after(P.t_flip * 1000, flip_psi)


if P.Flip_Roll:
    flip_phi()

if P.Flip_Pitch:
    flip_theta()

if P.Flip_Yaw:
    flip_psi()

# ----------------Start--------------------
root.mainloop()

# ----------------Close--------------------
ser.close()
