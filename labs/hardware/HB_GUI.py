import hummingbirdParam as P
import tkinter as tk
import numpy as np


class HB_GUI:
    def __init__(self, master, GUI_Frames):
        self.enabled_color = "#e8e8e8"  # White
        self.disabled_color = "#8e8e8e"  # Grey
        self.w_button_2 = 15
        self.w_button_3 = 9
        self.w_button_4 = 6
        self.h_button = 2

        Location_Km = P.find_grid_location(GUI_Frames, P.Frame_Km)
        Location_Roll = P.find_grid_location(GUI_Frames, P.Frame_Roll)
        Location_Encoder = P.find_grid_location(GUI_Frames, P.Frame_Encoder)
        Location_Yaw = P.find_grid_location(GUI_Frames, P.Frame_Yaw)
        Location_Action = P.find_grid_location(GUI_Frames, P.Frame_Action)
        Location_Action_Full = P.find_grid_location(GUI_Frames, P.Frame_Action_Full)
        Location_Pitch = P.find_grid_location(GUI_Frames, P.Frame_Pitch)
        Location_Control = P.find_grid_location(GUI_Frames, P.Frame_Control)

        self.Enabled_Km = False
        self.Enabled_Encoder = False
        self.Enabled_Roll = False
        self.Enabled_Pitch = False
        self.Enabled_Yaw = False

        if Location_Km != (-1, -1):
            self.pack_Km(master, Location_Km)
            self.Enabled_Km = True

        if Location_Encoder != (-1, -1):
            self.pack_Encoder(master, Location_Encoder)
            self.Enabled_Encoder = True

        if Location_Pitch != (-1, -1):
            self.pack_Pitch(master, Location_Pitch)
            P.Flip_Pitch = True
            self.Enabled_Pitch = True

        if Location_Roll != (-1, -1):
            self.pack_Roll(master, Location_Roll)
            P.Flip_Roll = True
            self.Enabled_Roll = True

        if Location_Yaw != (-1, -1):
            self.pack_Yaw(master, Location_Yaw)
            P.Flip_Yaw = True
            self.Enabled_Yaw = True

        if Location_Action != (-1, -1):
            self.set_roll()
            self.pack_Action(master, Location_Action)

        if Location_Control != (-1, -1):
            self.pack_Control(master, Location_Control)

        if Location_Action_Full != (-1, -1):
            self.pack_Action_Full(master, Location_Action_Full)
            self.phi_ref_slider["state"] = "disabled"
            self.phi_ref_slider["troughcolor"] = self.disabled_color

    # ----------------Km Slider Functions--------------------
    def pack_Km(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        km_frame = tk.LabelFrame(master, text="km", padx=10, pady=10)
        km_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        # Create and pack a label for the slider
        km_label = tk.Label(km_frame, text="km:")
        km_label.pack(side=tk.TOP)

        # Create a "Minus" button to decrease the slider value
        km_minus_button = tk.Button(km_frame, text="-", command=self.km_decrease)
        km_minus_button.pack(side=tk.LEFT)

        # Create a slider
        self.km_slider = tk.Scale(
            km_frame,
            length=200,
            from_=P.km_min,
            to=P.km_max,
            resolution=P.km_step,
            orient="horizontal",
        )
        self.km_slider.set(P.km)
        self.km_slider.pack(side=tk.LEFT)
        self.km_slider["troughcolor"] = self.enabled_color

        # Create a "Plus" button to increase the slider value
        km_plus_button = tk.Button(km_frame, text="+", command=self.km_increase)
        km_plus_button.pack(side=tk.LEFT)

    def km_decrease(self):
        current_value = self.km_slider.get()
        new_value = max(current_value - P.km_step, P.km_min)
        self.km_slider.set(new_value)

    def km_increase(self):
        current_value = self.km_slider.get()
        new_value = min(current_value + P.km_step, P.km_max)
        self.km_slider.set(new_value)

    # ----------------Encoder Display Functions--------------------
    def pack_Encoder(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        self.phi_val = tk.StringVar(master, value="0")
        self.theta_val = tk.StringVar(master, value="0")
        self.psi_val = tk.StringVar(master, value="0")

        labels = []
        values = []
        labels.append("Phi: ")
        values.append(self.phi_val)
        labels.append("Theta: ")
        values.append(self.theta_val)
        labels.append("Psi: ")
        values.append(self.psi_val)

        # A frame for the button
        Angle_frame = tk.LabelFrame(master, text="Encoder Values", padx=10, pady=10)
        Angle_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        # Use grid to place labels and their values
        for i, (label_text, var) in enumerate(zip(labels, values)):
            label = tk.Label(Angle_frame, text=label_text, anchor="e", width=8)
            label.grid(row=i, column=0, sticky="E")

            value_entry = tk.Entry(
                Angle_frame,
                textvariable=var,
                state="readonly",
                width=20,
                justify="center",
            )
            value_entry.grid(row=i, column=1, sticky="W")

        # Ensure columns are aligned
        Angle_frame.grid_columnconfigure(
            0, minsize=30
        )  # Adjust minsize for label column width
        Angle_frame.grid_columnconfigure(
            1, minsize=100
        )  # Adjust minsize for value column width

    # ----------------Pitch Slider Functions--------------------
    def pack_Pitch(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        pitch_frame = tk.LabelFrame(master, text="Pitch PID (theta)", padx=10, pady=10)
        pitch_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        self.theta_ref_slider = self.create_slider_control(
            pitch_frame,
            "Theta Ref:",
            self.theta_decrease,
            self.theta_increase,
            P.theta_ref_min,
            P.theta_ref_max,
            P.theta_ref_step,
            P.theta_ref,
        )
        self.theta_kp_slider = self.create_slider_control(
            pitch_frame,
            "Theta Kp:",
            self.theta_kp_decrease,
            self.theta_kp_increase,
            P.theta_kp_min,
            P.theta_kp_max,
            P.theta_kp_step,
            P.theta_kp,
        )
        self.theta_kd_slider = self.create_slider_control(
            pitch_frame,
            "Theta Kd:",
            self.theta_kd_decrease,
            self.theta_kd_increase,
            P.theta_kd_min,
            P.theta_kd_max,
            P.theta_kd_step,
            P.theta_kd,
        )
        self.theta_ki_slider = self.create_slider_control(
            pitch_frame,
            "Theta Ki:",
            self.theta_ki_decrease,
            self.theta_ki_increase,
            P.theta_ki_min,
            P.theta_ki_max,
            P.theta_ki_step,
            P.theta_ki,
        )

    def theta_decrease(self):
        current_value = self.theta_ref_slider.get()
        new_value = max(current_value - P.theta_ref_step, P.theta_ref_min)
        self.theta_ref_slider.set(new_value)

    def theta_increase(self):
        current_value = self.theta_ref_slider.get()
        new_value = min(current_value + P.theta_ref_step, P.theta_ref_max)
        self.theta_ref_slider.set(new_value)

    def theta_kp_decrease(self):
        current_value = self.theta_kp_slider.get()
        new_value = max(current_value - P.theta_kp_step, P.theta_kp_min)
        self.theta_kp_slider.set(new_value)

    def theta_kp_increase(self):
        current_value = self.theta_kp_slider.get()
        new_value = min(current_value + P.theta_kp_step, P.theta_kp_max)
        self.theta_kp_slider.set(new_value)

    def theta_kd_decrease(self):
        current_value = self.theta_kd_slider.get()
        new_value = max(current_value - P.theta_kd_step, P.theta_kd_min)
        self.theta_kd_slider.set(new_value)

    def theta_kd_increase(self):
        current_value = self.theta_kd_slider.get()
        new_value = min(current_value + P.theta_kd_step, P.theta_kd_max)
        self.theta_kd_slider.set(new_value)

    def theta_ki_decrease(self):
        current_value = self.theta_ki_slider.get()
        new_value = max(current_value - P.theta_ki_step, P.theta_ki_min)
        self.theta_ki_slider.set(new_value)

    def theta_ki_increase(self):
        current_value = self.theta_ki_slider.get()
        new_value = min(current_value + P.theta_ki_step, P.theta_ki_max)
        self.theta_ki_slider.set(new_value)

    # ----------------Roll Slider Functions--------------------
    def pack_Roll(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        phiPD_frame = tk.LabelFrame(master, text="Roll PD (phi)", padx=10, pady=10)
        phiPD_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        self.phi_ref_slider = self.create_slider_control(
            phiPD_frame,
            "Phi Ref:",
            self.phi_decrease,
            self.phi_increase,
            P.phi_ref_min,
            P.phi_ref_max,
            P.phi_ref_step,
            P.phi_ref,
        )
        self.phi_kp_slider = self.create_slider_control(
            phiPD_frame,
            "Phi Kp:",
            self.phi_kp_decrease,
            self.phi_kp_increase,
            P.phi_kp_min,
            P.phi_kp_max,
            P.phi_kp_step,
            P.phi_kp,
        )
        self.phi_kd_slider = self.create_slider_control(
            phiPD_frame,
            "Phi Kd:",
            self.phi_kd_decrease,
            self.phi_kd_increase,
            P.phi_kd_min,
            P.phi_kd_max,
            P.phi_kd_step,
            P.phi_kd,
        )

    def phi_decrease(self):
        current_value = self.phi_ref_slider.get()
        new_value = max(current_value - P.phi_ref_step, P.phi_ref_min)
        self.phi_ref_slider.set(new_value)

    def phi_increase(self):
        current_value = self.phi_ref_slider.get()
        new_value = min(current_value + P.phi_ref_step, P.phi_ref_max)
        self.phi_ref_slider.set(new_value)

    def phi_kp_decrease(self):
        current_value = self.phi_kp_slider.get()
        new_value = max(current_value - P.phi_kp_step, P.phi_kp_min)
        self.phi_kp_slider.set(new_value)

    def phi_kp_increase(self):
        current_value = self.phi_kp_slider.get()
        new_value = min(current_value + P.phi_kp_step, P.phi_kp_max)
        self.phi_kp_slider.set(new_value)

    def phi_kd_decrease(self):
        current_value = self.phi_kd_slider.get()
        new_value = max(current_value - P.phi_kd_step, P.phi_kd_min)
        self.phi_kd_slider.set(new_value)

    def phi_kd_increase(self):
        current_value = self.phi_kd_slider.get()
        new_value = min(current_value + P.phi_kd_step, P.phi_kd_max)
        self.phi_kd_slider.set(new_value)

    # ----------------Yaw Slider Functions--------------------
    def pack_Yaw(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        psiPID_frame = tk.LabelFrame(master, text="Yaw PID (psi)", padx=10, pady=10)
        psiPID_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        self.psi_ref_slider = self.create_slider_control(
            psiPID_frame,
            "Psi Ref:",
            self.psi_decrease,
            self.psi_increase,
            P.psi_ref_min,
            P.psi_ref_max,
            P.psi_ref_step,
            P.psi_ref,
        )
        self.psi_kp_slider = self.create_slider_control(
            psiPID_frame,
            "Psi Kp:",
            self.psi_kp_decrease,
            self.psi_kp_increase,
            P.psi_kp_min,
            P.psi_kp_max,
            P.psi_kp_step,
            P.psi_kp,
        )
        self.psi_kd_slider = self.create_slider_control(
            psiPID_frame,
            "Psi Kd:",
            self.psi_kd_decrease,
            self.psi_kd_increase,
            P.psi_kd_min,
            P.psi_kd_max,
            P.psi_kd_step,
            P.psi_kd,
        )
        self.psi_ki_slider = self.create_slider_control(
            psiPID_frame,
            "Psi Ki:",
            self.psi_ki_decrease,
            self.psi_ki_increase,
            P.psi_ki_min,
            P.psi_ki_max,
            P.psi_ki_step,
            P.psi_ki,
        )

    def psi_decrease(self):
        current_value = self.psi_ref_slider.get()
        new_value = max(current_value - P.psi_ref_step, P.psi_ref_min)
        self.psi_ref_slider.set(new_value)

    def psi_increase(self):
        current_value = self.psi_ref_slider.get()
        new_value = min(current_value + P.psi_ref_step, P.psi_ref_max)
        self.psi_ref_slider.set(new_value)

    def psi_kp_decrease(self):
        current_value = self.psi_kp_slider.get()
        new_value = max(current_value - P.psi_kp_step, P.psi_kp_min)
        self.psi_kp_slider.set(new_value)

    def psi_kp_increase(self):
        current_value = self.psi_kp_slider.get()
        new_value = min(current_value + P.psi_kp_step, P.psi_kp_max)
        self.psi_kp_slider.set(new_value)

    def psi_kd_decrease(self):
        current_value = self.psi_kd_slider.get()
        new_value = max(current_value - P.psi_kd_step, P.psi_kd_min)
        self.psi_kd_slider.set(new_value)

    def psi_kd_increase(self):
        current_value = self.psi_kd_slider.get()
        new_value = min(current_value + P.psi_kd_step, P.psi_kd_max)
        self.psi_kd_slider.set(new_value)

    def psi_ki_decrease(self):
        current_value = self.psi_ki_slider.get()
        new_value = max(current_value - P.psi_ki_step, P.psi_ki_min)
        self.psi_ki_slider.set(new_value)

    def psi_ki_increase(self):
        current_value = self.psi_ki_slider.get()
        new_value = min(current_value + P.psi_ki_step, P.psi_ki_max)
        self.psi_ki_slider.set(new_value)

    # ----------------Action Button Functions--------------------
    def pack_Action(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        action_frame = tk.LabelFrame(master, text="Controller", padx=10, pady=10)
        action_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        roll_button = tk.Button(
            action_frame,
            text="Roll",
            command=self.set_roll,
            width=self.w_button_2,
            height=self.h_button,
        )
        roll_button.pack(side=tk.LEFT, padx=5)

        yaw_button = tk.Button(
            action_frame,
            text="Yaw",
            command=self.set_yaw,
            width=self.w_button_2,
            height=self.h_button,
        )
        yaw_button.pack(side=tk.LEFT, padx=5)

    def set_roll(self):
        P.controller = P.controller_roll

        self.reset_references()

        enable_list = [self.phi_ref_slider, self.phi_kp_slider, self.phi_kd_slider]

        disable_list = [
            self.psi_ref_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    def set_yaw(self):
        P.controller = P.controller_yaw

        self.reset_references()

        enable_list = [
            self.psi_ref_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        disable_list = [self.phi_ref_slider, self.phi_kp_slider, self.phi_kd_slider]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    # ----------------Control Button Functions--------------------
    def pack_Control(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        HB_Control_frame = tk.LabelFrame(master, text="HB Control", padx=10, pady=10)
        HB_Control_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        zero_button = tk.Button(
            HB_Control_frame,
            text="Zero",
            command=self.send_zero,
            width=self.w_button_3,
            height=self.h_button,
        )
        zero_button.pack(side=tk.LEFT, padx=5)

        roll_button = tk.Button(
            HB_Control_frame,
            text="Pause",
            command=self.send_pause,
            width=self.w_button_3,
            height=self.h_button,
        )
        roll_button.pack(side=tk.LEFT, padx=5)

        yaw_button = tk.Button(
            HB_Control_frame,
            text="Play",
            command=self.send_play,
            width=self.w_button_3,
            height=self.h_button,
        )
        yaw_button.pack(side=tk.LEFT, padx=5)

    def send_zero(self):
        P.Zero = True

    def send_pause(self):
        P.Pause = True

    def send_play(self):
        P.phi_ref = 0.0
        P.theta_ref = 0.0
        P.psi_ref = 0.0
        P.Pause = False
        P.Play = True

    # ----------------Full Action Button Functions--------------------
    def pack_Action_Full(self, master, Location):
        Row = Location[0]
        Column = Location[1]

        full_action_frame = tk.LabelFrame(master, text="Controller", padx=10, pady=10)
        full_action_frame.grid(row=Row, column=Column, padx=10, pady=10, sticky="nsew")

        roll_button = tk.Button(
            full_action_frame,
            text="Roll",
            command=self.set_roll_full,
            width=self.w_button_4,
            height=self.h_button,
        )
        roll_button.pack(side=tk.LEFT, padx=5)

        pitch_button = tk.Button(
            full_action_frame,
            text="Pitch",
            command=self.set_pitch_full,
            width=self.w_button_4,
            height=self.h_button,
        )
        pitch_button.pack(side=tk.LEFT, padx=5)

        yaw_button = tk.Button(
            full_action_frame,
            text="Yaw",
            command=self.set_yaw_full,
            width=self.w_button_4,
            height=self.h_button,
        )
        yaw_button.pack(side=tk.LEFT, padx=5)

        full_button = tk.Button(
            full_action_frame,
            text="Full",
            command=self.set_full,
            width=self.w_button_4,
            height=self.h_button,
        )
        full_button.pack(side=tk.LEFT, padx=5)

    def set_roll_full(self):
        P.controller = P.controller_roll

        self.reset_references()

        enable_list = [self.phi_ref_slider, self.phi_kp_slider, self.phi_kd_slider]

        disable_list = [
            self.theta_ref_slider,
            self.theta_kp_slider,
            self.theta_ki_slider,
            self.theta_kd_slider,
            self.psi_ref_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    def set_pitch_full(self):
        P.controller = P.controller_pitch

        self.reset_references()

        enable_list = [
            self.theta_ref_slider,
            self.theta_kp_slider,
            self.theta_ki_slider,
            self.theta_kd_slider,
        ]

        disable_list = [
            self.phi_ref_slider,
            self.phi_kp_slider,
            self.phi_kd_slider,
            self.psi_ref_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    def set_yaw_full(self):
        P.controller = P.controller_yaw

        self.reset_references()

        enable_list = [
            self.psi_ref_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        disable_list = [
            self.phi_ref_slider,
            self.phi_kp_slider,
            self.phi_kd_slider,
            self.theta_ref_slider,
            self.theta_kp_slider,
            self.theta_ki_slider,
            self.theta_kd_slider,
        ]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    def set_full(self):
        P.controller = P.controller_full

        self.reset_references()

        enable_list = [
            self.theta_ref_slider,
            self.psi_ref_slider,
        ]

        disable_list = [
            self.phi_ref_slider,
            self.phi_kp_slider,
            self.phi_kd_slider,
            self.theta_kp_slider,
            self.theta_ki_slider,
            self.theta_kd_slider,
            self.psi_kp_slider,
            self.psi_ki_slider,
            self.psi_kd_slider,
        ]

        for slider in enable_list:
            slider["state"] = "normal"
            slider["troughcolor"] = self.enabled_color

        for slider in disable_list:
            slider["state"] = "disabled"
            slider["troughcolor"] = self.disabled_color

    ##########################################################
    ##################--Helper Function--#####################
    ##########################################################

    def create_slider_control(
        self,
        master,
        label_text,
        command_decrease,
        command_increase,
        from_,
        to_,
        resolution,
        start,
    ):
        frame = tk.Frame(master)
        frame.pack(side=tk.TOP, fill=tk.X, expand=True)

        label = tk.Label(frame, text=label_text)
        label.pack(side=tk.TOP)

        minus_button = tk.Button(frame, text="-", command=command_decrease)
        minus_button.pack(side=tk.LEFT)

        slider = tk.Scale(
            frame,
            length=200,
            from_=from_,
            to=to_,
            resolution=resolution,
            orient="horizontal",
        )
        slider.set(start)
        slider.pack(side=tk.LEFT)
        slider["troughcolor"] = self.enabled_color

        plus_button = tk.Button(frame, text="+", command=command_increase)
        plus_button.pack(side=tk.LEFT)

        return slider

    # ----------------Send/receive data--------------------
    def GUI_Values(self, phi, theta, psi):
        if self.Enabled_Km:
            P.km = float(self.km_slider.get())
        if self.Enabled_Pitch:
            P.theta_ref = float(self.theta_ref_slider.get())
            P.theta_kp = float(self.theta_kp_slider.get())
            P.theta_ki = float(self.theta_ki_slider.get())
            P.theta_kd = float(self.theta_kd_slider.get())
        if self.Enabled_Roll:
            P.phi_ref = float(self.phi_ref_slider.get())
            P.phi_kp = float(self.phi_kp_slider.get())
            P.phi_kd = float(self.phi_kd_slider.get())
        if self.Enabled_Yaw:
            P.psi_ref = float(self.psi_ref_slider.get())
            P.psi_kp = float(self.psi_kp_slider.get())
            P.psi_ki = float(self.psi_ki_slider.get())
            P.psi_kd = float(self.psi_kd_slider.get())
        if self.Enabled_Encoder:
            self.phi_val.set(str(round(phi * 180 / np.pi, 1)))
            self.theta_val.set(str(round(theta * 180 / np.pi, 1)))
            self.psi_val.set(str(round(psi * 180 / np.pi, 1)))

    def Flip_Theta(self):
        P.theta_ref = -P.theta_ref
        self.theta_ref_slider.set(P.theta_ref)

    def Flip_Psi(self):
        P.psi_ref = -P.psi_ref
        self.psi_ref_slider.set(P.psi_ref)

    def Flip_Phi(self):
        P.phi_ref = -P.phi_ref
        self.phi_ref_slider.set(P.phi_ref)

    def reset_references(self):
        try:
            P.phi_ref = 0
            self.phi_ref_slider.set(0)
        except:
            pass
        try:
            P.theta_ref = 0
            self.theta_ref_slider.set(0)
        except:
            pass
        try:
            P.psi_ref = 0
            self.psi_ref_slider.set(0)
        except:
            pass
