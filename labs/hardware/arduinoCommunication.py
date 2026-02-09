import hummingbirdParam as P
import numpy as np
import time


def read_from_arduino(ser):
    while True:
        if ser.inWaiting() > 0:
            try:
                x = ser.read().decode()
            except:
                print("Bad Input From Arduino")
                x = ""
            if x == chr(P.startMarker):
                P.reading = True
                P.bytesRec = 0
                P.tempBuffer = ""

            if P.reading:
                P.tempBuffer = P.tempBuffer + x
                P.bytesRec += 1

            if x == chr(P.endMarker):
                P.reading = False
                stringRec = ""
                for i in range(1, P.bytesRec - 1):
                    stringRec += P.tempBuffer[i]

                f_pos = stringRec.find("F")
                t_pos = stringRec.find("T")
                s_pos = stringRec.find("S")

                if f_pos != -1 and t_pos != -1 and s_pos != -1:
                    P.phi = float(stringRec[f_pos + 1 : t_pos])
                    P.theta = float(stringRec[t_pos + 1 : s_pos])
                    P.psi = float(stringRec[s_pos + 1 :])


def send_gains(ser):
    if P.Pause and not P.Play:
        GUI_send = (
            chr(P.startMarker)
            + P.char_pause
            + P.char_pause
            + P.char_pause
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)

    if P.Play and not P.Pause:
        GUI_send = (
            chr(P.startMarker)
            + P.char_play
            + P.char_play
            + P.char_play
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)
        P.Play = False

    if P.Zero:
        GUI_send = (
            chr(P.startMarker)
            + P.char_zero
            + P.char_zero
            + P.char_zero
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)
        P.Zero = False

    if P.controller == P.controller_km and not P.Pause:
        GUI_send = (
            chr(P.startMarker) + P.char_km + P.char_val1 + str(P.km) + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)

    if P.controller == P.controller_roll and not P.Pause:
        GUI_send = (
            chr(P.startMarker)
            + P.char_roll
            + P.char_val1
            + str(P.km)
            + P.char_val2
            + str(round(P.phi_ref * np.pi / 180, 4))
            + P.char_val3
            + str(P.phi_kp)
            + P.char_val4
            + str(P.phi_kd)
            + P.char_val5
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)

    if P.controller == P.controller_pitch and not P.Pause:
        GUI_send = (
            chr(P.startMarker)
            + P.char_pitch
            + P.char_val1
            + str(P.km)
            + P.char_val2
            + str(round(P.theta_ref * np.pi / 180, 4))
            + P.char_val3
            + str(P.theta_kp)
            + P.char_val4
            + str(P.theta_ki)
            + P.char_val5
            + str(P.theta_kd)
            + P.char_val6
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)

    if P.controller == P.controller_yaw and not P.Pause:
        GUI_send = (
            chr(P.startMarker)
            + P.char_yaw
            + P.char_val1
            + str(P.km)
            + P.char_val2
            + str(round(P.psi_ref * np.pi / 180, 4))
            + P.char_val3
            + str(P.psi_kp)
            + P.char_val4
            + str(P.psi_ki)
            + P.char_val5
            + str(P.psi_kd)
            + P.char_val6
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)

    if P.controller == P.controller_full and not P.Pause:
        GUI_send = (
            chr(P.startMarker)
            + P.char_full
            + P.char_val1
            + str(round(P.theta_ref * np.pi / 180, 4))
            + P.char_val2
            + str(round(P.psi_ref * np.pi / 180, 4))
            + P.char_val3
            + chr(P.endMarker)
        )
        send_msg(ser, GUI_send)


def send_msg(ser, GUI_msg):
    if GUI_msg != P.GUI_msg_d2:
        byte_msg = GUI_msg.encode()
        ser.write(byte_msg)
        P.GUI_msg_d2 = P.GUI_msg_d1
        P.GUI_msg_d1 = GUI_msg
        print(GUI_msg)
