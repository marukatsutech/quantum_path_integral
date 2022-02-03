# Quantum Mechanics, Path integral
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d


def update_curve_phi1():
    global curve_phi1
    y = amp1 * np.sin(phase1)
    z = amp1 * np.cos(phase1)
    curve_phi1.set_xdata(x)
    curve_phi1.set_ydata(y)
    curve_phi1.set_3d_properties(z)


def update_curve_phi0():
    global curve_phi0
    y = amp0 * np.sin(phase0)
    z = amp0 * np.cos(phase0)
    curve_phi0.set_xdata(x)
    curve_phi0.set_ydata(y)
    curve_phi0.set_3d_properties(z)


def update_curve_prob1():
    global curve_prob1
    curve_prob1.set_xdata(x)
    curve_prob1.set_ydata(x * 0.)
    curve_prob1.set_3d_properties(prob1 * prob_magnify)


def update_curve_prob0():
    global curve_prob0
    curve_prob0.set_xdata(x)
    curve_prob0.set_ydata(x * 0.)
    curve_prob0.set_3d_properties(prob0 * prob_magnify)


def set_phi0():
    global amp0, phase0, amp1
    amp0 = np.sqrt(prob0)
    if is_phase:
        if delta_x != 0.:
            phase0 = (x / (x_max - x_min) * 2 * np.pi) * ((x_max - x_min) / delta_x) * wave_number
        else:
            phase0 = x * 0.
    else:
        phase0 = x * 0.


def set_prob0():
    global curve_prob0, prob0
    if num_of_points_in_delta_x <= 1:
        for i in range(num_of_points):
            prob0[i] = 0.
        prob0[int((position - x_min) / (x_max - x_min) * num_of_points)] = 1.
    else:
        for i in range(num_of_points):
            if (position - delta_x / 2) <= x[i] <= (position + delta_x / 2):
                prob0[i] = 1. / num_of_points_in_delta_x
            else:
                prob0[i] = 0.


def update_additional_lines():
    global line_position, line_additional1, line_additional2, line_additional3
    line_position.set_data_3d([position, position], [0., 0.], [z_min * yz_scale, z_max * yz_scale])
    line_additional1.set_data_3d([position - delta_x / 2., position - delta_x / 2.], [0., 0.],
                                 [0., z_max * yz_scale * 1.2])
    line_additional2.set_data_3d([position + delta_x / 2., position + delta_x / 2.], [0., 0.],
                                 [0., z_max * yz_scale * 1.2])
    line_additional3.set_data_3d([position - delta_x / 2., position + delta_x / 2.], [0., 0.],
                                 [z_max * yz_scale * 1.15, z_max * yz_scale * 1.15])


def reset_phi1_prob1():
    global is_play, tx_step, cnt, prob1
    is_play = False
    cnt = 0
    tx_step.set_text("Step=" + str(cnt))
    curve_phi1.set_xdata(x * 0.)
    curve_phi1.set_ydata(x * 0.)
    curve_phi1.set_3d_properties(x * 0.)
    update_curve_phi0()
    prob1 = x * 0.
    update_curve_prob1()


def change_mass(m):
    global mass
    reset_phi1_prob1()
    mass = m


def change_wave_number(wn):
    global is_play, wave_number
    reset_phi1_prob1()
    wave_number = wn
    set_phi0()
    update_curve_phi0()


def switch_phase():
    global is_play, is_phase
    reset_phi1_prob1()
    if var_phase.get():
        is_phase = True
    else:
        is_phase = False
    set_phi0()
    update_curve_phi0()


def change_dt(dt):
    global is_play, delta_t
    reset_phi1_prob1()
    delta_t = dt


def change_dx(dx):
    global is_play, delta_x, num_of_points_in_delta_x, tx_num_pnt
    reset_phi1_prob1()
    delta_x = dx
    num_of_points_in_delta_x = int(num_of_points * delta_x / (x_max - x_min))
    if num_of_points_in_delta_x == 0:
        num_of_points_in_delta_x = 1
    tx_num_pnt.set_text("Number of points in dx=" + str(num_of_points_in_delta_x))
    tx_delta_x.set_text("dx=" + str(dx))
    update_additional_lines()
    set_prob0()
    update_curve_prob0()
    set_phi0()
    update_curve_phi0()


def change_position(pst):
    global is_play, position
    reset_phi1_prob1()
    position = pst
    update_additional_lines()
    xxz, yyz, _ = proj3d.proj_transform(position, 0., z_max * yz_scale * 1.2, ax1.get_proj())
    tx_delta_x.set_position((xxz, yyz))
    set_prob0()
    update_curve_prob0()
    set_phi0()
    update_curve_phi0()


def change_magnify(mg):
    global prob_magnify
    prob_magnify = mg
    update_curve_prob0()
    update_curve_prob1()


def change_yz_scale(sc):
    global y_min, y_max, z_min, z_max, yz_scale
    yz_scale = sc
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min * yz_scale, y_max * yz_scale)
    ax1.set_zlim(z_min * yz_scale, z_max * yz_scale)
    update_additional_lines()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt, tx_step, phase1, amp1, prob1, curve_phi1, curve_prob1
    tx_step.set_text("Step=" + str(cnt))
    if is_play:
        y1 = x * 0.
        z1 = x * 0.
        for i in range(num_of_points):
            if cnt != 0:
                # w = mx**2/2ht
                # Reference
                # THE QUANTUM UNIVERSE
                # (and why anything that can happen, does)
                # by Brian Cox and Jeff Forshaw
                m_per_2h = mass / (2 * h) * 10 ** (power_mass - power_h)
                w = m_per_2h * np.abs(x - x[i]) ** 2 / (cnt * delta_t)
                ph = (phase0 + 2. * np.pi * w) % (2. * np.pi)
                # Superpose the vectors of each points
                yy = amp0 * np.sin(ph) / num_of_points
                zz = amp0 * np.cos(ph) / num_of_points
                y1[i] = np.sum(yy)
                z1[i] = np.sum(zz)
        amp1 = np.sqrt(x * 0. + y1 ** 2 + z1 ** 2)
        phase1 = np.arctan2(y1, z1)
        update_curve_phi1()
        prob1 = amp1 ** 2.
        # Calibrate probability
        prob_cal = np.sum(prob1)
        if prob_cal != 0:
            prob1 = amp1 ** 2 / prob_cal
        update_curve_prob1()
        cnt += 1


# Global variables
num_of_points = 2000

x_min = 0.
x_max = 20.
y_min = -1.
y_max = 1.
z_min = -1.
z_max = 1.

is_play = False
cnt = 0
yz_scale_init = 0.1
yz_scale = yz_scale_init
prob_magnify_init = 1.
prob_magnify = prob_magnify_init

position_init = (x_max - x_min) / 2.
position = position_init
delta_x_init = 1.0
delta_x = delta_x_init
num_of_points_in_delta_x = int(num_of_points * delta_x / (x_max - x_min))
wave_number_init = 4.
wave_number = wave_number_init
is_phase = False

mass_init = 100.  # electron 9.10938 10**-31 kg
power_mass = -31
mass = mass_init
h = 6.626   # Planck constant 6.626 * 10**-34 kg*m**2 / s
power_h = -34
delta_t_init = 1000.
delta_t = delta_t_init

# Generate figure and axes
fig = Figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.set_box_aspect((2, 1, 1))
ax1.grid()
ax1.set_title('Quantum Mechanics, Path integral w = mx**2/2ht')
ax1.set_xlabel('x')
ax1.set_ylabel('phi (imaginary)')
ax1.set_zlabel('phi (real) or Probability')
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min * yz_scale, y_max * yz_scale)
ax1.set_zlim(z_min * yz_scale, z_max * yz_scale)

# Generate items
tx_step = ax1.text2D(x_min, y_max, "Step=" + str(0))
xz, yz, _ = proj3d.proj_transform(x_min, y_max * yz_scale, z_max * yz_scale * 1.2, ax1.get_proj())
tx_step.set_position((xz, yz))
tx_num_pnt = ax1.text2D(x_min, y_max, "Number of points in delta x=" + str(num_of_points_in_delta_x))
xz, yz, _ = proj3d.proj_transform(x_min, y_max * yz_scale, z_max * yz_scale, ax1.get_proj())
tx_num_pnt.set_position((xz, yz))
tx_delta_x = ax1.text2D(position, 0., "dx=" + str(delta_x_init))
xz, yz, _ = proj3d.proj_transform(position, 0., z_max * yz_scale * 1.2, ax1.get_proj())
tx_delta_x.set_position((xz, yz))

x = np.linspace(x_min, x_max, num_of_points)

prob0 = x * 0.
curve_prob0, = ax1.plot(x, x * 0., prob0, linewidth=2, label='Probability at t=0')
prob1 = x * 0.
curve_prob1, = ax1.plot(x, x * 0., prob1, linewidth=2, label='Probability at t')
amp0 = x * 0.
phase0 = x * 0.
y_phi0 = x * 0.
z_phi0 = x * 0.
curve_phi0, = ax1.plot(x, y_phi0, z_phi0, linewidth=1, label='phi at t=0')
amp1 = x * 0.
phase1 = x * 0.
y_phi1 = x * 0.
z_phi1 = x * 0.
curve_phi1, = ax1.plot(x, y_phi1, z_phi1, linewidth=1, label='phi at t')

ax1.legend()

set_prob0()
update_curve_prob0()
set_phi0()
update_curve_phi0()

line_position = art3d.Line3D([position, position], [0., 0.], [z_min * yz_scale, z_max * yz_scale],
                             color='red', ls="--", linewidth=1)
ax1.add_line(line_position)
line_additional1 = art3d.Line3D([position - delta_x / 2., position - delta_x / 2.], [0., 0.],
                                [0., z_max * yz_scale * 1.2], color='black', ls="--", linewidth=0.5)
ax1.add_line(line_additional1)
line_additional2 = art3d.Line3D([position + delta_x / 2., position + delta_x / 2.], [0., 0.],
                                [0., z_max * yz_scale * 1.2], color='black', ls="--", linewidth=0.5)
ax1.add_line(line_additional2)
line_additional3 = art3d.Line3D([position - delta_x / 2., position + delta_x / 2.], [0., 0.],
                                [z_max * yz_scale * 1.15, z_max * yz_scale * 1.15], color='black', ls="--", linewidth=1)
ax1.add_line(line_additional3)


# Embed in Tkinter
root = tk.Tk()
root.title("Quantum Mechanics, Path integral")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

frm1 = ttk.Labelframe(root, relief="ridge", text="Properties of a particle", labelanchor="n", width=100)
frm1.pack(side='left')
lbl_mass = tk.Label(frm1, text="mass(*10**"+str(power_mass)+")")
lbl_mass.pack(side='left')
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_position = tk.Spinbox(
    frm1, textvariable=var_mass, format="%.1f", from_=1, to=1000, increment=10.0,
    command=lambda: change_mass(float(var_position.get())), width=6
    )
spn_position.pack(side='left')
lbl_position = tk.Label(frm1, text=", x")
lbl_position.pack(side='left')
var_position = tk.StringVar(root)  # variable for spinbox-value
var_position.set(position_init)  # Initial value
spn_position = tk.Spinbox(
    frm1, textvariable=var_position, format="%.1f", from_=x_min, to=x_max, increment=0.1,
    command=lambda: change_position(float(var_position.get())), width=4
    )
spn_position.pack(side='left')
lbl_dx = tk.Label(frm1, text=", dx")
lbl_dx.pack(side='left')
var_dx = tk.StringVar(root)  # variable for spinbox-value
var_dx.set(delta_x_init)  # Initial value
spn_dx = tk.Spinbox(
    frm1, textvariable=var_dx, format="%.2f", from_=0.0, to=10.0, increment=0.01,
    command=lambda: change_dx(float(var_dx.get())), width=4
    )
spn_dx.pack(side='left')

frm2 = ttk.Labelframe(root, relief="ridge", text="Apply phase", labelanchor="n", width=100)
frm2.pack(side='left')
lbl_wn = tk.Label(frm2, text="wave number")
lbl_wn.pack(side='left')
var_wn = tk.StringVar(root)  # variable for spinbox-value
var_wn.set(wave_number_init)  # Initial value
spn_wn = tk.Spinbox(
    frm2, textvariable=var_wn, format="%.1f", from_=-10.0, to=10.0, increment=1.0,
    command=lambda: change_wave_number(float(var_wn.get())), width=4
    )
spn_wn.pack(side='left')
var_phase = tk.BooleanVar(root, value=False)
chk_phase = tk.Checkbutton(frm2, text="Phase", variable=var_phase, command=switch_phase)
chk_phase.pack(side='left')

lbl_dt = tk.Label(root, text="dt per step")
lbl_dt.pack(side='left')
var_dt = tk.StringVar(root)  # variable for spinbox-value
var_dt.set(delta_t_init)  # Initial value
spn_dt = tk.Spinbox(
    root, textvariable=var_dt, format="%.1f", from_=1, to=10000., increment=10.,
    command=lambda: change_dt(float(var_scale.get())), width=8
    )
spn_dt.pack(side='left')

lbl_scale = tk.Label(root, text=", Scale of y,z")
lbl_scale.pack(side='left')
var_scale = tk.StringVar(root)  # variable for spinbox-value
var_scale.set(yz_scale_init)  # Initial value
spn_scale = tk.Spinbox(
    root, textvariable=var_scale, format="%.2f", from_=0.01, to=1., increment=0.01,
    command=lambda: change_yz_scale(float(var_scale.get())), width=4
    )
spn_scale.pack(side='left')

lbl_magnify = tk.Label(root, text=", Magnify of probability")
lbl_magnify.pack(side='left')
var_magnify = tk.StringVar(root)  # variable for spinbox-value
var_magnify.set(prob_magnify_init)  # Initial value
spn_magnify = tk.Spinbox(
    root, textvariable=var_magnify, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_magnify(float(var_magnify.get())), width=4
    )
spn_magnify.pack(side='left')

btn_play = tk.Button(root, text="Play/Pause", command=switch)
btn_play.pack(side='left')
btn_reset = tk.Button(root, text="Reset", command=reset_phi1_prob1)
btn_reset.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=100)
root.mainloop()
