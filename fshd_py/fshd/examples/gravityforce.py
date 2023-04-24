#!/usr/bin/pyton3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fshd import FS_HydroDynamics


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=float, default=1.0,
                        help='(meters) sets the incident wave amplitude')
    args = parser.parse_args()

    A = args.a
    Tp = 5.
    omega = 2. * np.pi / Tp
    phase = 0. * np.pi / 180.

    modes = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg

    tf = 2. * Tp
    dt = 0.01

    BuoyA5 = FS_HydroDynamics()

    # Set area and 2nd moments of area for waterplane
    BuoyA5.SetWaterplane(5.47, 1.37, 1.37)
    BuoyA5.SetCOG(0., 0.,
                  -0.24)  # Set COG relative to waterplane coordinate system.
    BuoyA5.SetMass(buoy_mass)

    pts_t = dt * np.arange(0., tf / dt)
    pts_pos = A * np.cos(omega * pts_t + phase)
    pts_vel = -A * omega * np.sin(omega * pts_t + phase)
    pts_accel = -A * omega**2. * np.cos(omega * pts_t + phase)

    print('Computing Gravity Forces')
    # j denotes direction of resulting force
    for j in range(6):
        pts_x = []
        pts_F_G = []
        for k in range(pts_pos.shape[0]):
            x = np.zeros(6, dtype=np.float64)
            x[j] = pts_pos[k]
            pts_x.append(x[j])
            GravityForce = BuoyA5.GravityForce(x)
            pts_F_G.append(GravityForce[j])

        # plot
        fig, ax = plt.subplots(1)

        if j < 3:
            ylabel = 'F (N)'
        else:
            ylabel = 'M (N-m)'

        color = 'tab:blue'
        ax.plot(pts_t, pts_x, color=color)
        ax.set_ylabel('x(t) meters', color=color)
        ax.set_xlabel('time (s)')
        ax.tick_params(axis='y', labelcolor=color)

        color = 'tab:red'
        ax2 = ax.twinx()
        ax2.plot(pts_t, pts_F_G, color=color)
        ax2.set_ylabel(ylabel, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.suptitle(f'{modes[j]} Gravity Forces')
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
