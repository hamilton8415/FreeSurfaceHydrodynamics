#!/usr/bin/pyton3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fshd import FS_HydroDynamics


def main():
    import argparse
    import importlib.resources as ilr
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=float, default=1.0,
                        help='(meters) sets the incident wave amplitude')
    args = parser.parse_args()

    A = args.a
    Tp = 10.
    omega = 2. * np.pi / Tp
    phase = 0. * np.pi / 180.

    modes = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg

    tf = 3. * Tp
    dt = 0.015

    BuoyA5 = FS_HydroDynamics()

    # Set area and 2nd moments of area for waterplane
    BuoyA5.SetWaterplane(5.47, 1.37, 1.37)
    # Set COB relative to waterplane coordinate system.
    BuoyA5.SetCOB(0., 0., -0.22)
    BuoyA5.SetVolume(buoy_mass / rho)

    pts_t = dt * np.arange(0., tf / dt)
    pts_pos = A * np.cos(omega * pts_t + phase)
    pts_vel = -A * omega * np.sin(omega * pts_t + phase)
    pts_accel = -A * omega**2. * np.cos(omega * pts_t + phase)

    fig, ax = plt.subplots(1)

    ax.plot(pts_t, pts_vel, label='Vel (m/s)')
    ax.plot(pts_t, pts_accel, label='Acc (m/s**2)')
    ax.set_xlabel('time (s)')
    ax.grid()
    ax.legend()

    fig.suptitle(f'{A = :.1f}m {Tp = :.1f}s')
    fig.tight_layout()

    print('Computing Buoyancy Forces')
    # j denotes direction of resulting force
    for j in range(6):
        pts_x = []
        pts_F_B = []
        for k in range(pts_pos.shape[0]):
            x = np.zeros(6, dtype=np.float64)
            x[j] = pts_pos[k]
            pts_x.append(x[j])
            BuoyancyForce = BuoyA5.BuoyancyForce(x)
            pts_F_B.append(BuoyancyForce[j])

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
        ax2.plot(pts_t, pts_F_B, color=color)
        ax2.set_ylabel(ylabel, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.suptitle(f'{modes[j]} Buoyancy Forces')
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
