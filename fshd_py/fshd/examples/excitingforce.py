#!/usr/bin/pyton3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fshd import FS_HydroDynamics, LinearIncidentWave


def main():
    import argparse
    import importlib.resources as ilr
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=float, default=1.0,
                        help='(meters) sets the incident wave amplitude')
    parser.add_argument('-t', type=float, default=5.0,
                        help='(seconds) sets the incident wave period')
    parser.add_argument('-p', type=float, default=40.0,
                        help='(degrees) sets the incident wave phase angle')
    parser.add_argument('-b', type=float, default=180.0,
                        help='(degrees) sets the incident wave direction')
    args = parser.parse_args()

    A = args.a
    Tp = args.t
    phase = args.p * np.pi / 180.0
    beta = args.b * np.pi / 180.0

    modes = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    Inc = LinearIncidentWave()

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg
    BuoyA5 = FS_HydroDynamics()

    omega = 2. * np.pi / Tp
    tf = 3. * Tp
    dt = 0.015

    Inc.SetToMonoChromatic(A, Tp, phase, beta)

    BuoyA5.AssignIncidentWave(Inc)

    # Set area and 2nd moments of area for waterplane
    BuoyA5.SetWaterplane(5.47, 1.37, 1.37)
    # Set COB relative to waterplane coordinate system.
    BuoyA5.SetCOB(0., 0., -0.22)
    #Set COG relative to waterplane coordinate system.
    BuoyA5.SetCOG(0., 0., -0.24)
    BuoyA5.SetVolume(buoy_mass / rho)
    BuoyA5.SetMass(buoy_mass)

    HydrodynamicsBaseFilename = \
        ilr.files('fshd.examples').joinpath('example_hydrodynamic_coeffs',
                                            'BuoyA5')
    HydrodynamicsBaseFilename = str(HydrodynamicsBaseFilename)
    BuoyA5.ReadWAMITData_FD(HydrodynamicsBaseFilename)
    BuoyA5.ReadWAMITData_TD(HydrodynamicsBaseFilename)
    BuoyA5.SetTimestepSize(dt)

    pts_t = dt * np.arange(0., tf / dt)
    pts_pos = A * np.cos(omega * pts_t + phase)
    pts_vel = -A * omega * np.sin(omega * pts_t + phase)
    pts_accel = -A * omega**2. * np.cos(omega * pts_t + phase)

    # j denotes direction of resulting force
    for j in range(6):
        Chi = BuoyA5.WaveExcitingForceComponents(Inc.m_omega[0], j)
        pts_F_TD = []
        pts_F_FD = []
        pts_eta = []
        BuoyA5.SetTimestepSize(dt)
        for k in range(pts_t.shape[0]):
            pts_F_FD.append(
                Inc.m_A[0] * Chi.real *
                    np.cos(omega * pts_t[k] + phase * np.cos(180. * np.pi / 180.)) -
                Inc.m_A[0] * Chi.imag *
                    np.sin(omega * pts_t[k] + phase * np.cos(180. * np.pi / 180.)))
            pts_eta.append(Inc.eta(0., 0., pts_t[k]))
            ExtForce = BuoyA5.ExcitingForce()
            pts_F_TD.append(ExtForce[j])

        # plot
        fig, ax = plt.subplots(3)

        if j < 3:
            ylabel = 'F (N)'
        else:
            ylabel = 'M (N-m)'

        ax[0].plot(pts_t, pts_F_TD)
        ax[0].set_ylabel(ylabel)
        ax[0].set_xlabel('time (s)')
        ax[0].set_title('Time-Domain')

        ax[1].plot(pts_t, pts_F_FD)
        ax[1].set_ylabel(ylabel)
        ax[1].set_xlabel('time (s)')
        ax[1].set_title('Freq-Domain')

        ax[2].plot(pts_t, pts_eta)
        ax[2].set_ylabel(ylabel)
        ax[2].set_xlabel('time (s)')
        ax[2].set_title('eta(t)')

        fig.suptitle(f'{modes[j]} Exciting Forces')
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
