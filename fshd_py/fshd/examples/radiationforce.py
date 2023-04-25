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
    parser.add_argument('-t', type=float, default=5.0,
                        help='(seconds) sets the incident wave period')
    parser.add_argument('-p', type=float, default=40.0,
                        help='(degrees) sets the incident wave phase angle')
    args = parser.parse_args()

    A = args.a
    Tp = args.t
    phase = args.p * np.pi / 180.0

    modes = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg
    BuoyA5 = FS_HydroDynamics(1.0, g, rho)

    omega = 2. * np.pi / Tp
    tf = 3. * Tp
    dt = 0.01

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
    
    fig, ax = plt.subplots(1)

    ax.plot(pts_t, pts_vel, label='Vel (m/s)')
    ax.plot(pts_t, pts_accel, label='Acc (m/s**2)')
    ax.set_xlabel('time (s)')
    ax.grid()
    ax.legend()

    fig.suptitle(f'{A = :.1f}m {Tp = :.1f}s')
    fig.tight_layout()


    # Test Radiation Forces
    # Note:  This computes the radiation forces for each mode of motion
    # indivdually, so MemRadiation is called 6 times per timestep, once with each
    # acceleration set non-zero In use it will be called only once per time-step,
    # with all acclerations set.
    for idx in range(6):  # determines mode of motion.
        for jdx in range(6):  # denotes direction of resulting force
            am = BuoyA5.AddedMass(omega, idx, jdx)
            dmp = BuoyA5.RadiationDamping(omega, idx, jdx)
            am_inf = BuoyA5.fd_A_inf_freq[idx, jdx]

            pts_F_TD = []
            pts_F_FD = []
            last_accel = 0.
            F_max = -np.inf
            F_min = np.inf
            BuoyA5.SetTimestepSize(dt)  # Reset timestep to re-initialize storage in BuoyA5 Class
            for kdx in range(pts_accel.shape[0]):
                accel = pts_accel[kdx]
                xddot = np.zeros(6)
                xddot[idx] = last_accel

                MemForce = BuoyA5.RadiationForce(xddot)
                pts_F_TD.append(am_inf * accel + MemForce[jdx])
                last_accel = accel
                FD_Force = am * pts_accel[kdx] + dmp * pts_vel[kdx]
                pts_F_FD.append(FD_Force)
                if FD_Force > F_max:
                    F_max = FD_Force
                if FD_Force < F_min:
                    F_min = FD_Force

            if F_min < -1. and F_max > 1.:  # Don't plot near-zero forces
                fig, ax = plt.subplots(1)

                suptitle = f'\n{modes[idx]}(t) = {A = :.1f}cos(2 pi t/{Tp = :.1f})'
                axtitle = f'Radiation Forces: {modes[idx]} Motions' \
                          + f', {modes[jdx]} Forces:  {A = :.1f}m' \
                          + f' {Tp = :.1f}s'

                if jdx < 3:
                    ylabel = 'F (N)'
                else:
                    ylabel = 'M (N-m)'

                ax.plot(pts_t, pts_F_TD, 'r', label='Time-Domain')
                ax.plot(pts_t, pts_F_FD, '--g', label='Freq-Domain')
                ax.set_xlabel('time (s)')
                ax.set_ylabel(ylabel)
                ax.set_title(axtitle)
                ax.grid()
                ax.legend()

                fig.suptitle(suptitle)
                fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
