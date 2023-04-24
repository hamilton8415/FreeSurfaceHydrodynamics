#!/usr/bin/python3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import scipy.integrate as spinteg
import ode

from fshd import LinearIncidentWave, FS_HydroDynamics


# The rhs of x' = f(x) defined as a class
class SingleModeMotionRHS(object):
    def __init__(self, body):
        self.FloatingBody = body
        self.last_accel = 0
        self.mode = 0

    # x[0] = position
    # x[1] = velocity
    def __call__(self, t, x):
        dxdt = np.zeros_like(x, dtype=np.float64)
        pos = np.zeros((6, 1), dtype=np.float64)
        pos[self.mode] = x[0]
        vel = np.zeros((6, 1), dtype=np.float64)
        vel[self.mode] = x[1]
        F_LinDamping = self.FloatingBody.LinearDampingForce(vel)
        F_B = self.FloatingBody.BuoyancyForce(pos)
        F_G = self.FloatingBody.GravityForce(pos)
        accel = np.zeros((6, 1), dtype=np.float64)
        accel[self.mode] = self.last_accel
        F_R = -self.FloatingBody.RadiationForce(accel)
        F_E = self.FloatingBody.ExcitingForce()
        dxdt[0] = x[1]
        b = 0.1
        dxdt[1] = \
            (F_LinDamping[self.mode] + F_B[self.mode] + F_G[self.mode] \
             + F_R[self.mode] + F_E[self.mode]) \
            / (self.FloatingBody.M[self.mode, self.mode] \
               + self.FloatingBody.AddedMass(10000.0, self.mode, self.mode))
        self.last_accel = dxdt[1]

        return dxdt


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
    parser.add_argument('-b', type=float, default=20.0,
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
    tf = 2.0 * Tp
    dt = 0.005

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

    I = np.diag([1500., 1500., 650.]).astype(np.float64)
    BuoyA5.SetI(I)

    b = np.array([300., 300., 900., 400., 400., 100.], dtype=np.float64)
    BuoyA5.SetDampingCoeffs(b)

    pts_T = []
    pts_XiMod = [[] for _ in range(6)]
    pts_XiPh = [[] for _ in range(6)]

    dT = 0.01
    for T in np.arange(dT, 24, dT):
        pts_T.append(T)
        w = 2. * np.pi / T
        Xi = BuoyA5.ComplexAmplitude(w)
        for jdx in range(6):
            pts_XiMod[jdx].append(np.abs(Xi[jdx]))
            pts_XiPh[jdx].append(np.angle(Xi[jdx]) * 180. / np.pi)

    fig_amp, ax_amp = plt.subplots(2, 3)
    fig_amp.suptitle(f'RAO Amplitude')
    fig_ang, ax_ang = plt.subplots(2, 3)
    fig_ang.suptitle(f'RAO Phase Angle (deg)')
    for jdx in range(6):
        ax_amp_jdx = ax_amp[jdx // 3][jdx % 3]
        ax_amp_jdx.plot(pts_T, pts_XiMod[jdx])
        ax_amp_jdx.set_xlabel('Wave Period (s)')
        ax_amp_jdx.set_title(f'{modes[jdx]}')

        ax_ang_jdx = ax_ang[jdx // 3][jdx % 3]
        ax_ang_jdx.plot(pts_T, pts_XiPh[jdx])
        ax_ang_jdx.set_xlabel('Wave Period (s)')
        ax_ang_jdx.set_title(f'{modes[jdx]}')

        if jdx < 6:
            ax_amp_jdx.set_ylabel('meters')
            ax_ang_jdx.set_ylabel('meters')
        else:
            ax_amp_jdx.set_ylabel('degrees')
            ax_ang_jdx.set_ylabel('degrees')

    fig_amp.tight_layout()
    # fig_amp.savefig('motion_rao_amp.png')  # , bbox_inches='tight')
    fig_ang.tight_layout()
    # fig_ang.savefig('motion_rao_ang.png')  # , bbox_inches='tight')

    BuoyA5.SetTimestepSize(dt)

    for mode in range(len(modes)):
        Xi = BuoyA5.ComplexAmplitude(2. * np.pi / Tp, mode)

        x0 = np.array([0.,  # initial position
                       0.],  # initial velocity
                      dtype=np.float64)

        t_final = 10. * Tp
        t_span = (0.0, t_final)

        # integrate
        RHS = SingleModeMotionRHS(BuoyA5)
        RHS.mode = mode

        use_scipy = False
        if use_scipy:
            t_eval = np.arange(0.0, t_final, dt)
            soln = spinteg.solve_ivp(RHS, t_span=t_span, y0=x0, t_eval=t_eval,
                                     method='RK45', rtol=1e-6, atol=1e-6)
            pts_t = soln.t
            pts_x0, pts_x1 = soln.y[0, :], soln.y[1, :]

        else:
            tx = ode.Euler(RHS, x0, t_span, dt)
            t, x = list(zip(*tx))
            x = np.vstack(x)
            pts_t = np.array(t)
            pts_x0, pts_x1 = x[:, 0], x[:, 1]

        # Output
        pts_eta = np.array([Inc.eta(0.0, 0.0, t) for t in pts_t])
        pts_x = A * np.abs(Xi) * np.cos(2. * np.pi * pts_t / Tp - phase + np.angle(Xi))

        fig, ax = plt.subplots(1)

        color = 'tab:red'
        ln1 = ax.plot(pts_t, pts_x0, 'r', label='pos (td)')
        ln2 = ax.plot(pts_t, pts_x, '--g', label='pos (fd)')
        ax.set_ylabel('meters', color=color)
        ax.set_xlabel('t (s)')
        ax.tick_params(axis='y', labelcolor=color)

        color = 'tab:blue'
        ax2 = ax.twinx()
        ln3 = ax2.plot(pts_t, pts_eta, ':b', label='eta')
        ax2.set_ylabel('meters', color=color)
        ax2.set_xlabel('t (s)')
        ax2.tick_params(axis='y', labelcolor=color)

        lns = ln1 + ln2 + ln3
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, bbox_to_anchor=(1.12, 1), loc='upper left', borderaxespad=0)

        fig.suptitle(f'{modes[mode]} Motion Output')
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
