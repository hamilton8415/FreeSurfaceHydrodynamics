#!/usr/bin/python3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fshd import LinearIncidentWave, WaveSpectrumType


_wave_spectrum_type = dict(
    M=WaveSpectrumType.MonoChromatic,
    P=WaveSpectrumType.PiersonMoskowitz,
    B=WaveSpectrumType.Bretschneider,
    C=WaveSpectrumType.Custom
)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=float, default=1.0,
                        help='(meters) sets the incident wave amplitude')
    parser.add_argument('-t', type=float, default=12.0,
                        help='(seconds) sets the incident wave period')
    parser.add_argument('-p', type=float, default=0.0,
                        help='(degrees) sets the incident wave phase angle')
    parser.add_argument('-b', type=float, default=180.0,
                        help='(degrees) sets the incident wave direction')
    parser.add_argument('-s', type=int, default=0,
                        help='sets the random seed')
    parser.add_argument('-S', type=str, default='M',
                        choices=['M', 'P', 'B', 'C'],
                        help="Sets Spectrum:" \
                                 + " 'M'onoChromatic" \
                                 + ", 'P'iersonMoskowitz" \
                                 + ", 'B'retschneider" \
                                 + ", 'C'ustom")
    args = parser.parse_args()

    SpectrumType = _wave_spectrum_type[args.S]
    A = args.a
    T = args.t
    phase = args.p * np.pi / 180.0
    beta = args.b * np.pi / 180.0
    seed = args.s


    Inc = LinearIncidentWave()

    if seed > 0:
        Inc.SetSeed(seed)

    if SpectrumType == WaveSpectrumType.MonoChromatic:
        Inc.SetToMonoChromatic(A, T, phase, beta)
    elif SpectrumType == WaveSpectrumType.PiersonMoskowitz:
        Inc.SetToPiersonMoskowitzSpectrum(2.*A, beta)
        T = Inc.m_Tp
    elif SpectrumType == WaveSpectrumType.Bretschneider:
        Inc.SetToBretschneiderSpectrum(2.*A, T, beta)
    elif SpectrumType == WaveSpectrumType.Custom:
        grav = 9.81
        w0 = sqrt(0.21 * grav / (2. * A))
        a = 0.0081
        b = 0.74
        n_spectrum = 100
        d_omega = MAX_FREQ * 2 * M_PI / n_spectrum

        omega = d_omega * np.arange(1, n_spectrum + 1)
        S = (a * grav**2. / omega**5.) * exp(-b * (w0 / w)**4.)
        Inc.SetToCustomSpectrum(omega, S, beta)
        T = 2. * np.pi * np.sqrt((2. * A) / grav) / 0.4019

    k = ((2. * np.pi / T)**2.) / 9.81
    print(Inc)

    # plot
    pts_t = []
    pts_eta = []
    pts_eta_true = []
    x = 0.
    y = 0.
    xx = x * np.cos(beta) + y * np.sin(beta)
    for t in np.arange(0., 10.*T, 0.1):
        pts_t.append(t)
        eta = Inc.eta(x, y, t)
        pts_eta.append(eta)
        pts_eta_true.append(A * np.cos(k * xx - 2. * np.pi * t / T + phase))


    fig, ax = plt.subplots(1)

    ax.plot(pts_t, pts_eta, 'r', label='eta(t)')
    ax.plot(pts_t, pts_eta_true, '--g', label='eta_true(t)')
    ax.set_ylabel('meters')
    ax.set_xlabel('time (s)')
    ax.legend()

    fig.suptitle(f'Incident Wave Elevation at Origin')
    fig.tight_layout()

    dt = 0.05 * T
    y = 0.
    for t in np.arange(0., 4.0*dt, dt):
        pts_x = []
        pts_eta = []
        pts_eta_true = []
        pts_deta_dx = []
        pts_deta_dy = []
        for x in np.arange(-1.5 * 2. * np.pi / k,
                           1.5 * 2. * np.pi / k,
                           2.0):
            eta, deta_dx, deta_dy = Inc.eta(x, y, t, compute_deta=True)
            pts_x.append(x)
            pts_eta.append(eta)
            pts_deta_dx.append(deta_dx)
            pts_deta_dy.append(deta_dy)

            xx = x * np.cos(beta) + y * np.sin(beta)
            pts_eta_true.append(A * np.cos(k * xx - 2. * np.pi * t / T + phase))

        fig, ax = plt.subplots(1)

        color = 'tab:red'
        ln1 = ax.plot(pts_x, pts_eta, 'r', label='eta (m)')
        ax.set_ylabel('meters', color=color)
        ax.tick_params(axis='y', labelcolor=color)
        ax.set_xlabel('meters')

        ln2 = ax.plot(pts_x, pts_eta_true, '--g', label='eta true (m)')

        color = 'tab:blue'
        ax2 = ax.twinx()
        ln3 = ax2.plot(pts_x, pts_deta_dx, 'b', label='deta/dx (m)')
        ax2.set_ylabel('meters', color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        lns = ln1 + ln2 + ln3
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc=1)

        fig.suptitle(f'Incident Wave Elevation at {t = :.2f} (s)')
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
























