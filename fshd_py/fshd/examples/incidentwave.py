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
    parser.add_argument('-E', default=None,
                        help='specify json with S and f (spectrum and freq); defaults otherwise')
    parser.add_argument('-s', type=int, default=0,
                        help='sets the random seed')
    parser.add_argument('-S', type=str, default='M',
                        choices=['M', 'P', 'B', 'C'],
                        help='Sets Spectrum:' \
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
        # grav = 9.81
        # w0 = sqrt(0.21 * grav / (2. * A))
        # a = 0.0081
        # b = 0.74
        # n_spectrum = 100
        # d_omega = MAX_FREQ * 2 * M_PI / n_spectrum

        # omega = d_omega * np.arange(1, n_spectrum + 1)
        # S = (a * grav**2. / omega**5.) * exp(-b * (w0 / w)**4.)
        # Inc.SetToCustomSpectrum(omega, S, beta)
        # T = 2. * np.pi * np.sqrt((2. * A) / grav) / 0.4019
        freq = None
        S = None
        if args.E is not None:
            import json
            import pandas as pd
            with open(args.E, 'r') as fd:
                data = json.load(fd)
            dat = pd.DataFrame(data)
            try:
                df = dat.df.to_numpy()
                freq = dat.frequency.to_numpy()
                S = dat.varianceDensity.to_numpy()
            except AttributeError as err:
                print('Error: JSON file must have keys [df, frequency, varianceDensity];'
                      ' defaulting')

        if (freq is None) or (S is None):
            freq = np.array([0.029, 0.0341688345, 0.0390, 0.043932892, 0.0488, 0.0536951298,
                             0.0585, 0.0634555457, 0.0683, 0.0732257700, 0.0781, 0.0829861181,
                             0.0878, 0.0927488175, 0.0976, 0.102510986, 0.1074, 0.112273223,
                             0.1171, 0.122035460, 0.1269, 0.131797630, 0.1367, 0.141560329,
                             0.1464, 0.151320677, 0.1562, 0.161090902, 0.1660, 0.170851250,
                             0.1757, 0.180613949, 0.1855, 0.190376118, 0.1953, 0.200138288,
                             0.2050, 0.209900987, 0.2148, 0.219661335, 0.2246, 0.2294315, 0.2343,
                             0.239191908, 0.2441, 0.248954607, 0.2539, 0.258716776, 0.2636,
                             0.268479013, 0.2734, 0.278241250, 0.283, 0.288003420, 0.2929,
                             0.297766119, 0.3027, 0.307460319, 0.312, 0.317483218, 0.3222,
                             0.326489863, 0.3320, 0.340055272, 0.3515, 0.3655377, 0.3808,
                             0.395586045, 0.4101, 0.42375818, 0.4394, 0.457258768, 0.4687,
                             0.473450624, 0.4980, 0.562050220, 0.654, 0.7461795179])
            df = np.append([freq[1]-freq[0]], np.diff(freq))
            S = np.array([0.000999424000, 0.00499205992, 0.00900300, 0.0219735333, 0.03500851,
                          0.0524594296, 0.07001702, 0.0605521481, 0.05101158, 0.0629676605,
                          0.0750182, 0.172313322, 0.27056537, 0.65807237, 1.04975564, 1.42013208,
                          1.7949388, 1.70947731, 1.6228966, 1.27950984, 0.93122764, 0.836434556,
                          0.74017996, 0.642456728, 0.54313164, 0.469537840, 0.39459635,
                          0.372046441, 0.34908569, 0.306504292, 0.26306355, 0.222481950,
                          0.181043, 0.164968711, 0.14853529, 0.138648986, 0.12853043,
                          0.118648956, 0.10852556, 0.09717288, 0.08552038, 0.0941588507,
                          0.1030246, 0.0973540982, 0.09152307, 0.0833895727, 0.0750182,
                          0.067875086, 0.06051430, 0.0592834484, 0.05801369, 0.0557992646,
                          0.05351219, 0.0503162942, 0.0470118, 0.0457830270, 0.04451123,
                          0.0418478599, 0.03901030, 0.0341634123, 0.02950758, 0.02950758,
                          0.02950758, 0.0276578248, 0.0250060, 0.0228586044, 0.02050457,
                          0.02050457, 0.02050457, 0.0198085867, 0.0190054, 0.0177892797,
                          0.01700454, 0.0162021945, 0.01200332, 0.00852036595, 0.00350003,
                          0.001500151916])

        Inc.SetToCustomSpectrum(freq, S, beta)
        m0 = S @ (freq**0. * df)
        m1 = S @ (freq**1. * df)
        Hm0 = 4.*np.sqrt(m0)
        A = Hm0/2.
        T = m0/m1
        print(f'{m0=}\n{m1=}\n{Hm0=}\n{A=}\n{T=}')

    k = ((2. * np.pi / T)**2.) / 9.81
    # print(Inc)

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
