#!/usr/bin/python3

import numpy as np

from fshd import LinearIncidentWave, FS_HydroDynamics


M_PI = 3.141592653589793238462643383


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', nargs=1, type=float, default=None,
                        help='(meters) sets the incident wave amplitude')
    parser.add_argument('-t', nargs=1, type=float, default=None,
                        help='(seconds) sets the incident wave period')
    parser.add_argument('-p', nargs=1, type=float, default=None,
                        help='(degrees) sets the incident wave phase angle')
    parser.add_argument('-b', nargs=1, type=float, default=None,
                        help='(degrees) sets the incident wave direction')
    args = parser.parse_args()

    # Defaults
    A = 1.0
    Tp = 5.0
    phase = 40.0 * M_PI / 180.0
    beta = 20.0 * M_PI / 180.0

    if args.a:
        A = args.a
    if args.t:
        Tp = args.t
    if arg.p:
        phase = args.p * M_PI / 180.0
    if args.b:
        beta = args.b * M_PI / 180.0

    modes = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    Inc = LinearIncidentWave()

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg
    BuoyA5 = FS_HydroDynamics()

    omega = 2 * M_PI / Tp
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

    # TODO(andermi) use importlib.resource
    HydrodynamicsBaseFilename =
      './example_hydrodynamic_coeffs/BuoyA5'
    BuoyA5.ReadWAMITData_FD(HydrodynamicsBaseFilename)
    BuoyA5.ReadWAMITData_TD(HydrodynamicsBaseFilename)
    BuoyA5.SetTimestepSize(dt)

    pts_t = []
    pts_pos = []
    pts_vel = []
    pts_accel = []
    k = 0
    while k < tf / dt:
        tt = dt * k
        k += 1
        pts_t.append(tt)
        pts_pos.append(A * cos(omega * tt + phase))
        pts_vel.append(-A * omega * sin(omega * tt + phase))
        pts_accel.append(-A * pow(omega, 2) * cos(omega * tt + phase))

    I = np.diag([1500., 1500., 650.], dtype=np.float64)
    BuoyA5.SetI(I)

    b = np.array([300., 300., 900., 400., 400., 100.], dtype=np.float64)
    BuoyA5.SetDampingCoeffs(b);

    """ # TODO(andermi) plotting

        std::vector<double> pts_T;
        Eigen::Matrix<std::vector<double>, 6, 1> pts_XiMod, pts_XiPh;

        double dT = .01;
        for (double T = dT; T < 24; T += dT) {
          pts_T.push_back(T);
          double w = 2 * M_PI / T;
          auto Xi = BuoyA5.ComplexAmplitude(w);
          for (int j = 0; j < 6; j++) {
            pts_XiMod(j).push_back(std::abs(Xi(j)));
            pts_XiPh(j).push_back(std::arg(Xi(j)) * 180 / M_PI);
          }
        }

        for (int j = 0; j < 6; j++) {
          Gnuplot gp1;
          gp1 << "set term qt title  '" << modes[j] << " RAO Amplitude'\n";
          gp1 << "set grid\n";
          gp1 << "set xlabel 'Wave Period (s)'\n";
          gp1 << "plot '-' w l \n";
          gp1.send1d(boost::make_tuple(pts_T, pts_XiMod(j)));

          Gnuplot gp2;
          gp2 << "set term qt title  '" << modes[j] << " RAO Phase Angle (deg)'\n";
          gp2 << "set grid\n";
          gp2 << "set xlabel 'Wave Period (s)'\n";
          gp2 << "plot '-' w l \n";
          gp2.send1d(boost::make_tuple(pts_T, pts_XiPh(j)));
        }

    """

    BuoyA5.SetTimestepSize(dt);

    for mode in range(len(modes)):
        Xi = BuoyA5.ComplexAmplitude(2. * M_PI / Tp, mode)

      std::vector<double> x(2);
      x[0] = 0.0;  # initial position
      x[1] = 0.0;  # initial velocity

      double t_final = 10 * Tp;
      # integrate_observ
      std::vector<std::vector<double>> x_vec;
      std::vector<double> times;
      boost::numeric::odeint::euler<std::vector<double>> stepper;
      SingleModeMotionRHS RHS(&BuoyA5);
      RHS.mode = mode;
      int steps = boost::numeric::odeint::integrate_const(
          stepper, RHS, x, 0.0, t_final, dt,
          push_back_state_and_time(x_vec, times));

      # Output
      std::vector<double> pts_t, pts_x0, pts_x1;
      std::vector<double> pts_x;
      std::vector<double> pts_eta;
      for (size_t i = 0; i <= steps; i++) {
        pts_t.push_back(times[i]);
        pts_eta.push_back(Inc->eta(0.0, 0.0, times[i]));
        pts_x0.push_back(x_vec[i][0]);
        pts_x1.push_back(x_vec[i][1]);
        pts_x.push_back(A * std::abs(Xi) *
                        cos(2 * M_PI * times[i] / Tp - phase + std::arg(Xi)));
      }

      Gnuplot gp;
      gp << "set term qt title  '" << modes[mode] << " Motion Output'\n";
      gp << "set grid\n";
      gp << "set xlabel 't (s)'\n";
      gp << "plot '-' w l title 'eta'"
         << ",'-' w l title 'pos (td)'"
         << ",'-' w l title 'pos (fd)'\n";
      gp.send1d(boost::make_tuple(pts_t, pts_eta));
      gp.send1d(boost::make_tuple(pts_t, pts_x0));
      gp.send1d(boost::make_tuple(pts_t, pts_x));
    }




def __name__=='__main__':
    main()










