// Copyright 2022 Monterey Bay Aquarium Research Institute
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <gnuplot-iostream.h>
#include <iostream>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>
#include <csignal>

#include <FreeSurfaceHydrodynamics/config.h>
#include <FreeSurfaceHydrodynamics/FS_Hydrodynamics.hpp>
#include <FreeSurfaceHydrodynamics/LinearIncidentWave.hpp>


void signal_callback_handler(int signum) {
   std::string s = "pkill gnuplot_qt";
   int ret = system(s.c_str());
   exit(signum);
}

/* The rhs of x' = f(x) defined as a class */
class SingleModeMotionRHS {
public:
  FS_HydroDynamics *FloatingBody = NULL;
  double last_accel = 0;
  int mode = 0;

  explicit SingleModeMotionRHS(FS_HydroDynamics *Body) : FloatingBody(Body) {}

  // x[0] = position
  // x[1] = velocity
  void operator()(const std::vector<double> &x, std::vector<double> &dxdt,
                  const double t) {
    Eigen::VectorXd pos(6);
    pos(0) = 0; pos(1) = 0; pos(2) = 0; pos(3) = 0; pos(4) = 0; pos(5) = 0;
    pos(mode) = x[0];
    Eigen::VectorXd vel(6);
    vel(0) = 0; vel(1) = 0; vel(2) = 0; vel(3) = 0; vel(4) = 0; vel(5) = 0;
    vel(mode) = x[1];
    Eigen::VectorXd F_LinDamping(6);
    F_LinDamping = FloatingBody->LinearDampingForce(vel);
    Eigen::VectorXd F_B(6);
    F_B = FloatingBody->BuoyancyForce(pos);
    Eigen::VectorXd F_G(6);
    F_G = FloatingBody->GravityForce(pos);
    Eigen::VectorXd F_R(6);
    Eigen::VectorXd accel(6);
    accel(0) = 0;
    accel(1) = 0;
    accel(2) = 0;
    accel(3) = 0;
    accel(4) = 0;
    accel(5) = 0;
    accel(mode) = last_accel;
    F_R = -FloatingBody->RadiationForce(accel);
    Eigen::VectorXd F_E(6);
    F_E = FloatingBody->ExcitingForce();
    dxdt[0] = x[1];
    double b = 0.1;
    dxdt[1] =
        (F_LinDamping(mode) + F_B(mode) + F_G(mode) + F_R(mode) + F_E(mode)) /
        (FloatingBody->M(mode, mode) +
         FloatingBody->AddedMass(10000.0, mode, mode));
    last_accel = dxdt[1];
  }
};

//[ integrate_observer
struct push_back_state_and_time {
  std::vector<std::vector<double>> &m_states;
  std::vector<double> &m_times;

  push_back_state_and_time(std::vector<std::vector<double>> &states,
                           std::vector<double> &times)
      : m_states(states), m_times(times) {}

  void operator()(const std::vector<double> &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

int main(int argc, char **argv) {
  {
    std::string s = "pkill gnuplot_qt";
    int ret = system(s.c_str());
  }
  signal(SIGINT, signal_callback_handler);

  // Defaults
  double A = 1.0;
  double Tp = 5.0;
  double phase = 40.0* M_PI / 180.0;
  double beta = 20.0 * M_PI /180.0;

  int c;
  while ((c = getopt(argc, argv, ":a:t:p:b:h")) != -1) {
    switch (c) {
    case 'a':
      A = atof(optarg);
      break;
    case 't':
      Tp = atof(optarg);
      break;
    case 'p':
      phase = atof(optarg)*M_PI/180.0;
      break;
    case 'b':
      beta = atof(optarg)*M_PI/180.0;
      break;
    case 'h':
      std::cout << "Version: " << PROJECT_VER << std::endl;
      std::cout << "Usage: MotionExample [-atpbh]" << std::endl;
      std::cout << " For example:" << std::endl;
      std::cout << "  [-a 2.0] sets the incident wave amplitude to 2.0 meters" <<std::endl;
      std::cout << "  [-t 6.0] sets the incident wave period to 6.0 seconds" <<std::endl;
      std::cout << "  [-p 45.0] sets the incident wave phase angle to 45.0 degrees" <<std::endl;
      std::cout << "  [-b 30.0] sets the incident wave direction to 30.0 degrees" <<std::endl;
      return 0;
      break;
    }
  }

  const char *modes[6] = {"Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"};

  std::shared_ptr<LinearIncidentWave> Inc = std::make_shared<LinearIncidentWave> ();

  double rho = 1025;
  double g = 9.81;
  double buoy_mass = 1400; // kg
  FS_HydroDynamics BuoyA5;

  double omega = 2 * M_PI / Tp;
  double tf = 2.0 * Tp;
  double dt = 0.005;

  Inc->SetToMonoChromatic(A, Tp, phase, beta);

  BuoyA5.AssignIncidentWave(Inc);

  BuoyA5.SetWaterplane(5.47, 1.37,
                       1.37); // Set area and 2nd moments of area for waterplane
  BuoyA5.SetCOB(0, 0,
                -.22); // Set COB relative to waterplane coordinate system.
  BuoyA5.SetCOG(0, 0,
                -.24); // Set COG relative to waterplane coordinate system.
  BuoyA5.SetVolume(buoy_mass / rho);
  BuoyA5.SetMass(buoy_mass);

  std::string HydrodynamicsBaseFilename =
      "./example_hydrodynamic_coeffs/mbari_snl";
//      "./example_hydrodynamic_coeffs/BuoyA5";
  BuoyA5.ReadWAMITData_FD(HydrodynamicsBaseFilename);
  BuoyA5.ReadWAMITData_TD(HydrodynamicsBaseFilename);
  BuoyA5.SetTimestepSize(dt);

  std::vector<double> pts_t;
  std::vector<double> pts_pos;
  std::vector<double> pts_vel;
  std::vector<double> pts_accel;
  for (int k = 0; k < tf / dt; k++) {
    double tt = dt * k;
    pts_t.push_back(tt);
    pts_pos.push_back(A * cos(omega * tt + phase));
    pts_vel.push_back(-A * omega * sin(omega * tt + phase));
    pts_accel.push_back(-A * pow(omega, 2) * cos(omega * tt + phase));
  }

    Eigen::Matrix<double, 3, 3> I;
    I << 1500, 0, 0, 0, 1500, 0, 0, 0, 650;
    BuoyA5.SetI(I);

    Eigen::VectorXd b(6);
    b(0) = 300.0;
    b(1) = 300.0;
    b(2) = 900.0;
    b(3) = 400.0;
    b(4) = 400.0;
    b(5) = 100.0;
    BuoyA5.SetDampingCoeffs(b);

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

    BuoyA5.SetTimestepSize(dt);

    for (int mode = 0; mode < 6; mode++) {
      auto Xi = BuoyA5.ComplexAmplitude(2 * M_PI / Tp, mode);

      std::vector<double> x(2);
      x[0] = 0.0; // initial position
      x[1] = 0.0; // initial velocity

      double t_final = 10 * Tp;
      // integrate_observ
      std::vector<std::vector<double>> x_vec;
      std::vector<double> times;
      boost::numeric::odeint::euler<std::vector<double>> stepper;
      SingleModeMotionRHS RHS(&BuoyA5);
      RHS.mode = mode;
      int steps = boost::numeric::odeint::integrate_const(
          stepper, RHS, x, 0.0, t_final, dt,
          push_back_state_and_time(x_vec, times));

      /* output */
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

  std::cout << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear plots if necessary)"  << std::endl;
  while(1);

}
