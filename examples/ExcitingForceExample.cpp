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
#include <signal.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <memory>

#include "config.h"
#include "FS_Hydrodynamics.hpp"
#include "LinearIncidentWave.hpp"

void signal_callback_handler(int signum) {
  std::string s = "pkill gnuplot_qt";
  int ret = system(s.c_str());
  exit(signum);
}

int main(int argc, char **argv) {
  {
    std::string s = "pkill gnuplot_qt";
    int ret = system(s.c_str());
  }
  signal(SIGINT, signal_callback_handler);


  // Defaults
  double A = 1.0;
  double Tp = 5.0;
  double phase = 40.0 * M_PI / 180;
  double beta = 180.0 * M_PI / 180;

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
      std::cout << "Usage: ExcitingForceExample [-atpbh]" << std::endl;
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


  double tf = 3.0 * Tp;
  double omega = 2.0 * M_PI / Tp;
  double dt = 0.015;

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
      "./example_hydrodynamic_coeffs/BuoyA5";
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


  for (int j = 0; j < 6; j++) { // j denotes direction of resulting force
    std::complex<double> Chi =
        BuoyA5.WaveExcitingForceComponents(Inc->m_omega[0], j);
    std::vector<double> pts_F_TD, pts_F_FD, pts_eta;
    BuoyA5.SetTimestepSize(dt);
    for (int k = 0; k < pts_t.size(); k++) {
      pts_F_FD.push_back(
          Inc->m_A[0] * Chi.real() *
              cos(omega * pts_t[k] + phase * cos(180 * M_PI / 180)) -
          Inc->m_A[0] * Chi.imag() *
              sin(omega * pts_t[k] + phase * cos(180 * M_PI / 180)));
      pts_eta.push_back(Inc->eta(0, 0, pts_t[k]));
      Eigen::VectorXd ExtForce(6);
      ExtForce = BuoyA5.ExcitingForce();
      pts_F_TD.push_back(ExtForce(j));
    }
    Gnuplot gp;
    gp << "set term qt title  '" << modes[j] << " Exciting Forces'\n";
    gp << "set grid\n";
    gp << "set xlabel 'time (s)'\n";
    if (j < 3) {
      gp << "set ylabel 'F (N)'\n";
    } else {
      gp << "set  ylabel 'M (N-m)'\n";
    }
    gp << "plot '-' w l title 'Time-Domain'"
       << ",'-' w l title 'Freq-Domain'"
       << ",'-' w l title 'eta(t)'\n";
    gp.send1d(boost::make_tuple(pts_t, pts_F_TD));
    gp.send1d(boost::make_tuple(pts_t, pts_F_FD));
    gp.send1d(boost::make_tuple(pts_t, pts_eta));
    gp << "set xlabel 'time (s)'\n";
    if (j < 3) {
      gp << "set ylabel 'F (N)'\n";
    } else {
      gp << "set  ylabel 'M (N-m)'\n";
    }
    gp << "replot\n";
  }

  std::cout << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear "
               "plots if necessary)"
            << std::endl;
  while (1);
}
