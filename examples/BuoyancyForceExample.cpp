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
  double A = 1;
  double Tp = 10;
  double omega = 2 * M_PI / Tp;
  double phase = 0 * M_PI / 180;

  int c;
  while ((c = getopt(argc, argv, ":a:h")) != -1) {
    switch (c) {
    case 'a':
      A = atof(optarg);
      break;
    case 'h':
      std::cout << "Version: " << PROJECT_VER << std::endl;
      std::cout << "Usage: BuoyancyForceExample [-ah]" << std::endl;
      std::cout << " For example:" << std::endl;
      std::cout << "  [-a 2.0] sets the body motion amplitude to 2.0 meters" <<std::endl;
      return 0;
      break;
    }
  }

  const char *modes[6] = {"Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"};
  LinearIncidentWave Inc;
  LinearIncidentWave &IncRef = Inc;
  double rho = 1025;
  double g = 9.81;
  double buoy_mass = 1400; // kg
  FS_HydroDynamics BuoyA5(IncRef, 1.0, g, rho);

  double tf = 2.0 * Tp;
  double dt = 0.1;

  BuoyA5.SetWaterplane(5.47, 1.37,
                       1.37); // Set area and 2nd moments of area for waterplane
  BuoyA5.SetCOB(0, 0,
                -.22); // Set COB relative to waterplane coordinate system.
  BuoyA5.SetVolume(buoy_mass / rho);

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

  Gnuplot gp;
  char Amp[10];
  snprintf(Amp, sizeof(Amp), "%.1f", A);
  char Per[10];
  snprintf(Per, sizeof(Per), "%.1f", Tp);
  gp << "set term qt title  'A = " << Amp << "m  T = " << Per << "s'\n";
  gp << "set grid\n";
  gp << "set xlabel 'time (s)'\n";
  gp << "plot '-' w l title 'Vel'"
     << ",'-' w l title 'Accel'\n";
  gp.send1d(boost::make_tuple(pts_t, pts_vel));
  gp.send1d(boost::make_tuple(pts_t, pts_accel));

  std::cout << "Computing Buoyancy Forces" << std::endl;
  for (int j = 0; j < 6; j++) { // j denotes direction of resulting force
    std::vector<double> pts_x, pts_F_B;
    for (int k = 0; k < pts_pos.size(); k++) {
      Eigen::VectorXd x(6);
      x(0) = 0;
      x(1) = 0;
      x(2) = 0;
      x(3) = 0;
      x(4) = 0;
      x(5) = 0;
      x(j) = pts_pos[k];
      pts_x.push_back(x(j));
      Eigen::VectorXd BuoyancyForce(6);
      BuoyancyForce = BuoyA5.BuoyancyForce(x);
      pts_F_B.push_back(BuoyancyForce(j));
    }
    Gnuplot gp;
    gp << "set term qt title  '" << modes[j] << " Buoyancy Forces'\n";
    gp << "set grid\n";
    gp << "set xlabel 'time (s)'\n";
    if (j < 3) {
      gp << "set ylabel 'F (N)'\n";
    } else {
      gp << "set  ylabel 'M (N-m)'\n";
    }
    gp << "plot '-' w l title 'x(t)'"
       << ",'-' w l title 'Buoyancy Force'\n";
    gp.send1d(boost::make_tuple(pts_t, pts_x));
    gp.send1d(boost::make_tuple(pts_t, pts_F_B));
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
