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

  const char *modes[6] = {"Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"};
  double rho = 1025;
  double g = 9.81;
  double buoy_mass = 1400; // kg

  // Defaults
  double A = 1;  // .5 + ((float)(std::rand() % 20) / 10);
  double Tp = 5; // 3.0 + (std::rand() % 9);
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
      std::cout << "Usage: GravityForceExample [-ah]" << std::endl;
      std::cout << " For example:" << std::endl;
      std::cout << "  [-a 2.0] sets the body motion amplitude to 2.0 meters" <<std::endl;
      return 0;
      break;
    }
  }

  double tf = 2.0 * Tp;
  double dt = 0.01;

  FS_HydroDynamics BuoyA5;

  BuoyA5.SetCOG(0, 0,
                -.24); // Set COG relative to waterplane coordinate system.
                       //  BuoyA5.SetVolume(buoy_mass / rho);
  BuoyA5.SetMass(buoy_mass);

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

  std::cout << "Computing Gravity Forces" << std::endl;
  for (int j = 0; j < 6; j++) { // j denotes direction of resulting force
    std::vector<double> pts_x, pts_F_G;
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
      Eigen::VectorXd GravityForce(6);
      GravityForce = BuoyA5.GravityForce(x);
      pts_F_G.push_back(GravityForce(j));
    }
    Gnuplot gp;
    gp << "set term qt title  '" << modes[j] << " Gravity Forces'\n";
    gp << "set grid\n";
    gp << "set xlabel 'time (s)'\n";
    if (j < 3) {
      gp << "set ylabel 'F (N)'\n";
    } else {
      gp << "set  ylabel 'M (N-m)'\n";
    }
    gp << "plot '-' w l title 'x(t)'"
       << ",'-' w l title 'Gravity Force'\n";
    gp.send1d(boost::make_tuple(pts_t, pts_x));
    gp.send1d(boost::make_tuple(pts_t, pts_F_G));
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
  while (1)
    ;
}
