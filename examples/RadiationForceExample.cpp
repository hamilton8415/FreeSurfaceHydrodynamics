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
  double Tp = 5;
  double phase = 40 * M_PI / 180;

  int c;
  while ((c = getopt(argc, argv, ":a:t:p:h")) != -1) {
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
    case 'h':
      std::cout << "Version: " << PROJECT_VER << std::endl;
      std::cout << "Usage: RadiationForceExample [-atph]" << std::endl;
      std::cout << " For example:" << std::endl;
      std::cout << "  [-a 2.0] sets the body motoin amplitude to 2.0 meters" <<std::endl;
      std::cout << "  [-t 6.0] sets the body motion period to 6.0 seconds" <<std::endl;
      std::cout << "  [-p 45.0] sets the body motion phase angle to 45.0 degrees" <<std::endl;
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

  double omega = 2.0 * M_PI / Tp;
  double tf = 3.0 * Tp;
  double dt = 0.01;

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

  // Test Radiation Forces
  // Note:  This computes the radiation forces for each mode of motion
  // indivdually, so MemRadiation is called 6 times per timestep, once with each
  // acceleration set non-zero In use it will be called only once per time-step,
  // with all acclerations set.
  for (int i = 0; i < 6; i++) {   // i determines mode of motion.
    for (int j = 0; j < 6; j++) { // j denotes direction of resulting force
      double am = BuoyA5.AddedMass(omega, i, j);
      double dmp = BuoyA5.RadiationDamping(omega, i, j);
      double am_inf = BuoyA5.fd_A_inf_freq(i, j);

      std::vector<double> pts_F_TD, pts_F_FD;
      double last_accel = 0;
      double F_max = -std::numeric_limits<double>::max();
      double F_min = std::numeric_limits<double>::max();
      BuoyA5.SetTimestepSize(
          dt); // Reset timestep to re-initialize storage in BuoyA5 Class
      for (int k = 0; k < pts_accel.size(); k++) {
        double accel = pts_accel[k];
        Eigen::VectorXd xddot(6);
        for (int n = 0; n < 6; n++) {
          xddot(n) = 0;
        }
        xddot(i) = last_accel;

        Eigen::VectorXd MemForce(6);
        MemForce = BuoyA5.RadiationForce(xddot);
        pts_F_TD.push_back(am_inf * accel + MemForce(j));
        last_accel = accel;
        double FD_Force = am * pts_accel[k] + dmp * pts_vel[k];
        pts_F_FD.push_back(FD_Force);
        if (FD_Force > F_max) {
          F_max = FD_Force;
        }
        if (FD_Force < F_min) {
          F_min = FD_Force;
        }
      }

      if ((F_min < -1) && (F_max > 1)) { // Don't plot near-zero forces
        Gnuplot gp;
        char Amp[10];
        snprintf(Amp, sizeof(Amp), "%.1f", A);
        char Per[10];
        snprintf(Per, sizeof(Per), "%.1f", Tp);
        gp << "set term qt title 'Radiation Forces: " << modes[i]
           << " Motions, " << modes[j] << " Forces:  A = " << Amp
           << "m  T = " << Per << "s  \n";
        gp << "set grid\n";
        gp << "set xlabel 'time (s)'\n";
        if (j < 3) {
          gp << "set ylabel 'F (N)'\n";
        } else {
          gp << "set  ylabel 'M (N-m)'\n";
        }
        gp << "plot '-' w l title 'Time-Domain'"
           << ",'-' w l title 'Freq-Domain'\n";
        gp.send1d(boost::make_tuple(pts_t, pts_F_TD));
        gp.send1d(boost::make_tuple(pts_t, pts_F_FD));
        gp << "set xlabel 'time (s)'\n";
        if (j < 3) {
          gp << "set ylabel 'F (N)'\n";
        } else {
          gp << "set  ylabel 'M (N-m)'\n";
        }
        gp << "set title '" << modes[i] << "(t) = " << std::fixed
           << std::setprecision(1) << Amp << "cos(2 pi t/" << Per << ")'  \n";
        gp << "replot\n";
      }
    }
  }

  std::cout
      << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear plots "
         "if necessary)"
      << std::endl;
  while (1)
    ;
}
