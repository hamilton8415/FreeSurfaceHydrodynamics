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

  const char *modes[6] = {"Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"};
  LinearIncidentWave Inc;
  LinearIncidentWave &IncRef = Inc;
  double rho = 1025;
  double g = 9.81;
  double buoy_mass = 1400; // kg
  FS_HydroDynamics BuoyA5(IncRef, 1.0, g, rho);

  double dt = 0.005;

  int c;
  while ((c = getopt(argc, argv, "h")) != -1) {
    switch (c) {
    case 'h':
      std::cout << "Version " << BuoyA5.Version() << std::endl;
      std::cout << "Usage: PlotCoeffsExample" << std::endl;
      return 0;
      break;
    }
  }

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

  BuoyA5.Plot_FD_Coeffs();
  BuoyA5.Plot_TD_Coeffs();

  std::cout << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear "
               "plots if necessary)"
            << std::endl;
  while (1)
    ;
}
