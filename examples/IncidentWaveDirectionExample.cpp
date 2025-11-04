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
#include <cstdlib>
#include <fstream>
#include <gnuplot-iostream.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <numeric>
#include <boost/tuple/tuple.hpp>

#include <FreeSurfaceHydrodynamics/config.h>
#include <FreeSurfaceHydrodynamics/LinearIncidentWave.hpp>


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

  WaveSpectrumType SpectrumType = WaveSpectrumType::MonoChromatic;
  double A = 1;
  double T = 12;
  double phase = 0 * M_PI / 180.0;
  double beta = 180 * M_PI / 180.0;
  unsigned int seed = 0;

  int c;
  while ((c = getopt(argc, argv, "ha:t:p:b:s:S:")) != -1) {
    switch (c) {
    case 'a':
      A = atof(optarg);
      break;
    case 't':
      T = atof(optarg);
      break;
    case 'p':
      phase = atof(optarg)*M_PI/180.0;
      break;
    case 'b':
      beta = atof(optarg)*M_PI/180.0;
      break;
    case 's':
      seed = atof(optarg);
      break;
    case 'S':
      if(*optarg == 'P')
        SpectrumType = WaveSpectrumType::PiersonMoskowitz;
      if(*optarg == 'B')
        SpectrumType = WaveSpectrumType::Bretschneider;
      if(*optarg == 'C')
        SpectrumType = WaveSpectrumType::Custom;
      break;
    case 'h':
      std::cout << "Version: " << PROJECT_VER << std::endl;
      std::cout << "Usage: IncidentWaveExample [-atpbch]" << std::endl;
      std::cout << " For example:" << std::endl;
      std::cout << "  [-a 2.0] sets the incident wave amplitude to 2.0 meters" <<std::endl;
      std::cout << "  [-t -6.0] sets the monochromatic incident wave period to 6.0 seconds" <<std::endl;
      std::cout << "  [-t 8.0] sets the peak period of an Pierson Moskowitz incident wave to 8.0 seconds" <<std::endl;
      std::cout << "  [-p 45.0] sets the incident wave phase angle to 45.0 degrees" <<std::endl;
      std::cout << "  [-b 30.0] sets the incident wave direction to 30.0 degrees" <<std::endl;
      std::cout << "  [-s 30.0] sets the random seed to 42" <<std::endl;
      std::cout << "  [-S char] Sets Spectrum, 'M' = MonoChromatic, 'P' = PiersonMoskwitz, 'B' = Bretschneider, 'C' = Custom" <<std::endl;
      return 0;
      break;
    }
  }

LinearIncidentWave Inc;

if(seed > 0)
  Inc.SetSeed(seed);

switch (SpectrumType) {
  case WaveSpectrumType::MonoChromatic :
      Inc.SetToMonoChromatic(A, T, phase, beta);
    break;
  case WaveSpectrumType::PiersonMoskowitz :
    Inc.SetToPiersonMoskowitzSpectrum(2*A,beta);
    T = Inc.m_Tp;
    break;
  case WaveSpectrumType::Bretschneider :
      Inc.SetToBretschneiderSpectrum(2*A,T,beta);
    break;
  case WaveSpectrumType::Custom :
    std::vector<double> freq{ 0.029, 0.0341688345, 0.0390, 0.043932892, 0.0488, 0.0536951298, 0.0585, 0.0634555457, 0.0683, 0.0732257700, 0.0781, 0.0829861181, 0.0878, 0.0927488175, 0.0976, 0.102510986, 0.1074, 0.112273223, 0.1171, 0.122035460, 0.1269, 0.131797630, 0.1367, 0.141560329, 0.1464, 0.151320677, 0.1562, 0.161090902, 0.1660, 0.170851250, 0.1757, 0.180613949, 0.1855, 0.190376118, 0.1953, 0.200138288, 0.2050, 0.209900987, 0.2148, 0.219661335, 0.2246, 0.2294315, 0.2343, 0.239191908, 0.2441, 0.248954607, 0.2539, 0.258716776, 0.2636, 0.268479013, 0.2734, 0.278241250, 0.283, 0.288003420, 0.2929, 0.297766119, 0.3027, 0.307460319, 0.312, 0.317483218, 0.3222, 0.326489863, 0.3320, 0.340055272, 0.3515, 0.3655377, 0.3808, 0.395586045, 0.4101, 0.42375818, 0.4394, 0.457258768, 0.4687, 0.473450624, 0.4980, 0.562050220, 0.654, 0.7461795179};
    std::vector<double> S{0.000999424000, 0.00499205992, 0.00900300, 0.0219735333, 0.03500851, 0.0524594296, 0.07001702, 0.0605521481, 0.05101158, 0.0629676605, 0.0750182, 0.172313322, 0.27056537, 0.65807237, 1.04975564, 1.42013208, 1.7949388, 1.70947731, 1.6228966, 1.27950984, 0.93122764, 0.836434556, 0.74017996, 0.642456728, 0.54313164, 0.469537840, 0.39459635, 0.372046441, 0.34908569, 0.306504292, 0.26306355, 0.222481950, 0.181043, 0.164968711, 0.14853529, 0.138648986, 0.12853043, 0.118648956, 0.10852556, 0.09717288, 0.08552038, 0.0941588507, 0.1030246, 0.0973540982, 0.09152307, 0.0833895727, 0.0750182, 0.067875086, 0.06051430, 0.0592834484, 0.05801369, 0.0557992646, 0.05351219, 0.0503162942, 0.0470118, 0.0457830270, 0.04451123, 0.0418478599, 0.03901030, 0.0341634123, 0.02950758, 0.02950758, 0.02950758, 0.0276578248, 0.0250060, 0.0228586044, 0.02050457, 0.02050457, 0.02050457, 0.0198085867, 0.0190054, 0.0177892797, 0.01700454, 0.0162021945, 0.01200332, 0.00852036595, 0.00350003, 0.001500151916};
    Inc.SetToCustomSpectrum(freq,S,beta);
    break;
}

  double k = pow(2 * M_PI / T, 2) / 9.81;
  std::cout << Inc << std::endl;

  double dt = .05 * T;
  double y = 0;
  for (double t = 0; t <= 4 * dt; t += dt) {
    char time[10];
    snprintf(time, sizeof(time), "%.2f", t);

    char outputfilename[64];
    sprintf(outputfilename, "elev_t%.2f.dat", t);

    char gnuplot_cmd[1024];
    sprintf(gnuplot_cmd,"splot '%s' w l,'%s' with vectors", outputfilename,outputfilename);

    std::ofstream outputFile(outputfilename);

    double n_wavelengths = 1;
    int n_divisions = 10;

    for (double x = -(n_wavelengths/2) * 2 * M_PI / k; x <= (n_wavelengths/2) * 2 * M_PI / k; x += (n_wavelengths/n_divisions)*2*M_PI/k) {
     for (double y = -(n_wavelengths/2) * 2 * M_PI / k; y <= (n_wavelengths/2) * 2 * M_PI / k; y += (n_wavelengths/n_divisions)*2*M_PI/k) {
      double detadx,detady;
      double u_east,v_north;
      double eta = Inc.eta(x,y,t,nullptr,&detady,&u_east,&v_north);
      double eta_dot = Inc.etadot(x,y,t);
      outputFile << x << " " << y << " " << eta << " " << 100*u_east << " " << 100*v_north << " " << eta_dot << std::endl;
     }
     outputFile << std::endl;
    }
    outputFile.close();
    Gnuplot gp;
    //gp << "set term qt title  'Incident Wave Elevation at t = " << time
    //   << " s'\n";
    gp << "set term qt title  'Incident Wave Elevation at t = " << time
       << " s'\n";
    gp << "set grid\n";
    gp << "set xlabel '(m))'\n";
    gp << "set ylabel '(m)'\n";
    gp << "set contour base\n";
    gp << gnuplot_cmd << std::endl;
    std::cout << "For interactive plot, use \""  << gnuplot_cmd << "\" in gnuplot from command line" << std::endl;
  }
  std::cout << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear "
               "plots if necessary)"
            << std::endl;
  while (1);
}
