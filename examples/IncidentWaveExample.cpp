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
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>


#include "config.h"
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
/*
    std::vector<double> freq;
    std::vector<double> S;
    double grav = 9.81;
    double w0 = sqrt(.21 * grav / (2*A));
    double a = 0.0081;
    double b = 0.74;
    int n_spectrum = 10;
    double d_omega = MAX_FREQ * 2 * M_PI / n_spectrum;

    for (int i = 0; i < n_spectrum; i++) {
      double w = d_omega* (i + 1);
      S.push_back(2*M_PI*(a * grav * grav / pow(w, 5)) * exp(-b * pow(w0 / w, 4)));
      freq.push_back(w/(2*M_PI));
      }
 */
    std::vector<double> freq{0.0, 0.2/(2*M_PI), 0.4/(2*M_PI), 0.6/(2*M_PI), 2.0/(2*M_PI)};
    std::vector<double> S{0.0, 0.4*2*M_PI, 1.0*2*M_PI, 1.0*2*M_PI, 0.0}; 
    Inc.SetToCustomSpectrum(freq,S,beta);
    //T = 2*M_PI*sqrt((2*A)/grav)/0.4019;
    break;
}

  double k = pow(2 * M_PI / T, 2) / 9.81;
  std::cout << Inc << std::endl;

  {
    std::vector<double> pts_t;
    std::vector<double> pts_eta, pts_eta_true;
    double x = 0;
    double y = 0;
    double xx = x * cos(beta) + y * sin(beta);
    for (double t = 0; t < 10 * T; t += .1) {
      pts_t.push_back(t);
      pts_eta.push_back(Inc.eta(x, y, t));
      pts_eta_true.push_back(A * cos(k * xx - 2 * M_PI * t / T + phase));
    }
    Gnuplot gp;
    gp << "set term qt title  'Incident Wave Elevation at Origin'\n";
    gp << "set grid\n";
    gp << "set xlabel 'time (s)'\n";
    gp << "set ylabel '(m)'\n";
    gp << "plot '-' w l title 'eta'"
       << ",'-' w l title 'eta\\_true'\n";
    gp.send1d(boost::make_tuple(pts_t, pts_eta));
    gp.send1d(boost::make_tuple(pts_t, pts_eta_true));
  }

  double dt = .05 * T;
  double y = 0;
  for (double t = 0; t <= 3 * dt; t += dt) {
    std::vector<double> pts_x;
    std::vector<double> pts_eta, pts_eta_true;
    std::vector<double> pts_deta_dx, pts_deta_dy;
    for (double x = -1.5 * 2 * M_PI / k; x < 1.5 * 2 * M_PI / k; x += 0.5) {
      double eta;
      double detadx;
      double detady;
      eta = Inc.eta(x,y,t,&detadx,&detady);
      pts_x.push_back(x);
      pts_eta.push_back(eta);
      pts_deta_dx.push_back(detadx);
      pts_deta_dy.push_back(detady);

      double xx = x * cos(beta) + y * sin(beta);
      pts_eta_true.push_back(A * cos(k * xx - 2 * M_PI * t / T + phase));
    }
    char time[10];
    snprintf(time, sizeof(time), "%.2f", t);
    Gnuplot gp;
    gp << "set term qt title  'Incident Wave Elevation at t = " << time
       << " s'\n";
    gp << "set grid\n";
    gp << "set xlabel '(m))'\n";
    gp << "set ylabel '(m)'\n";
    gp << "plot '-' w l title 'eta'"
       << ",'-' w l title 'eta\\_true'"
       << ",'-' w l title 'deta/dx'\n";
    gp.send1d(boost::make_tuple(pts_x, pts_eta));
    gp.send1d(boost::make_tuple(pts_x, pts_eta_true));
    gp.send1d(boost::make_tuple(pts_x, pts_deta_dx));
  }
  std::cout << "Enter Ctrl-C to quit.  (Enter 'pkill gnuplot_qt' to clear "
               "plots if necessary)"
            << std::endl;
  while (1);
}
