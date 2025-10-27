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
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <FreeSurfaceHydrodynamics/config.h>
#include <FreeSurfaceHydrodynamics/interp1d.hpp>
#include <FreeSurfaceHydrodynamics/LinearIncidentWave.hpp>



/// \brief Constructor, defaults to monotchromatic wave and default gravity and density
LinearIncidentWave::LinearIncidentWave() : m_grav(9.81), m_rho(1025)
{
    std::srand(time(0));
}

/// \brief Constructor, sets seed and defaults to monotchromatic wave and default gravity and density
LinearIncidentWave::LinearIncidentWave(unsigned int seed) : m_grav(9.81), m_rho(1025)
{
  SetSeed(seed);
}

/// \brief Sets seed to specified value
void LinearIncidentWave::SetSeed(unsigned int seed) 
{
  if(seed == 0)
    std::srand(time(0));
  else
    std::srand(seed);
}

/// \brief Select single frequency wave
void LinearIncidentWave::SetToMonoChromatic(double A, double T, double phase, double beta)
{
  std::cout << "Mono A = " << A << std::endl;
  m_SpectrumType = WaveSpectrumType::MonoChromatic;
  m_Hs = 2 * A;
  m_Tp = T;
  m_beta = beta;
  m_omega.resize(1);
  m_k.resize(1);
  m_phases.resize(1);
  m_Spectrum.resize(1);
  m_A.resize(1);
  m_phases(0) = phase;
  m_omega(0) = 2 * M_PI / T;
  m_k(0) = m_omega(0) * m_omega(0) / m_grav;
  m_A(0) = A;
}

/// \brief Select PM-Spectrum (default num of phases)
void LinearIncidentWave::SetToBretschneiderSpectrum(double Hs, double Tp, double beta)
{
  SetToBretschneiderSpectrum(Hs, Tp, beta, DEFAULT_N_PHASES);
}

/// \brief Select Bretschneider Spectrum (set num of phases)
void LinearIncidentWave::SetToBretschneiderSpectrum(
  double Hs, double Tp, double beta,
  int n_phases)
{
  m_SpectrumType = WaveSpectrumType::Bretschneider;
  m_beta = beta;
  m_Hs = Hs;
  
  double wp = 2.0*M_PI/(1.2957*Tp);

  m_omega.resize(n_phases);
  m_k.resize(n_phases);
  m_phases.resize(n_phases);
  m_Spectrum.resize(n_phases);
  m_A.resize(n_phases);

  double d_omega = MAX_FREQ * 2 * M_PI / n_phases;

  for (int i = 0; i < m_k.size(); i++) {
    m_omega(i) = d_omega * (i + 1) + (0.25*d_omega*(std::rand()-RAND_MAX/2))/(RAND_MAX/2);
    m_k(i) = m_omega(i) * m_omega(i) / m_grav;
    m_Spectrum(i) = 5.0*Hs*Hs*pow(wp,4) * exp(-1.25*pow(wp/m_omega(i),4)) / (16.0*pow(m_omega(i),5));
    if(i == 0)
      m_A(i) = sqrt(2.0*m_omega(0) * m_Spectrum(i));  // Precompute components once here.
    else
      m_A(i) = sqrt(2.0*(m_omega(i)-m_omega(i-1)) * m_Spectrum(i));  // Precompute components once here.
    m_phases(i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
}


/// \brief Select PM-Spectrum (default num of phases)  
/// [DEPRECATED - This version including the unused Tp specificatoin may be removed in the future]
void LinearIncidentWave::SetToPiersonMoskowitzSpectrum(double Hs, double UnusedTp, double beta)
{
  SetToPiersonMoskowitzSpectrum(Hs, beta, DEFAULT_N_PHASES);
}

/// \brief Select PM-Spectrum (set num of phases) 
/// [DEPRECATED - This version including the unused Tp specificatoin may be removed in the future]
void LinearIncidentWave::SetToPiersonMoskowitzSpectrum(
 double Hs, double UnusedTp, double beta, int n_phases)
 {
 SetToPiersonMoskowitzSpectrum(Hs, beta, n_phases);
 }


/// \brief Select PM-Spectrum (default num of phases)  
void LinearIncidentWave::SetToPiersonMoskowitzSpectrum(double Hs, double beta)
{
  SetToPiersonMoskowitzSpectrum(Hs, beta, DEFAULT_N_PHASES);
}

/// \brief Select PM-Spectrum (set num of phases) 
void LinearIncidentWave::SetToPiersonMoskowitzSpectrum(
  double Hs, double beta, int n_phases)
{
  m_SpectrumType = WaveSpectrumType::PiersonMoskowitz;
  m_beta = beta;
  m_Hs = Hs;
  m_Tp = 2*M_PI*sqrt(m_Hs/m_grav)/0.4019;

  m_omega.resize(n_phases);
  m_k.resize(n_phases);
  m_phases.resize(n_phases);
  m_Spectrum.resize(n_phases);
  m_A.resize(n_phases);

  double w0 = sqrt(.21 * m_grav / Hs);
  double a = 0.0081;
  double b = 0.74;

  double d_omega = MAX_FREQ * 2 * M_PI / n_phases;

  for (int i = 0; i < m_k.size(); i++) {
    m_omega(i) = d_omega * (i + 1) + (0.25*d_omega*(std::rand()-RAND_MAX/2))/(RAND_MAX/2);
    m_k(i) = m_omega(i) * m_omega(i) / m_grav;
    m_Spectrum(i) = (a * m_grav * m_grav / pow(m_omega(i), 5)) * exp(-b * pow(w0 / m_omega(i), 4));
    if(i == 0)
      m_A(i) = sqrt(2.0*m_omega(0) * m_Spectrum(i));  // Precompute components once here.
    else
      m_A(i) = sqrt(2.0*(m_omega(i)-m_omega(i-1)) * m_Spectrum(i));  // Precompute components once here.
    m_phases(i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
}

/// \brief Specify Custom Spectrum (default num of phases)
void LinearIncidentWave::SetToCustomSpectrum(std::vector<double> freq, std::vector<double> S, double beta)
{
  SetToCustomSpectrum(freq, S, beta, DEFAULT_N_PHASES);
}

/// \brief Specify Custom Spectrum (set num of phases)
/// freq[Hz], S[m^2/Hz]
void LinearIncidentWave::SetToCustomSpectrum(std::vector<double> freq, std::vector<double> S, double beta, int n_phases)
{

  simple_interp::Interp1d CustomSpectrum(freq,S);  // Subsequent calculations are done in ang freq

  m_SpectrumType = WaveSpectrumType::Custom;
  m_beta = beta;

  m_omega.resize(n_phases);
  m_k.resize(n_phases);
  m_phases.resize(n_phases);
  m_Spectrum.resize(n_phases);
  m_A.resize(n_phases);

  double d_freq = MAX_FREQ / n_phases;
  Eigen::VectorXd f;
  f.resize(n_phases);

  for (int i = 0; i < m_k.size(); i++) {
    f(i) = d_freq * (i + 1) + (0.25*d_freq*(std::rand()-RAND_MAX/2))/(RAND_MAX/2);
    m_omega(i) = 2*M_PI*f(i);
    m_k(i) = m_omega(i) * m_omega(i) / m_grav;
    m_Spectrum(i) = CustomSpectrum(f(i)); //Interpolate from supplied spectrum
    if(i == 0)
      m_A(i) = sqrt(2.0 * f(0) * m_Spectrum(i));  // Precompute components once here.
    else
      m_A(i) = sqrt(2.0 * (f(i)-f(i-1)) * m_Spectrum(i));  // Precompute components once here.
    m_phases(i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
}

std::ostream & operator<<(std::ostream & out, const LinearIncidentWave & IncWave)
{
  // Since operator<< is a friend of the LinearIncidentWave class, we can access members directly.
  switch (IncWave.m_SpectrumType) {
    case WaveSpectrumType::MonoChromatic:
      std::cout << "# IncidentWave Type = Mono-Chromatic" << std::endl;
      std::cout << "# Amplitude = " << IncWave.m_Hs / 2 << std::endl;
      std::cout << "# Period = " << IncWave.m_Tp << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum.size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega.transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k.transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases.transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A.transpose() << std::endl;
      break;
    case WaveSpectrumType::PiersonMoskowitz:
      std::cout << "# IncidentWave Type = Pierson Moskowitz" << std::endl;
      std::cout << "# Hs = " << IncWave.m_Hs << std::endl;
      std::cout << "# Tp = " << IncWave.m_Tp << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum.size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega.transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k.transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases.transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum.transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A.transpose() << std::endl;
      break;
    case WaveSpectrumType::Bretschneider:
      std::cout << "# IncidentWave Type = Bretschneider" << std::endl;
      std::cout << "# Hs = " << IncWave.m_Hs << std::endl;
      std::cout << "# Tp = " << IncWave.m_Tp << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum.size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega.transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k.transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases.transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum.transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A.transpose() << std::endl;
      break;
    case WaveSpectrumType::Custom:
      std::cout << "# IncidentWave Type = Custom Spectrum";
      std::cout << "# Num Phases = " << IncWave.m_Spectrum.size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega.transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k.transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases.transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum.transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A.transpose() << std::endl;
      break;
  }
  return out;  // return std::ostream so we can chain calls to operator<<
}

double LinearIncidentWave::eta(double x, double y, double t) const
{
  double xx = x * cos(m_beta) + y * sin(m_beta);

  double eta = 0;
  for (int i = 0; i < m_A.size(); i++) {
    eta = eta + m_A(i) * cos(m_k(i) * xx - m_omega(i) * t + m_phases(i));
  }
  return eta;

// Eigen::VectorXd temp;
// temp.array()  = (xx*m_k-t*m_omega+m_phases).array().sin();
// return -m_A.dot(temp);  //With -03, these two approaches are equally fast/slow
}

double LinearIncidentWave::eta(double x, double y, double t, double *deta_dx, double *deta_dy) const
{
  return eta(x, y, t, deta_dx, deta_dy, nullptr, nullptr);
}

double LinearIncidentWave::eta(double x, double y, double t,
                               double *deta_dx, double *deta_dy,
                               double *u_east, double *v_north) const
{
  double xx = x * cos(m_beta) + y * sin(m_beta);

  double eta = 0.0;
  double deta_dxx = 0.0;

  // Eulerian along-wave contribution
  double u_along = 0.0;
  // Stokes along-wave contribution
  double u_stokes_along = 0.0;

  for (int i = 0; i < m_A.size(); i++) {
    double k = m_k(i);
    double omega = m_omega(i);
    double a = m_A(i);
    double phase = m_phases(i);

    double arg = k * xx - omega * t + phase;
    double cosarg = cos(arg);
    double sinarg = sin(arg);

    // freesurface heave
    eta += a * cosarg;

    // water plane slope
    deta_dxx -= k * a * sinarg;

    // Eulerian along-wave surface velocity (assume deep water)
    // can just use omega directly since deep water dispersion w^2 = gk
    u_along += a * omega * cosarg;
  }

  // water plane slope
  if (deta_dx) *deta_dx = deta_dxx*cos(m_beta);  // deta/dx
  if (deta_dy) *deta_dy = deta_dxx*sin(m_beta);  // deta/dy

  // u/v Eulerian surface velocities
  if (u_east) *u_east = u_along * cos(m_beta);
  if (v_north) *v_north = u_along * sin(m_beta);

  return eta;
}

double LinearIncidentWave::etadot(double x, double y, double t) const
{
  double xx = x * cos(m_beta) + y * sin(m_beta);

  double etadot = 0;
  for (int i = 0; i < m_A.size(); i++) {
    etadot = etadot + m_omega(i) * m_A(i) * sin(m_k(i) * xx - m_omega(i) * t + m_phases(i));
  }

  return etadot;
}

/// \brief Returns Version String
std::string LinearIncidentWave::Version() { return PROJECT_VER; }

/// \brief Returns Major Version Number
int LinearIncidentWave::MajorVersionNumber() { return PROJECT_VER_MAJOR; }

/// \brief Returns Minor Version Number
int LinearIncidentWave::MinorVersionNumber() { return PROJECT_VER_MINOR; }

/// \brief Returns Major Version Number
int LinearIncidentWave::PatchVersionNumber() { return PROJECT_VER_PATCH; }
