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
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <FreeSurfaceHydrodynamics/config.h>
#include <FreeSurfaceHydrodynamics/interp1d.hpp>
#include <FreeSurfaceHydrodynamics/LinearIncidentWave.hpp>



// Define a functor for the system of nonlinear equations
struct DispersionFunctor {
  //Constructor
  DispersionFunctor(double omega, double depth, double grav) : w(omega), h(depth),g(grav){};
  // Number of variables and equations (n x n)
  int inputs() const { return 1; }
  int values() const { return 1; }

    // Compute the residual vector fvec given x
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
        // Example system:
        // f1(x1, x2) = x1^2 + x2^2 - 2 = 0
        // f2(x1, x2) = x1 - x2 = 0
        fvec(0) = x(0)*tanh(x(0)*h)-w*w/g;
        return 0;
    }
    private:
    double h;
    double g;
    double w;
};


/// \brief Constructor, defaults to monotchromatic wave and default gravity and density
LinearIncidentWave::LinearIncidentWave() : m_grav(9.81), m_rho(1025), m_depth(-1.0)
{
  std::srand(time(0));
}

/// \brief Constructor, sets seed and defaults to monotchromatic wave and default gravity and density
LinearIncidentWave::LinearIncidentWave(unsigned int seed) : m_grav(9.81), m_rho(1025), m_depth(-1.0)
{
  SetSeed(seed);
}

/// \brief Sets seed to specified value
void LinearIncidentWave::SetSeed(unsigned int seed)
{
  if (seed == 0)
    std::srand(time(0));
  else
    std::srand(seed);
}

/// \brief Sets depth, negative values indicate infinite depth
void LinearIncidentWave::SetDepth(double h)
{
  m_depth = h;
}

/// \brief Select single frequency wave
void LinearIncidentWave::SetToMonoChromatic(double A, double T, double phase, double beta)
{
  m_SpectrumType.push_back(WaveSpectrumType::MonoChromatic);
  m_Hs.push_back(2 * A);
  m_Tp.push_back(T);
  m_beta.push_back(beta);
  Eigen::VectorXd tmp_array(1); // Create array and push a copy and put a copy in each relevant std::vector
  m_omega.push_back(tmp_array);
  m_k.push_back(tmp_array);
  m_phases.push_back(tmp_array);
  m_Spectrum.push_back(tmp_array);
  m_A.push_back(tmp_array);
  m_phases[NumWaveComponents](0) = phase;
  m_omega[NumWaveComponents](0) = 2 * M_PI / T;
  m_k[NumWaveComponents](0) = WaveNumber(m_omega[NumWaveComponents](0));
  m_A[NumWaveComponents](0) = A;
  NumWaveComponents++; // Adding a wave component
}

/// \brief Select PM-Spectrum (default num of phases)
void LinearIncidentWave::SetToBretschneiderSpectrum(double Hs, double Tp, double beta)
{
  SetToBretschneiderSpectrum(Hs, Tp, beta, DEFAULT_N_PHASES);
}

void LinearIncidentWave::SetToBretschneiderSpectrumWithCos2Spreading(double Hs, double Tp, double beta_0, int spreading_factor, int n_sectors)
{
  SetToBretschneiderSpectrumWithCos2Spreading(Hs, Tp, beta_0, DEFAULT_N_PHASES, spreading_factor, n_sectors);
}

void LinearIncidentWave::SetToBretschneiderSpectrumWithCos2Spreading(double Hs, double Tp, double beta_0, int n_phases, int spreading_factor, int n_sectors)
{
  double d_beta = 2.0 * M_PI / n_sectors;
  for (int n = 1; n < n_sectors; n++) // start with n = 1 b/c the reciprocal heading wave is of identically zero  (cos(pi/2) = 0)
  {
    double beta = beta_0 - (n * d_beta - M_PI);
    double D = std::pow(cos((beta - beta_0) / 2), 2 * spreading_factor) * std::tgamma(1.0 + spreading_factor) / (2 * sqrt(M_PI) * std::tgamma(0.5 + spreading_factor));
    SetToBretschneiderSpectrum(Hs * sqrt(D), Tp, beta, n_phases); // Spectrum energy scales as the square of Hs, so sqrt(D) introduces a factor of D into the Spectrum
  }
}

/// \brief Select Bretschneider Spectrum (set num of phases)
void LinearIncidentWave::SetToBretschneiderSpectrum(
    double Hs, double Tp, double beta,
    int n_phases)
{
  m_SpectrumType.push_back(WaveSpectrumType::Bretschneider);
  m_Hs.push_back(Hs);
  m_Tp.push_back(Tp);
  m_beta.push_back(beta);
  Eigen::VectorXd tmp_array(n_phases); // Create array and push a copy and put a copy in each relevant std::vector
  m_omega.push_back(tmp_array);
  m_k.push_back(tmp_array);
  m_phases.push_back(tmp_array);
  m_Spectrum.push_back(tmp_array);
  m_A.push_back(tmp_array);

  double wp = 2.0 * M_PI / (1.2957 * Tp);
  double d_omega = MAX_FREQ * 2 * M_PI / n_phases;

  for (int i = 0; i < m_k[NumWaveComponents].size(); i++)
  {
    m_omega[NumWaveComponents](i) = d_omega * (i + 1) + (0.25 * d_omega * (std::rand() - RAND_MAX / 2)) / (RAND_MAX / 2);
    m_k[NumWaveComponents](i) = WaveNumber(m_omega[NumWaveComponents](i));
    m_Spectrum[NumWaveComponents](i) = 5.0 * Hs * Hs * pow(wp, 4) * exp(-1.25 * pow(wp / m_omega[NumWaveComponents](i), 4)) / (16.0 * pow(m_omega[NumWaveComponents](i), 5));
    if (i == 0)
      m_A[NumWaveComponents](i) = sqrt(2.0 * m_omega[NumWaveComponents](0) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    else
      m_A[NumWaveComponents](i) = sqrt(2.0 * (m_omega[NumWaveComponents](i) - m_omega[NumWaveComponents](i - 1)) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    m_phases[NumWaveComponents](i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
  NumWaveComponents++; // Adding a wave component
}

void LinearIncidentWave::SetToPiersonMoskowitzSpectrumWithCos2Spreading(double Hs, double beta_0, int spreading_factor, int n_sectors)
{
  SetToPiersonMoskowitzSpectrumWithCos2Spreading(Hs, beta_0, DEFAULT_N_PHASES, spreading_factor, n_sectors);
}

void LinearIncidentWave::SetToPiersonMoskowitzSpectrumWithCos2Spreading(double Hs, double beta_0, int n_phases, int spreading_factor, int n_sectors)
{
  double d_beta = 2.0 * M_PI / n_sectors;
  for (int n = 1; n < n_sectors; n++) // start with n = 1 b/c the reciprocal heading wave is of identically zero  (cos(pi/2) = 0)
  {
    double beta = beta_0 - (n * d_beta - M_PI);
    double D = std::pow(cos((beta - beta_0) / 2), 2 * spreading_factor) * std::tgamma(1.0 + spreading_factor) / (2 * sqrt(M_PI) * std::tgamma(0.5 + spreading_factor));
    SetToPiersonMoskowitzSpectrum(Hs * sqrt(D), beta, n_phases); // Spectrum energy scales as the square of Hs, so sqrt(D) introduces a factor of D into the Spectrum
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
  m_SpectrumType.push_back(WaveSpectrumType::PiersonMoskowitz);
  m_Hs.push_back(Hs);
  m_Tp.push_back(2 * M_PI * sqrt(Hs / m_grav) / 0.4019);
  m_beta.push_back(beta);
  Eigen::VectorXd tmp_array(n_phases); // Create array and push a copy and put a copy in each relevant std::vector
  m_omega.push_back(tmp_array);
  m_k.push_back(tmp_array);
  m_phases.push_back(tmp_array);
  m_Spectrum.push_back(tmp_array);
  m_A.push_back(tmp_array);

  double w0 = sqrt(.21 * m_grav / Hs);
  double a = 0.0081;
  double b = 0.74;

  double d_omega = MAX_FREQ * 2 * M_PI / n_phases;

  for (int i = 0; i < m_k[NumWaveComponents].size(); i++)
  {
    m_omega[NumWaveComponents](i) = d_omega * (i + 1) + (0.25 * d_omega * (std::rand() - RAND_MAX / 2)) / (RAND_MAX / 2);
    m_k[NumWaveComponents](i) = WaveNumber(m_omega[NumWaveComponents](i));
    m_Spectrum[NumWaveComponents](i) = (a * m_grav * m_grav / pow(m_omega[NumWaveComponents](i), 5)) * exp(-b * pow(w0 / m_omega[NumWaveComponents](i), 4));
    if (i == 0)
      m_A[NumWaveComponents](i) = sqrt(2.0 * m_omega[NumWaveComponents](0) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    else
      m_A[NumWaveComponents](i) = sqrt(2.0 * (m_omega[NumWaveComponents](i) - m_omega[NumWaveComponents](i - 1)) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    m_phases[NumWaveComponents](i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
  NumWaveComponents++; // Adding a wave component
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

  simple_interp::Interp1d CustomSpectrum(freq, S); // Subsequent calculations are done in ang freq

  m_SpectrumType.push_back(WaveSpectrumType::Custom);
  m_Hs.push_back(0.0); // Not defined, could be computed from supplied spectrum
  m_Tp.push_back(0.0); // Not defined, could be computed from supplied spectrum
  m_beta.push_back(beta);
  Eigen::VectorXd tmp_array(n_phases); // Create array and push a copy and put a copy in each relevant std::vector
  m_omega.push_back(tmp_array);
  m_k.push_back(tmp_array);
  m_phases.push_back(tmp_array);
  m_Spectrum.push_back(tmp_array);
  m_A.push_back(tmp_array);

  double d_freq = MAX_FREQ / n_phases;
  Eigen::VectorXd f;
  f.resize(n_phases);

  for (int i = 0; i < m_k[NumWaveComponents].size(); i++)
  {
    f(i) = d_freq * (i + 1) + (0.25 * d_freq * (std::rand() - RAND_MAX / 2)) / (RAND_MAX / 2);
    m_omega[NumWaveComponents](i) = 2 * M_PI * f(i);
    m_k[NumWaveComponents](i) = WaveNumber(m_omega[NumWaveComponents](i));
    m_Spectrum[NumWaveComponents](i) = CustomSpectrum(f(i)); // Interpolate from supplied spectrum
    if (i == 0)
      m_A[NumWaveComponents](i) = sqrt(2.0 * f(0) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    else
      m_A[NumWaveComponents](i) = sqrt(2.0 * (f(i) - f(i - 1)) * m_Spectrum[NumWaveComponents](i)); // Precompute components once here.
    m_phases[NumWaveComponents](i) = (2 * M_PI * std::rand()) / RAND_MAX;
  }
  NumWaveComponents++; // Adding a wave component
}

/*
void PrintWaveInfo()
{
  std::cout << "# IncidentWave consists of " << IncWave.NumWaveComponents << " Components" << std::endl;
  for (int n = 0; n < IncWave.NumWaveComponents; n++)
  {
    switch (IncWave.m_SpectrumType[n])
    {
    case WaveSpectrumType::MonoChromatic:
      std::cout << "# IncidentWave Type = Mono-Chromatic" << std::endl;
      std::cout << "# Amplitude = " << IncWave.m_Hs[n] / 2 << std::endl;
      std::cout << "# Period = " << IncWave.m_Tp[n] << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::PiersonMoskowitz:
      std::cout << "# IncidentWave Type = Pierson Moskowitz" << std::endl;
      std::cout << "# Hs = " << IncWave.m_Hs[n] << std::endl;
      std::cout << "# Tp = " << IncWave.m_Tp[n] << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::Bretschneider:
      std::cout << "# IncidentWave Type = Bretschneider" << std::endl;
      std::cout << "# Hs = " << IncWave.m_Hs[n] << std::endl;
      std::cout << "# Tp = " << IncWave.m_Tp[n] << std::endl;
      std::cout << "# Beta = " << IncWave.m_beta[n] << std::endl;
      std::cout << "# Num Phases = " << IncWave.m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::Custom:
      std::cout << "# IncidentWave Type = Custom Spectrum";
      std::cout << "# Num Phases = " << IncWave.m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << IncWave.m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << IncWave.m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << IncWave.m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << IncWave.m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << IncWave.m_A[n].transpose() << std::endl;
      break;
    }
    std::cout << std::endl;
  }
    */

void LinearIncidentWave::PrintWaveInfo() const
{
  std::cout << "# IncidentWave consists of " << NumWaveComponents << " Components" << std::endl;
  for (int n = 0; n < NumWaveComponents; n++)
  {
    switch (m_SpectrumType[n])
    {
    case WaveSpectrumType::MonoChromatic:
      std::cout << "# IncidentWave Type = Mono-Chromatic" << std::endl;
      std::cout << "# Amplitude = " << m_Hs[n] / 2 << std::endl;
      std::cout << "# Period = " << m_Tp[n] << std::endl;
      std::cout << "# Num Phases = " << m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << m_phases[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::PiersonMoskowitz:
      std::cout << "# IncidentWave Type = Pierson Moskowitz" << std::endl;
      std::cout << "# Hs = " << m_Hs[n] << std::endl;
      std::cout << "# Tp = " << m_Tp[n] << std::endl;
      std::cout << "# Num Phases = " << m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::Bretschneider:
      std::cout << "# IncidentWave Type = Bretschneider" << std::endl;
      std::cout << "# Hs = " << m_Hs[n] << std::endl;
      std::cout << "# Tp = " << m_Tp[n] << std::endl;
      std::cout << "# Beta = " << m_beta[n] << std::endl;
      std::cout << "# Num Phases = " << m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << m_A[n].transpose() << std::endl;
      break;
    case WaveSpectrumType::Custom:
      std::cout << "# IncidentWave Type = Custom Spectrum";
      std::cout << "# Num Phases = " << m_Spectrum[n].size() << std::endl;
      std::cout << "# Wave Freq = " << m_omega[n].transpose() << std::endl;
      std::cout << "# Wave Numbers = " << m_k[n].transpose() << std::endl;
      std::cout << "# Phases = " << m_phases[n].transpose() << std::endl;
      std::cout << "# Spectrum = " << m_Spectrum[n].transpose() << std::endl;
      std::cout << "# Component Amplitudes = " << m_A[n].transpose() << std::endl;
      break;
    }
    std::cout << std::endl;
  }
}

std::ostream &operator<<(std::ostream &out, const LinearIncidentWave &IncWave)
{
  IncWave.PrintWaveInfo();
  return out; // return std::ostream so we can chain calls to operator<<
}

double LinearIncidentWave::eta(double x, double y, double t,
                               double *deta_dx, double *deta_dy,
                               double *u_east, double *v_north, int n) const
{
  double eta = 0.0;
  if (deta_dx)
    *deta_dx = 0.0;
  if (deta_dy)
    *deta_dy = 0.0;
  if (u_east)
    *u_east = 0.0;
  if (v_north)
    *v_north = 0.0;

  if ((n < 0) || (n >= NumWaveComponents)) // Test to be sure wave component argument is valid, return zero otherwise
    return eta;

  double xx = x * cos(M_PI+m_beta[n]) + y * sin(M_PI+m_beta[n]);

  // Per-direction-component slope / along-wave velocity contributions.
  // These must be reset for each directional component before projecting
  // into global East/North. Accumulating them across `n` and then applying
  // the current component heading mixes previous sectors into the wrong
  // direction, which especially corrupts multi-sector directional seas.
  double deta_dxx = 0.0;
  double u_along = 0.0;

  for (int i = 0; i < m_A[n].size(); i++)
  {
    double k = m_k[n](i);
    double omega = m_omega[n](i);
    double a = m_A[n](i);
    double phase = m_phases[n](i);

    double arg = k * xx - omega * t + phase;
    double cosarg = cos(arg);

    // freesurface heave
    eta += a * cosarg;

    // water plane slope
    if (deta_dx || deta_dy) // Only compute this if needed
    {
      double sinarg = sin(arg);
      deta_dxx -= k * a * sinarg;
    }

    // Eulerian along-wave surface velocity
    u_along += (m_grav * k * a / omega) * cosarg;
  }

  // water plane slope
  if (deta_dx)
    *deta_dx += deta_dxx * cos(m_beta[n]); // deta/dx
  if (deta_dy)
    *deta_dy += deta_dxx * sin(m_beta[n]); // deta/dy

  // u/v Eulerian surface velocities
  if (u_east)
    *u_east += u_along * cos(m_beta[n]);
  if (v_north)
    *v_north += u_along * sin(m_beta[n]);

  return eta;
}

double LinearIncidentWave::eta(double x, double y, double t, double *deta_dx, double *deta_dy, int n) const
{
  return eta(x, y, t, deta_dx, deta_dy, nullptr, nullptr, n);
}

double LinearIncidentWave::eta(double x, double y, double t, int n) const
{
  return eta(x, y, t, nullptr, nullptr, nullptr, nullptr, n);
}

double LinearIncidentWave::eta(double x, double y, double t,
                               double *deta_dx, double *deta_dy,
                               double *u_east, double *v_north) const
{
  double eta_sum = 0;
  double deta_dx_sum = 0;
  double deta_dy_sum = 0;
  double u_east_sum = 0;
  double v_north_sum = 0;

  for (int n = 0; n < NumWaveComponents; n++)
  {
    eta_sum += eta(x, y, t, deta_dx, deta_dy, u_east, v_north, n);
    if (deta_dx)
      deta_dx_sum += *deta_dx;
    if (deta_dy)
      deta_dy_sum += *deta_dy;
    if (u_east)
      u_east_sum += *u_east;
    if (v_north)
      v_north_sum += *v_north;

    if (deta_dx)
      *deta_dx = deta_dx_sum;
    if (deta_dy)
      *deta_dy = deta_dy_sum;
    if (u_east)
      *u_east = u_east_sum;
    if (v_north)
      *v_north = v_north_sum;
  }
  return eta_sum;
}

double LinearIncidentWave::eta(double x, double y, double t, double *deta_dx, double *deta_dy) const
{
  return eta(x, y, t, deta_dx, deta_dy, nullptr, nullptr);
}

double LinearIncidentWave::eta(double x, double y, double t) const
{
  return eta(x, y, t, nullptr, nullptr, nullptr, nullptr);
}

double LinearIncidentWave::etadot(double x, double y, double t, int n) const
{
  double etadot = 0;

  if ((n < 0) || (n >= NumWaveComponents)) // Test to be sure wave component argument is valid, return zero otherwise
    return etadot;

  double xx = x * cos(m_beta[n]) + y * sin(m_beta[n]);

  for (int i = 0; i < m_A[n].size(); i++)
  {
    etadot += m_omega[n](i) * m_A[n](i) * sin(m_k[n](i) * xx - m_omega[n](i) * t + m_phases[n](i));
  }
  return etadot;
}

double LinearIncidentWave::etadot(double x, double y, double t) const
{
  double etadot_sum = 0;
  for (int n = 0; n < NumWaveComponents; n++)
    etadot_sum += etadot(x, y, t, n);

  return etadot_sum;
}

void LinearIncidentWave::vel(double x, double y, double z, double t, double *vel_x, double *vel_y, double *vel_z) const
{
 *vel_x = 0;
 *vel_y = 0;
 *vel_z = 0;

  for (int n = 0; n < NumWaveComponents; n++)
  {
    double loc_vel_x;
    double loc_vel_y;
    double loc_vel_z;
    vel(x, y, z, t, &loc_vel_x, &loc_vel_y, &loc_vel_z , n);
    *vel_x += loc_vel_x;
    *vel_y += loc_vel_y;
    *vel_z += loc_vel_z;
  }
  return;
}

void LinearIncidentWave::vel(double x, double y, double z, double t, double *vel_x, double *vel_y, double *vel_z, int n) const
{

  double v_z = 0.0;
  double v_along = 0.0;

  if ((n < 0) || (n >= NumWaveComponents)) // Test to be sure wave component argument is valid, return zero otherwise
  {
    *vel_x = v_along;
    *vel_y = v_along;
    *vel_z = v_z;
    return;
  } 

  double xx = x * cos(m_beta[n]) + y * sin(m_beta[n]);

  for (int i = 0; i < m_A[n].size(); i++)
  {
    double arg = m_k[n](i) * xx - m_omega[n](i) * t + m_phases[n](i);
    double sinarg = sin(arg);
    double cosarg = cos(arg);
    if (m_depth < 0) // Infinite depth
    {
      v_z += m_omega[n](i) * m_A[n](i) * exp(m_k[n](i) * z) * sinarg;
      v_along += m_omega[n](i) * m_A[n](i) * exp(m_k[n](i) * z) * cosarg;
    }
    else
    {
      v_z += m_grav*m_A[n](i)*m_k[n](i) * (sinh(m_k[n](i)*(z+m_depth)) / cosh(m_k[n](i)*m_depth)) * sinarg / m_omega[n](i);
      v_along += m_grav*m_A[n](i)*m_k[n](i) * (cosh(m_k[n](i)*(z+m_depth)) / cosh(m_k[n](i)*m_depth)) * cosarg / m_omega[n](i);
    }
  }
  // u/v Eulerian surface velocities
    *vel_x = v_along * cos(m_beta[n]);
    *vel_y = v_along * sin(m_beta[n]);
    *vel_z = v_z;
  
  return;
}

  double LinearIncidentWave::WaveNumber(double w)
  {
    DispersionFunctor functor(w, m_depth, m_grav);
    // HybridNonLinearSolver uses Powell's hybrid method
    Eigen::HybridNonLinearSolver<DispersionFunctor> solver(functor);

    Eigen::VectorXd k(1);
    k(0) = w*w/m_grav;  // Infinite depth result, used as initial guess if needed
    if(m_depth >=0)
    {
    int status = solver.hybrd1(k);   // Solve using hybrd1 (numerical differentiation)
    //std::cout << "status = " << status << "  k = " << k(0) <<  "  sqrt(g*k*tanh(kh) = " << sqrt(m_grav*k(0)*tanh(k(0)*m_depth)) << std::endl;
    }
   return k(0);
  }

/// \brief Returns Version String
std::string LinearIncidentWave::Version() { return PROJECT_VER; }

/// \brief Returns Major Version Number
int LinearIncidentWave::MajorVersionNumber() { return PROJECT_VER_MAJOR; }

/// \brief Returns Minor Version Number
int LinearIncidentWave::MinorVersionNumber() { return PROJECT_VER_MINOR; }

/// \brief Returns Major Version Number
int LinearIncidentWave::PatchVersionNumber() { return PROJECT_VER_PATCH; }
