#include <gtest/gtest.h>
#include "FreeSurfaceHydrodynamics/LinearIncidentWave.hpp"
#include <iostream>
#include <random>

TEST(CheckConstantsAndRelations, CheckConstants)
{
  LinearIncidentWave Inc;

  // Check Gravity Constant
  double g = Inc.m_grav;
  ASSERT_EQ(9.81, g);

  // Check density Constant
  double rho = Inc.m_rho;
  ASSERT_EQ(1025.0, rho);
}

TEST(CheckConstantsAndRelations, CheckDispersionRelations)
{

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> ddist(.1,4);
  double beta = 0;
  double phase = 0;
  double spreading_factor = 45 * M_PI / 180;
  std::uniform_int_distribution<int> dist(10,100);

  for (int k = -40; k <= 40; k += 40) // Run for three depths
  {
    LinearIncidentWave Inc;
    double H = (double) k;
    if (k != 0)
      Inc.SetDepth(H);
    double A = ddist(gen);
    double T = 4.0*ddist(gen);
    int n_phases;
    int n_sectors;
    int ExpectedWaveComponents = 0;

    Inc.SetToMonoChromatic(A, T, phase, beta); // Amp, Period, phase, wave_direction
    ExpectedWaveComponents += 1;

    Inc.SetToPiersonMoskowitzSpectrum(2 * A, beta); // Hs, wave_direction
    ExpectedWaveComponents += 1;
   
    n_phases = dist(gen);
    Inc.SetToPiersonMoskowitzSpectrum(2 * A, beta, n_phases); // Hs, wave_direction, number of phases
    ExpectedWaveComponents += 1;

    Inc.SetToPiersonMoskowitzSpectrum(2 * A, T, beta); // Hs, Unused Tp, wave_direction
    ExpectedWaveComponents += 1;

    n_phases = dist(gen);
    Inc.SetToPiersonMoskowitzSpectrum(2 * A, T, beta, n_phases); // Hs, Unused Tp, wave_direction, number of phases
    ExpectedWaveComponents += 1;

    n_phases = dist(gen);
    n_sectors = dist(gen);
    Inc.SetToPiersonMoskowitzSpectrumWithCos2Spreading(2 * A, beta, n_phases, spreading_factor, n_sectors); // Hs, wave_direction, number of phases, spreading angle, number of sectors
    ExpectedWaveComponents += n_sectors - 1;

    n_sectors = dist(gen);
    Inc.SetToPiersonMoskowitzSpectrumWithCos2Spreading(2 * A, beta, spreading_factor, n_sectors); // Hs, wave_direction, spreading angle, number of sectors
    ExpectedWaveComponents += n_sectors - 1;

    Inc.SetToBretschneiderSpectrum(2 * A, T, beta); // Hs, Tp, wave_direction
    ExpectedWaveComponents += 1;

    n_phases = dist(gen);
    Inc.SetToBretschneiderSpectrum(2 * A, T, beta, n_phases); // Hs, Tp, wave_direction, number of phases
    ExpectedWaveComponents += 1;

    n_sectors = dist(gen);
    Inc.SetToBretschneiderSpectrumWithCos2Spreading(2 * A, T, beta, spreading_factor, n_sectors); // Hs, Tp, wave_direction, spreading factor, number of sectors
    ExpectedWaveComponents += n_sectors - 1;

    n_phases = dist(gen);
    n_sectors = dist(gen);
    Inc.SetToBretschneiderSpectrumWithCos2Spreading(2 * A, T, beta, n_phases, spreading_factor, n_sectors); // Hs, Tp, wave_direction, n_phases, spreading factor, number of sectors
    ExpectedWaveComponents += n_sectors - 1;

    EXPECT_EQ(Inc.NumWaveComponents, ExpectedWaveComponents); // Total number of wave-trains added, depends on waves added...

    EXPECT_EQ(Inc.NumWaveComponents, Inc.m_beta.size());
    EXPECT_EQ(Inc.NumWaveComponents, Inc.m_Hs.size());
    EXPECT_EQ(Inc.NumWaveComponents, Inc.m_Tp.size());
    EXPECT_EQ(Inc.NumWaveComponents, Inc.m_SpectrumType.size());

    for (int n = 0; n < Inc.NumWaveComponents; n++)
    {
      int n_phases = Inc.m_phases[n].size();
      EXPECT_EQ(n_phases, Inc.m_Spectrum[n].size());
      EXPECT_EQ(n_phases, Inc.m_A[n].size());
      EXPECT_EQ(n_phases, Inc.m_omega[n].size());
      EXPECT_EQ(n_phases, Inc.m_k[n].size());
      double g = Inc.m_grav;
      for (int i = 0; i < n_phases; i++)
      {
        double w = Inc.m_omega[n](i);
        double k = Inc.m_k[n](i);
        double dispersion_error;
        if (H <= 0) // Infinite Depth
          dispersion_error = k - w * w / g;
        else
          dispersion_error = k * tanh(k * H) - w * w / g;
        EXPECT_NEAR(0.0, dispersion_error, 1e-10); // Computed wave-number is a numerical solution of an implicit relation, so may not be perfect.
      }
    }
  }
}
