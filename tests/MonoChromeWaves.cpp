#include <gtest/gtest.h>
#include "FreeSurfaceHydrodynamics/LinearIncidentWave.hpp"
#include <iostream>
#include <random>

// compute wave elevation for one single-frequency wave.
double confirm_eta_ENU(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx = -y * cos(beta) - x * sin(beta);
  return A * cos(k * xx - w * t + phi);
}

// compute wave elevation for one single-frequency wave.
double confirm_eta_NED(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx = -x * cos(beta) - y * sin(beta);
  return -A * cos(k * xx - w * t + phi);
}

// compute wave elevation slope for one single-frequency wave.
double confirm_deta_dx_ENU(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx = -y * cos(beta) - x * sin(beta);
  return A * k * sin(k * xx - w * t + phi) * sin(beta);
}

// compute wave elevation slope for one single-frequency wave.
double confirm_deta_dy_ENU(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx = -y * cos(beta) - x * sin(beta);
  return A * k * sin(k * xx - w * t + phi) * cos(beta);
}

// compute surface velocity for one single-frequency wave.
double confirm_ueast_ENU(double A, double w, double k, double phi, double beta, double g, double x, double y, double t)
{
  double xx = -y * cos(beta) - x * sin(beta);
  return -(g * k * A / w) * cos(k * xx - w * t + phi) * sin(beta);
}

// compute surface velocity for one single-frequency wave.
double confirm_unorth_ENU(double A, double w, double k, double phi, double beta, double g, double x, double y, double t)
{
  double xx = -y * cos(beta) - x * sin(beta);
  return -(g * k * A / w) * cos(k * xx - w * t + phi) * cos(beta);
}

TEST(MonoChrome, CheckElevations)
{
  std::random_device rd;
  std::mt19937 gen(rd());

  LinearIncidentWave Inc;

  double A = 1.1;
  double T = 11.1;
  double phase = 33 * M_PI / 180.0;
  double beta = 22 * M_PI / 180.0;
  double H = 10; // Set depth, negative means infinite depth.
  // Inc.SetDepth(H);
  Inc.SetToMonoChromatic(A, T, phase, beta);
  double w = 2.0 * M_PI / T;
  double k = Inc.m_k[0](0);
  std::uniform_real_distribution<double> dist(-100, 100);
  for (int i = 0; i < 100; i++)
  {
    double x = dist(gen);
    double y = dist(gen);
    double t = dist(gen);
    std::cout << t << std::endl;
    Inc.SetCoordSysType(CoordSysType::ENU);
    double eta = Inc.eta(x, y, t);
    double eta0 = confirm_eta_ENU(A, w, k, phase, beta, x, y, t);
    std::cout << eta << "  " << eta0 << std::endl;
    EXPECT_DOUBLE_EQ(eta, eta0);

    double deta_dx, deta_dy;
    double u_east, u_north;
    eta = Inc.eta(x, y, t, &deta_dx, &deta_dy, &u_east, &u_north);
    double deta0_dx = confirm_deta_dx_ENU(A, w, k, phase, beta, x, y, t);
    double deta0_dy = confirm_deta_dy_ENU(A, w, k, phase, beta, x, y, t);
    double u0_east = confirm_ueast_ENU(A, w, k, phase, beta, Inc.m_grav, x, y, t);
    double u0_north = confirm_unorth_ENU(A, w, k, phase, beta, Inc.m_grav, x, y, t);
    EXPECT_DOUBLE_EQ(eta, eta0);
    EXPECT_DOUBLE_EQ(deta_dx, deta0_dx);
    EXPECT_DOUBLE_EQ(deta_dy, deta0_dy);
    EXPECT_DOUBLE_EQ(u_east, u0_east);
    EXPECT_DOUBLE_EQ(u_north, u0_north);

    eta = Inc.eta(x, y, t, &deta_dx, &deta_dy);
    deta0_dx = confirm_deta_dx_ENU(A, w, k, phase, beta, x, y, t);
    deta0_dy = confirm_deta_dy_ENU(A, w, k, phase, beta, x, y, t);
    EXPECT_DOUBLE_EQ(eta, eta0);
    EXPECT_DOUBLE_EQ(deta_dx, deta0_dx);
    EXPECT_DOUBLE_EQ(deta_dy, deta0_dy);

    Inc.SetCoordSysType(CoordSysType::NED);
    eta = Inc.eta(y, x, t);
    eta0 = confirm_eta_NED(A, w, k, phase, beta, y, x, t);
    EXPECT_DOUBLE_EQ(eta, eta0);


  }
}

TEST(MonoChrome, CheckVelocities)
{
  EXPECT_EQ(1 + 1 + 1, 3);
}
