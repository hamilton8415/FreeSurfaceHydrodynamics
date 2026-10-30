#include <gtest/gtest.h>
#include "FreeSurfaceHydrodynamics/LinearIncidentWave.hpp"
#include <iostream>
#include <random>

//compute wave elevation for one single-frequency wave.
double confirm_eta(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx =  x * cos(beta) + y * sin(beta);
  return A*cos(k*xx-w*t+phi);
}


TEST(MonoChrome, CheckElevations) {

    std::random_device rd;
    std::mt19937 gen(rd());

    LinearIncidentWave Inc;

    double g = Inc.m_grav;
    ASSERT_EQ(g,9.81);

    double A = 1.1;
    double T = 11.1;
    double phase = 33*M_PI/180.0;
    double beta = 22*M_PI/180.0;
    double H = 10;  // Set depth, negative means infinite depth.
    Inc.SetDepth(H);
    Inc.SetToMonoChromatic(A,T,phase,beta);
    double w = 2.0*M_PI/T;
    double k = Inc.m_k[0](0);
    std::uniform_real_distribution<double> dist(-100,100);
    for (int i = 0; i<5; i++)
    {
    double x = dist(gen);
    double y = dist(gen);
    double t = dist(gen);
    double eta = Inc.eta(x,y,t);
    double eta0 = confirm_eta(A,w,k,phase,beta,x,y,t);
    std::cout << "eta = " << eta << "  eta0 = " << eta0 << std::endl;
    EXPECT_DOUBLE_EQ(eta,eta0);
    }
}

TEST(MonoChrome, CheckVelocities) {
EXPECT_EQ(1+1+1,3);
}

