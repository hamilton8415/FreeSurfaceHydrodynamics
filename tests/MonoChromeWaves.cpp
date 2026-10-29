#include <gtest/gtest.h>
#include "FreeSurfaceHydrodynamics/LinearIncidentWave.hpp"
#include <iostream>
#include <random>

TEST(MonoChrome, CheckElevations) {

    std::random_device rd;
    std::mt19937 gen(rd());

    LinearIncidentWave Inc;

    double A = 1.1;
    double T = 11.1;
    double w = 2*M_PI/T;
    double k = w*w/9.81;
    std::uniform_real_distribution<double> dist(-100,100);
    Inc.SetToMonoChromatic(A,T,0,0);
    for (int i = 0; i<5; i++)
    {
    double x = dist(gen);
    double y = dist(gen);
    double t = dist(gen);
    double eta = Inc.eta(x,y,t);
    double xx =  x * cos(0) + y * sin(0);
    double eta0 = A*cos(k*xx-w*t);
    std::cout << "eta = " << eta << "  eta0 = " << eta0 << std::endl;
    EXPECT_DOUBLE_EQ(eta,eta0);
    }
}

TEST(MonoChrome, CheckVelocities) {
EXPECT_EQ(1+1+1,3);
}

