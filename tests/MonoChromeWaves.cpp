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

//compute wave elevation slope for one single-frequency wave.
double confirm_deta_dx(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx =  x * cos(beta) + y * sin(beta);
  return -A*k*sin(k*xx-w*t+phi)*cos(beta);
}

//compute wave elevation slope for one single-frequency wave.
double confirm_deta_dy(double A, double w, double k, double phi, double beta, double x, double y, double t)
{
  double xx =  x * cos(beta) + y * sin(beta);
  return -A*k*sin(k*xx-w*t+phi)*sin(beta);
}

//compute surface velocity for one single-frequency wave.
double confirm_ueast(double A, double w, double k, double phi, double beta, double g, double x, double y, double t)
{
  double xx =  x * cos(beta) + y * sin(beta);
  return (g*k*A/w)*cos(k*xx-w*t+phi)*cos(beta);
}

//compute surface velocity for one single-frequency wave.
double confirm_unorth(double A, double w, double k, double phi, double beta, double g, double x, double y, double t)
{
  double xx =  x * cos(beta) + y * sin(beta);
  return (g*k*A/w)*cos(k*xx-w*t+phi)*sin(beta);
}

TEST(MonoChrome, CheckElevations) {

    std::random_device rd;
    std::mt19937 gen(rd());

    LinearIncidentWave Inc;

    double A = 1.1;
    double T = 11.1;
    double phase = 33*M_PI/180.0;
    double beta = 22*M_PI/180.0;
    double H = 10;  // Set depth, negative means infinite depth.
    //Inc.SetDepth(H);
    Inc.SetToMonoChromatic(A,T,phase,beta);
    double w = 2.0*M_PI/T;
    double k = Inc.m_k[0](0);
    std::uniform_real_distribution<double> dist(-100,100);
    for (int i = 0; i<1000; i++)
    {
    double x = dist(gen);
    double y = dist(gen);
    double t = dist(gen);
    std::cout << t << std::endl;
    double eta = Inc.eta(x,y,t);
    double eta0 = confirm_eta(A,w,k,phase,beta,x,y,t);
    EXPECT_DOUBLE_EQ(eta,eta0);
  
    double deta_dx, deta_dy;
    double u_east,u_north;
    eta = Inc.eta(x,y,t,&deta_dx,&deta_dy,&u_east,&u_north);
    double deta0_dx = confirm_deta_dx(A,w,k,phase,beta,x,y,t);
    double deta0_dy = confirm_deta_dy(A,w,k,phase,beta,x,y,t);
    double u0_east = confirm_ueast(A,w,k,phase,beta,Inc.m_grav,x,y,t);
    double u0_north = confirm_unorth(A,w,k,phase,beta,Inc.m_grav,x,y,t);
    EXPECT_DOUBLE_EQ(eta,eta0);
    EXPECT_DOUBLE_EQ(deta_dx,deta0_dx);
    EXPECT_DOUBLE_EQ(deta_dy,deta0_dy);
    EXPECT_DOUBLE_EQ(u_east,u0_east);
    EXPECT_DOUBLE_EQ(u_north,u0_north);
 
    eta = Inc.eta(x,y,t,&deta_dx,&deta_dy);
    deta0_dx = confirm_deta_dx(A,w,k,phase,beta,x,y,t);
    deta0_dy = confirm_deta_dy(A,w,k,phase,beta,x,y,t);
    EXPECT_DOUBLE_EQ(eta,eta0);
    EXPECT_DOUBLE_EQ(deta_dx,deta0_dx);
    EXPECT_DOUBLE_EQ(deta_dy,deta0_dy);
    }
}

TEST(MonoChrome, CheckVelocities) {
EXPECT_EQ(1+1+1,3);
}

