close all
clear all


addpath('../../matlab');
A = 1;
T = 12;
phase = 0 * pi / 180.0;
beta = 180 * pi / 180.0;
seed = 0;

Inc = LinearIncidentWave();

Inc.SetSeed(seed);
Inc.SetToMonoChromatic(A,T,phase,beta);
Inc.SetToMonoChromatic(A,.3*T,phase,beta);

x = 0; y = 0;

dt = .01;
t = 0:dt:2*T;

eta = NaN(size(t));
eta_0 = NaN(size(t));
eta_1 = NaN(size(t));

for i = 1:length(t)
    deta_dx = 0;
    deta_dy = 0;
  eta(i) = Inc.eta(x,y,t(i),deta_dx,deta_dy);
  deta_dx
  eta_0(i) = Inc.eta(x,y,t(i),0);
  eta_1(i) = Inc.eta(x,y,t(i),1);

end

figure
plot(t,eta);
hold all;
plot(t,eta_0);
plot(t,eta_1);
grid on;

xlabel('time (s)');
ylabel('eta (m)');


