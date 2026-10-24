close all
clear all
clear mex

%setenv('DYLD_LIBRARY_PATH', '/System/Volumes/Data/opt/FreeSurfaceHydrodynamics/lib');

addpath('../build');

A = 1;
T = 12;
phase = 0 * pi / 180.0;
beta = 180 * pi / 180.0;
seed = 0;

Inc = LinearIncidentWave();

Inc.SetSeed(seed);
%Inc.SetToMonoChromatic(A,T,phase,beta);
%Inc.SetToMonoChromatic(A,.3*T,phase,beta);

%Inc.SetToPiersonMoskowitzSpectrum(2*A,beta);
%Inc.SetToPiersonMoskowitzSpectrumWithCos2Spreading(2*A,beta, 34.0, 10);

%Inc.SetToBretschneiderSpectrum(2*A,T,beta);
Inc.SetToBretschneiderSpectrumWithCos2Spreading(2*A,T,beta, 50, 34.0, 10);

Inc.PrintWaveInfo()

x = 0; y = 0;

dt = .01;
t = 0:dt:2*T;
eta = NaN(size(t));
deta_dx = NaN(size(t));
deta_dy = NaN(size(t));

eta_0 = NaN(size(t));
eta_1 = NaN(size(t));

for i = 1:length(t)
  [eta(i) deta_dx(i) deta_dy(i)] = Inc.eta(x,y,t(i));
  eta_0(i) = Inc.eta(x,y,t(i),0);
  eta_1(i) = Inc.eta(x,y,t(i),1);

end

figure
ah1 = subplot(2,1,1);
plot(t,eta);
hold all;
plot(t,eta_0);
plot(t,eta_1);
grid on;
xlabel('time (s)');
ylabel('eta (m)');

ah2 = subplot(2,1,2);
plot(t,deta_dx);
hold all;
plot(t,deta_dy);
legend('deta_dx','deta_dy');
%plot(t,eta_0);
%plot(t,eta_1);
grid on;
xlabel('time (s)');
ylabel('(m/m)');

linkaxes([ah1 ah2],'x');


etadot = NaN(size(t));

etadot_0 = NaN(size(t));
etadot_1 = NaN(size(t));

for i = 1:length(t)
  etadot(i) = Inc.etadot(x,y,t(i));
%  etadot_0(i) = Inc.etadot(x,y,t(i),0);
%  etadot_1(i) = Inc.etadot(x,y,t(i),1);
end

figure
plot(t,etadot);
hold all;
plot(t,etadot_0);
plot(t,etadot_1);
grid on;
xlabel('time (s)');
ylabel('etadot (m/s)');


xx = -100:.2:100;
yy = -100:.2:100;
dt = .1;

[X, Y] = meshgrid(-xx, yy);

Z = NaN(size(X));

t = 3;

Z = NaN(length(xx),length(yy));

    for m = 1:length(xx)
        for n = 1:length(yy)
            Z(m,n) = Inc.eta(X(m,n),Y(m,n),t);
        end
    end

figure
s = surf(X, Y, Z);
s.EdgeColor = "none";
