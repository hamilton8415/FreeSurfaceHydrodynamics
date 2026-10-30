close all
clear all
clear mex

addpath('../build');

H = 2.2;
T = 12;
beta = 30 * pi / 180.0;

Inc = LinearIncidentWave();

%Inc.SetSeed(seed);


%Inc.SetToPiersonMoskowitzSpectrum(H,beta);
Inc.SetToBretschneiderSpectrum(H,T,beta);

x = 0; 
y = 0;


dt = .2;
t = 0:dt:30*60;
eta = NaN(size(t));


deta_dx = NaN(size(t));
deta_dy = NaN(size(t));

vel_x = NaN(size(t));
vel_y = NaN(size(t));


for i = 1:length(t)
  [eta(i) deta_dx(i) deta_dy(i) vel_x(i) vel_y(i)] = Inc.eta(x,y,t(i));
end

figure
ah1 = subplot(2,1,1);
plot(t,eta);
grid on;
xlabel('time (s)');
ylabel('eta (m)');

ah2 = subplot(2,1,2);
plot(t,vel_x);
hold all;
plot(t,vel_y);
legend('vel_x','vel_y');
grid on;
xlabel('time (s)');
ylabel('(m/s)');

linkaxes([ah1 ah2],'x');


fs = 1/dt;

[ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check ] = SBGwaves(vel_x,vel_y,-eta,fs);

disp(["Hs_input = " num2str(H) "  Hs_output = " num2str(Hs)]);
disp(["Tp_input = " num2str(T) "  Tp_output = " num2str(Tp)]);
disp(["Dp_input = " num2str(beta*180/pi) "  Dp_output = " num2str(Dp)]);


figure
plot(f,E);