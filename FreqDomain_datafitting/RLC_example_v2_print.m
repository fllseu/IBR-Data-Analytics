%% RLC parameter specification
clear; clc;
R = 0.1; L = 0.5/377; C = 1/0.2/377;

%% create frequency-domain data
f = logspace(-1, 3, 30); s = 1i*2*pi*f; 
sys = 1./(R*C*s+L*C*s.^2+ 1); 

%% add noise
[n,m ]= size(sys); sys1 = (1+ (rand(n,m)-0.5)*0.2).*sys;

%% tfest
Measurement_F = idfrd(sys1,f*2*pi,0); 
Estimation_F = tfest(Measurement_F, 2);

%% 

%% Bode plots
figure(1000)
bode(Measurement_F, Estimation_F); hold on; 
legend('Measurement', 'Estimation');
xlim([0.1, 1000]); 

