
clear 
close all

%% utils
addpath utils
deg=pi/180; kph=1/3.6;g=9.81; 
small=1e-3;tiny=1e-6;

%% Track definition
TRACK=trackgen4(2,50);
track=TRACK; %used in results file
vconst=20*kph; %alternatively prescribe set speed

%% Low level control param 
Tmc=0.01;
Td=0.01; %driving actuator time constant
Tb=0.002; %braking time constant
brake_lag=0.002; %pure delay on brakes
brake_rate=5000; %Nm/sec rate limit on each brake channel
Tst=0.02; %steer actuator time constant
RW=0.499;%wheel rolling radius (m)

%apportionment of drive and brake torque
DRIVE=[0.25 0.25  0 0 0 0 0.25 0.25]/1';
BRAKE=[1 1 1 1 1 1 1 1]/8';

%% ART data
Q_ART_par

%% driver model parameters
Cparam14  

% KPspd = 20;
% KIspd = 2;
% KDspd = 2;

KPspd = 20;
KIspd = 0;
KDspd = 0;

% Q0=zeros(Nu+2,1);
% QD0=Q0;
% QD0(1)=vconst;

% sign = [0, 0];


% %% SIM
% sim_model = 'AFG_s1v1.mdl'; 
% sim(sim_model);
% 
% return
%% post proc 
% post_proc_AFGs1_debug;