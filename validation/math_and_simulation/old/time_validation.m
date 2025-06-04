
clear 
close all

%% utils
addpath utils

deg=pi/180; kph=1/3.6; g=9.81; 
small=1e-3; tiny=1e-6;

%% Track definition

% 如果有track则需要在这里选择track，否则应置为空
% TRACK=trackgen4(2,50);
TRACK=[];
track=TRACK; %used in results file

%% set simulation conditions

duration = 100;
steercase = 4; %  1=step 4=sinus. See open_loop_steer.m for other options
vconst = 20*kph; %alternatively prescribe set speed

%% Low level control param 
% Tmc=0.01;
% Td=0.01; %driving actuator time constant
% Tb=0.002; %braking time constant
% brake_lag=0.002; %pure delay on brakes
% brake_rate=5000; %Nm/sec rate limit on each brake channel
% Tst=0.02; %steer actuator time constant
% RW=0.499;%wheel rolling radius (m)
% 
% %apportionment of drive and brake torque
% DRIVE=[0.25 0.25  0 0 0 0 0.25 0.25]/1';
% BRAKE=[1 1 1 1 1 1 1 1]/8';

%% ART data
Q_ART_par

%% normally don't need to edit the following
v_ref = vconst;
sim_model = 't_validation_V2'; %small mod in version 02,  A\ instead of inv(A)
%open-loop steer, no control on following axles
tstr_in = (0:0.01:100)';

% important
[tref,swaref] = open_loop_steer(steercase); 
str_in = interp1(tref, swaref, tstr_in,'linear','extrap');

%% driver model parameters
Cparam14  

% KPspd = 20;
% KIspd = 2;
% KDspd = 2;

KPspd = 20;
KIspd = 0;
KDspd = 0;

%% Simulation (with initial conditions)
Q0 = zeros(Nu+2,1);
QD0 = Q0;
QD0(1) = v_ref;
sim(sim_model, duration)
disp('simulation finished')

% post_proc_qqd
time_final_proc

return
