
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

duration = 500;
vconst = 20*kph; %alternatively prescribe set speed

%% ART data
Q_ART_par

%% normally 
v_ref = vconst;
sim_model = 'f_validation_V1'; 
% tstr_in = (0:0.01:1000)';

% important
f = 0.05;
a = 3;
[tref,swaref] = open_loop_steer(f, a); 
% str_in = interp1(tref, swaref, tstr_in,'linear','extrap');
str_in = swaref';
tstr_in = tref';
duration = length(str_in) / 100;

%% driver model parameters
Cparam14  

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
