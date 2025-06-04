
% clear 
close all

%% utils
addpath utils
deg=pi/180; kph=1/3.6;g=9.81; 
small=1e-3;tiny=1e-6;
cutoff_freq = 60;

%% Track definition
TRACK=trackgen4(2,50);
% TRACK=trackgen4(3,30);
% load('smooth.mat');
% points = [x_smooth,y_smooth];
% load('originMap.mat');
% points = [x_ref,y_ref];
% TRACK = pts2trk_v2(points);
% load('track_circle15.mat');
% load('track_straight_2km.mat');
% load('track_angle90turn_r25.mat');
% load('track_angle90turn_r15.mat');
% load('track_double_lane_change.mat');
% load('track_single_lane_change.mat');

% TRACK = temptrack;
track=TRACK; %used in results file
% vconst=70*kph; %alternatively prescribe set speed
vconst=15*kph;
%% Low level control param 
% Tmc=0.01;
Tmc=0.02;
% Td=0.01; %driving actuator time constant
% Tb=0.002; %braking time constant
% brake_lag=0.002; %pure delay on brakes
% brake_rate=5000; %Nm/sec rate limit on each brake channel
Tst=0.05; %steer actuator time constant
RW=0.499;%wheel rolling radius (m)

%apportionment of drive and brake torque
DRIVE=[0.25 0.25  0 0 0 0 0.25 0.25]/1';
BRAKE=[1 1 1 1 1 1 1 1]/8';

%% ART data
Q_ART_par

%% driver model parameters
% Cparam14  

% KPspd = 20;
% KIspd = 2;
% KDspd = 2;
Nspd = 10;
Uspd = 5;
Lspd = -7;

axle_rate_limit=1; %rad/sec at all 6 axles

KPspd = 10;
KIspd = 0.1;
KDspd = 0.1;
% KPspd.CoderInfo.StorageClass = 'ExportedGlobal' ;

% KPspd = 0.5;
% KIspd = 0;
% KDspd = 0;


Kap=0.1;Kai=0.05;Kad=0.0;
% Kap=0.3;
% Kai=0.1;
% Kad=0;


% vref_mat=[ones(1,1000),linspace(1,vconst,3000),ones(1,1000)*vconst,linspace(vconst,0,3000),zeros(1,100000)];
% vref_mat=[linspace(20/3.6,vconst,4861),ones(1,1000)*vconst,linspace(vconst,0,4861),zeros(1,100000)];

% vref_mat=[ones(1,1000)/36,linspace(1/36,vconst,4861),ones(1,1000)*vconst,linspace(vconst,0,4861),zeros(1,100000)];
vref_mat=[zeros(1,1000),linspace(0,vconst,4861),ones(1,1000)*vconst,linspace(vconst,0,4861),zeros(1,100000)];


%% State Estimate  (Params needs to tune)
% %%%%%%%%%%%noise psd%%%%%%%%%%%%%%%
% 
% % PSD_PX=1e-4/2*10^(50/20);
% % PSD_PY=1e-4/2*10^(50/20);
% % PSD_VX=4e-5/2*10^(5/20);
% % PSD_VY=4e-5/2*10^(5/20);
% % % PSD_ARTIC_1=4e-5/2*deg^2;
% % % PSD_ARTIC_2=4e-5/2*deg^2*10^(10/20);
% % PSD_ARTIC_1=4e-5/2*deg^2*10^(10/20);
% % PSD_ARTIC_2=4e-5/2*deg^2*10^(10/20);
% % PSD_PSI=1e-8/2*10^(5/20);
% % PSD_R=5e-8/2;PSD_R=PSD_R*10^(-5/20);
% 
% %KF TUNING
% 
% % %previous values for process noise
% % Q1=1e-8*eye(10);
% % Q1(1,1)=1e-8;Q1(2,2)=Q1(1,1);
% % Q1(3,3)=1e-9;Q1(4,4)=Q1(3,3);Q1(5,5)=Q1(3,3);
% % Q1(6,6)=5e-5;Q1(7,7)=Q1(6,6);
% % Q1(8,8)=5e-7;Q1(9,9)=Q1(8,8);Q1(10,10)=Q1(8,8);
% 
% q=1e-4*ones(1,10); %default for rms process noise
% q(1:2)=q(1:2) * 3;
% q(3:5)=q(3:5)/5; %yaw angle kinematics
% q(6:7)=q(6:7)*10; %vx vy velocities
% q(8:10)=q(8:10)*2; %yaw rate equations
% 
% Q1=diag(q.^2);
% Q2=Q1;
% Q3=Q1;
% 
% % R1=diag([0.0456,0.0456,6.7708e-08,2.2371e-04,2.2371e-04,6.0907e-06,6.0907e-06,8.7317e-08,8.7317e-08]);%x,y,psi,vx,vy,r,r,theta,theta
% %set R1 based on set noise levels in simulink model
% 
% % psd_all=[PSD_PX,PSD_PY,PSD_PSI,PSD_VX,PSD_VY,PSD_R,PSD_R,PSD_ARTIC_1,PSD_ARTIC_2];
% % bw=[5 5 10 10 10 50 50 10 10];
% psd_all=[0.00006, 0.00006, 0.00001,0,0,0,0,0,0];
% bw=[1 1 1 1 1 1 1 1 1];
% R1=diag(psd_all.*bw);
% R2=R1;
% R3=diag([0.0456,0.0456,6.7708e-08,2.2371e-04,2.2371e-04,6.0907e-06,8.7317e-08,8.7317e-08]);%x,y,psi,vx,vy,r,theta,theta
% 
% % Q0 = [0 0 0 0 0];
% Q0 = [32.2 0 0 0 0];
% QD0 = [0 0 0 0 0];
% zz_gps = [36.334837 120.2745019];
% % Dgps = 0.5;
% Dgps = 11.4;

%% Test data load
% load('data_qqd.mat');
% load('test4SE_J2.mat');

% %% SIM
% sim_model = 'AFG_s1v1.mdl'; 
% sim(sim_model);
% 
% return
%% post proc 
% post_proc_AFGs1_debug;