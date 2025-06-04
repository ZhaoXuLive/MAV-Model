%Cparam
%driver model parameters

Tmc=0.01; %Motion Controller subsystem time step


%speed controller
KPspd = 10;
KIspd = 0.5;
KDspd = 1;
Nspd = 10;
Uspd = 5;
Lspd = -7;
maxspeed=30; %SET NEGATIVE FOR CONSTANT REFEENCE SPEED
ayref=2; %max lateral accel tolerance (approx) for speed planner

%driver path angle feedback (higher level controller for path angle error to yaw rate correction)
% KPpae=5; %when flow acc absent
KPpae=2;
KIpae=0;
KDpae=0;
SATpae=1.0;%output saturation for PID control



%driver steering control (version 13 => PID control)
KPstr = 40;
KIstr = 10;
KDstr = 0;
Nstr = 10;
Ustr = 270*deg;
Lstr = -270*deg;


SWAmax=180*deg; %rad at steering wheel
SWRateMax=10; %rad/sec at steering wheel


%axle feedback steering control (per unit angle error)
% Kap=0.25;Kai=0;Kad=0.0;
Kap=0.1;Kai=0.05;Kad=0.0;
% Kap=0.0;Kai=0.0;Kad=0.0;

axle_rate_limit=1.0; %rad/sec at all 6 axles


