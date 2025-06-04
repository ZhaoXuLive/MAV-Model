clear 
close all

%% utils
addpath utils

deg = pi/180;
kph = 1/3.6; 
g = 9.81; 
small = 1e-3; 
tiny = 1e-6;

%% set simulation conditions
duration = 100;
steercase = 4; 
vconst = 20*kph; 

%% ART data
Q_ART_par

%% normally don't need to edit the following
v_ref = vconst;
sim_model = 'manual_steer'; 
[tstr_in, str_in] = get_input(steercase);

%% driver model parameters
Cparam14  

KPspd = 20;
KIspd = 0;
KDspd = 0;

%% Simulation 
set_param('manual_steer','AlgebraicLoopSolver','Auto');
sim(sim_model, duration)
disp('simulation finished')

% post_proc_qqd
% time_final_proc

%% plot result

figure(1)
plot(ans.time(1:2000,1),ans.steer_out(1:2000,1) * 180 / pi)
title('Vehicle Steering Angle')
xlabel('t (s)')
ylabel('steering angle (deg)')
legend("steer angle");

% plot(ans.time(1:2000,1)_single(1:2000,1),ans.steer_out_single(1:2000,1) * 180 / pi)
% plot(ans.time(1:2000,1)_single_littles(1:2000,1),ans.steer_out_single_littles(1:2000,1) * 180 / pi)
% plot(ans.time(1:2000,1)_newton(1:2000,1),ans.steer_out_newton(1:2000,1) * 180 / pi)

figure(2)
x_newton = ans.q_out_newton(1, :, 1:2000);
x_newton = x_newton(:);
y_newton = ans.q_out_newton(2, :, 1:2000);
y_newton = y_newton(:);
plot(x_newton, y_newton, 'y')
hold on
x_kinetic = ans.q_out_kinetic(1, :, 1:2000);
x_kinetic = x_kinetic(:);
y_kinetic = ans.q_out_kinetic(2, :, 1:2000);
y_kinetic = y_kinetic(:);
plot(x_kinetic, y_kinetic, 'k');
hold on
plot(ans.q_out_single_littles(1:2000,1), ans.q_out_single_littles(1:2000,2), 'b--')
hold on
plot(ans.q_out_single(1:2000,1), ans.q_out_single(1:2000,2), 'r--')
hold on
plot(ans.q_out(1:2000,1), ans.q_out(1:2000,2), 'g--')
title('Vehicle Trajectory')
xlabel('X (m)')
ylabel('Y (m)')
legend("newton", "kinetic", "lagrange single little steer", "lagrange single", "lagrange");

figure(3)
vx_newton = ans.qd_out_newton(1, :, 1:2000);
vx_newton = vx_newton(:);
plot(ans.time(1:2000,1), vx_newton, 'y')
hold on
vx_kinetic = ans.qd_out_kinetic(1, :, 1:2000);
vx_kinetic = vx_kinetic(:);
plot(ans.time(1:2000,1), vx_kinetic, 'k');
hold on
plot(ans.time(1:2000,1), ans.qd_out_single_littles(1:2000, 1), 'b--')
hold on
plot(ans.time(1:2000,1), ans.qd_out_single(1:2000, 1), 'r--')
hold on
plot(ans.time(1:2000,1), ans.qd_out(1:2000, 1), 'g--')
title('Vehicle longitudinal velocity')
xlabel('T (s)')
ylabel('Vx (m/s)')
legend("newton", "kinetic", "lagrange single little steer", "lagrange single", "lagrange");

figure(4)
vy_newton = ans.qd_out_newton(2, :, 1:2000);
vy_newton = vy_newton(:);
plot(ans.time(1:2000,1), vy_newton, 'y')
hold on
vy_kinetic = ans.qd_out_kinetic(2, :, 1:2000);
vy_kinetic = vy_kinetic(:);
plot(ans.time(1:2000,1), vy_kinetic, 'k');
hold on
plot(ans.time(1:2000,1), ans.qd_out_single_littles(1:2000, 2), 'b--')
hold on
plot(ans.time(1:2000,1), ans.qd_out_single(1:2000, 2), 'r--')
hold on
plot(ans.time(1:2000,1), ans.qd_out(1:2000, 2), 'g--')
title('Vehicle lateral velocity')
xlabel('T (s)')
ylabel('Vx (m/s)')
legend("newton", "kinetic", "lagrange single little steer", "lagrange single", "lagrange");

