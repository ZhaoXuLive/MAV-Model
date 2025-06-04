close all

%% the vehicle parameters

% the number of vehicle units
unit_n = 3;
% the number of vehicle axles
axle_m = 8;
% the number of axles of every vehicle unit 
getNumAxle = [0, 3, 5, 8];
% the unit of the respective axle
getNumUnit = [1, 1, 1, 2, 2, 3, 3, 3];

% the mass of every vehicle unit
m = [11142, 7651, 11142];
% the rotational inertia of every vehicle unit
Iz = [143570, 37800, 143570];
% the tire cornering stiffness
Ca = 3065 * 180 / pi;
% the articulated angle damping
art_d = [10000, 10000];
% art_d = [0, 0];

% the length of every vehicle unit
lu = [11.4, 9.4, 11.4];
% the length of every vehicle unit from the front joint to the mass center
lg = [5.862, 4.643, 5.538];
% the length of every vehicle unit from the front joint to the axle center
la = [2.5, 4, 9.9, 1.5, 7.9, 1.5, 7.4, 8.9];
% the lateral length of every axle from the axle centor to the tyre 
lt = [0.9625, 0.9625, 0.98, 0.98, 0.98, 0.98, 0.9625, 0.9625];
deg = pi / 180;

sat = 10;

%% real vehicle data

% 3 无意�?
% 4 直线 运动学最�?
% 1267 LS > LSL > L > N > K
% 5 L > LS > LSL > N > K

data_mode = 1;

% 1 single 2 multi
mode = 2;

dt = 0.01;
if data_mode == 1
    load('quarter_turn.mat')
    load('quarter_turn_steer.mat')
    startpoint = 5000;
    endpoint = 7000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 2
    load('circle.mat')
    load('circle_steer.mat')
    startpoint = 26000;
    endpoint = 28000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 3
    load('double_change.mat')
    load('double_change_steer.mat')
    startpoint = 1200;
    endpoint = 5000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 4
    load('2318.mat')
    startpoint = 4000;
    endpoint = 6000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 5
    load('2122.mat')
    startpoint = 4000;
    endpoint = 6000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 6
    load('2146.mat')
    startpoint = 2000;
    endpoint = 4000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 7
    load('2148.mat')
    startpoint = 3000;
    endpoint = 5000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end

if mode == 1
    qd_lagrange_single_littles = zeros(data_len, unit_n+2);
    q_lagrange_single_littles = zeros(data_len, unit_n+2);
    qd_lagrange_single = zeros(data_len, unit_n+2);
    q_lagrange_single = zeros(data_len, unit_n+2);
    qd_lagrange = zeros(data_len, unit_n+2);
    q_lagrange = zeros(data_len, unit_n+2);
    qd_newton = zeros(data_len, unit_n+2);
    q_newton = zeros(data_len, unit_n+2);
    qd_kinetic = zeros(data_len, unit_n+2);
    q_kinetic = zeros(data_len, unit_n+2);
end
if mode == 2
    qd_lagrange_single_littles_multi = zeros(data_len, unit_n+2);
    q_lagrange_single_littles_multi = zeros(data_len, unit_n+2);
    qd_lagrange_single_multi = zeros(data_len, unit_n+2);
    q_lagrange_single_multi = zeros(data_len, unit_n+2);
    qd_lagrange_multi = zeros(data_len, unit_n+2);
    q_lagrange_multi = zeros(data_len, unit_n+2);
    qd_newton_multi = zeros(data_len, unit_n+2);
    q_newton_multi = zeros(data_len, unit_n+2);
    qd_kinetic_multi = zeros(data_len, unit_n+2);
    q_kinetic_multi = zeros(data_len, unit_n+2);
end

for i = 1:data_len
    q_i = q_debug(i, :)';
    qd_i = qd_debug(i, :)';
    delta_i(1) = steer(i, 1);
    delta_i(2) = 0.7 * steer(i, 1);
    delta_i(3:6) = steer(i, 2:5)';
    delta_i(7) = 0.7 * steer(i, 6);
    delta_i(8) = steer(i, 6);

    if mode == 1
        % lagrange nonlinear double 
        qdd_lagrange = lagrange_fix_model(q_i, qd_i, delta_i, sat);
        qd_lagrange(i, :) = qdd_lagrange * dt + qd_i;
        q_lagrange(i, :) = (qd_lagrange(i, :))' * dt + q_i; 

        % lagrange nonliear single
        qdd_lagrange_single = lagrange_single_fix_model(q_i, qd_i, delta_i, sat);
        qd_lagrange_single(i, :) = qdd_lagrange_single * dt + qd_i;
        q_lagrange_single(i, :) = (qd_lagrange_single(i, :))' * dt + q_i; 

        % lagrange nonliear single little steer
        qdd_lagrange_single_littles = lagrange_littlesteer_single_fix_model(q_i, qd_i, delta_i, sat);
        qd_lagrange_single_littles(i, :) = qdd_lagrange_single_littles * dt + qd_i;
        q_lagrange_single_littles(i, :) = (qd_lagrange_single_littles(i, :))' * dt + q_i; 

        % newton linear 
        qdd_newton = newton_model(unit_n, axle_m, ..., 
        getNumUnit, m, Iz, Ca, lu, lg, la, q_i, qd_i, delta_i, sat);
        qd_newton(i, :) = qdd_newton * dt + qd_i;
        q_newton(i, :) = (qd_newton(i, :))' * dt + q_i;

        % kinetic
        vxg = qd_i(1) + lg(1) * sin(q_i(3)) * qd_i(3);
        vyg = qd_i(2) - lg(1) * cos(q_i(3)) * qd_i(3);
        qd_kinetic(i, :) = kinetic_model(q_i, vxg, vyg, delta_i);
        q_kinetic(i, :) = (qd_kinetic(i, :))' * dt + q_i;
        qdd_lagrange = lagrange_fix_model(q_i, qd_i, delta_i, sat);
    end
    if mode == 2
        if i == 1
            qd_lagrange_single_littles_multi(i, :) = qd_i;
            q_lagrange_single_littles_multi(i, :) = q_i;
            qd_lagrange_single_multi(i, :) = qd_i;
            q_lagrange_single_multi(i, :) = q_i;
            qd_lagrange_multi(i, :) = qd_i;
            q_lagrange_multi(i, :) = q_i;
            qd_newton_multi(i, :) = qd_i();
            q_newton_multi(i, :) = q_i();
            qd_kinetic_multi(i, :) = qd_i();
            q_kinetic_multi(i, :) = q_i();
        else
            % lagrange nonlinear double 
            qdd_lagrange = lagrange_fix_model(q_lagrange_multi(i-1, :)', qd_lagrange_multi(i-1, :)', delta_i, sat);
            qd_lagrange_multi(i, :) = qdd_lagrange * dt + qd_lagrange_multi(i-1, :)';
            q_lagrange_multi(i, :) = (qd_lagrange_multi(i, :))' * dt + q_lagrange_multi(i-1, :)';

            % lagrange nonliear single
            qdd_lagrange_single = lagrange_single_fix_model(q_lagrange_multi(i-1, :)', qd_lagrange_multi(i-1, :)', delta_i, sat);
            qd_lagrange_single_multi(i, :) = qdd_lagrange_single * dt + qd_lagrange_single_multi(i-1, :)';
            q_lagrange_single_multi(i, :) = (qd_lagrange_single_multi(i, :))' * dt + q_lagrange_single_multi(i-1, :)';

            % lagrange nonliear single little steer
            qdd_lagrange_single_littles = lagrange_single_fix_model(q_lagrange_multi(i-1, :)', qd_lagrange_multi(i-1, :)', delta_i, sat);
            qd_lagrange_single_littles_multi(i, :) = qdd_lagrange_single_littles * dt + qd_lagrange_single_littles_multi(i-1, :)';
            q_lagrange_single_littles_multi(i, :) = (qd_lagrange_single_littles_multi(i, :))' * dt + q_lagrange_single_littles_multi(i-1, :)';

            % newton linear 
            qdd_newton = newton_model(unit_n, axle_m, ..., 
                    getNumUnit, m, Iz, Ca, lu, lg, la, q_newton_multi(i-1, :)', qd_newton_multi(i-1, :)', delta_i, sat);
            qd_newton_multi(i, :) = qdd_newton * dt + qd_newton_multi(i-1, :)';
            q_newton_multi(i, :) = (qd_newton_multi(i, :))' * dt + q_newton_multi(i-1, :)';

            % kinetic
            vxg = qd_debug(1, 1) + lg(1) * sin(q_debug(1, 3)) * qd_debug(1, 3);
            vyg = qd_debug(1, 2) - lg(1) * cos(q_debug(1, 3)) * qd_debug(1, 3);
            qd_kinetic_multi(i, :) = kinetic_model(q_kinetic_multi(i-1, :)', vxg, vyg, delta_i);
            q_kinetic_multi(i, :) = (qd_kinetic_multi(i, :))' * dt + q_kinetic_multi(i-1, :)';
        end
    end
end

%% mode == 1
if mode == 1
    figure(1)
    plot(time(1:end-1), qd_lagrange(1:end-1, 1), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single(1:end-1, 1), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles(1:end-1, 1), 'm--');
    hold on
    plot(time(1:end-1), qd_newton(1:end-1, 1), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic(1:end-1, 1), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 1));
    axis equal
    
    title('Vehicle x speed compare')
    xlabel('t (s)')
    ylabel('Y (m/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(2)
    plot(time(1:end-1), qd_lagrange(1:end-1, 2), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single(1:end-1, 2), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles(1:end-1, 2), 'm--');
    hold on
    plot(time(1:end-1), qd_newton(1:end-1, 2), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic(1:end-1, 2), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 2));
    axis equal
    
    title('Vehicle y speed compare')
    xlabel('t (s)')
    ylabel('Y (m/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(3)
    plot(time(1:end-1), qd_lagrange(1:end-1, 3), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single(1:end-1, 3), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles(1:end-1, 3), 'm--');
    hold on
    plot(time(1:end-1), qd_newton(1:end-1, 3), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic(1:end-1, 3), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 3));
    axis equal
    
    title('Vehicle yawrate1 compare')
    xlabel('t (s)')
    ylabel('yawrate1 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(4)
    plot(time(1:end-1), qd_lagrange(1:end-1, 4), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single(1:end-1, 4), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles(1:end-1, 4), 'm--');
    hold on
    plot(time(1:end-1), qd_newton(1:end-1, 4), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic(1:end-1, 4), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 4));
    axis equal
    
    title('Vehicle yawrate2 compare')
    xlabel('t (s)')
    ylabel('yawrate2 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(5)
    plot(time(1:end-1), qd_lagrange(1:end-1, 5), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single(1:end-1, 5), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles(1:end-1, 5), 'm--');
    hold on
    plot(time(1:end-1), qd_newton(1:end-1, 5), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic(1:end-1, 5), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 5));
    axis equal
    
    title('Vehicle yawrate3 compare')
    xlabel('t (s)')
    ylabel('yawrate3 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(6)
    plot(q_lagrange(1:end-1, 1), q_lagrange(1:end-1, 2), 'g--');
    hold on
    plot(q_lagrange_single(1:end-1, 1), q_lagrange_single(1:end-1, 2), 'y--');
    hold on
    plot(q_lagrange_single_littles(1:end-1, 1), q_lagrange_single_littles(1:end-1, 2), 'm--');
    hold on
    plot(q_newton(1:end-1, 1), q_newton(1:end-1, 2), 'b--');
    hold on
    plot(q_kinetic(1:end-1, 1), q_kinetic(1:end-1, 2), 'r--');
    hold on
    plot(q_debug(2:data_len, 1), q_debug(2:data_len, 2));
    axis equal
    
    title('Vehicle trajectory compare')
    xlabel('X(m)')
    ylabel('Y(m)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    MSE_L = 0;
    MSE_LS = 0;
    MSE_LSL = 0;
    MSE_N = 0;
    MSE_K = 0;
    RMSE_L = 0;
    RMSE_LS = 0;
    RMSE_LSL = 0;
    RMSE_N = 0;
    RMSE_K = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data_len_l = data_len-1;
    data_len_ls = data_len-1;
    data_len_lsl = data_len-1;
    data_len_n = data_len-1;
    data_len_k = data_len-1;
    for i = 1:data_len-1
        if isnan(q_lagrange(i, 1)) || isnan(q_lagrange(i, 2))
            data_len_l = data_len_l - 1;
        else
            RMSE_L = RMSE_L + hypot(q_debug(i+1, 1) - q_lagrange(i, 1), q_debug(i+1, 2) - q_lagrange(i, 2))^2;
            MSE_L = MSE_L + hypot(q_debug(i+1, 1) - q_lagrange(i, 1), q_debug(i+1, 2) - q_lagrange(i, 2)); 
        end
        if isnan(q_lagrange_single(i, 1)) || isnan(q_lagrange_single(i, 2))
            data_len_ls = data_len_ls - 1;
        else
            RMSE_LS = RMSE_LS + hypot(q_debug(i+1, 1) - q_lagrange_single(i, 1), q_debug(i+1, 2) - q_lagrange_single(i, 2))^2;
            MSE_LS = MSE_LS + hypot(q_debug(i+1, 1) - q_lagrange_single(i, 1), q_debug(i+1, 2) - q_lagrange_single(i, 2)); 
        end
        if isnan(q_lagrange_single_littles(i, 1)) || isnan(q_lagrange_single_littles(i, 2))
            data_len_lsl = data_len_lsl - 1;
        else
            RMSE_LSL = RMSE_LSL + hypot(q_debug(i+1, 1) - q_lagrange_single_littles(i, 1), q_debug(i+1, 2) - q_lagrange_single_littles(i, 2))^2;
            MSE_LSL = MSE_LSL + hypot(q_debug(i+1, 1) - q_lagrange_single_littles(i, 1), q_debug(i+1, 2) - q_lagrange_single_littles(i, 2)); 
        end
        if isnan(q_newton(i, 1)) || isnan(q_newton(i, 2))
            data_len_n = data_len_n - 1;
        else
            RMSE_N = RMSE_N + hypot(q_debug(i+1, 1) - q_newton(i, 1), q_debug(i+1, 2) - q_newton(i, 2))^2;
            MSE_N = MSE_N + hypot(q_debug(i+1, 1) - q_newton(i, 1), q_debug(i+1, 2) - q_newton(i, 2)); 
        end
        if isnan(q_kinetic(i, 1)) || isnan(q_kinetic(i, 2))
            data_len_k = data_len_k - 1;
        else
            RMSE_K = RMSE_K + hypot(q_debug(i+1, 1) - q_kinetic(i, 1), q_debug(i+1, 2) - q_kinetic(i, 2))^2;
            MSE_K = MSE_K + hypot(q_debug(i+1, 1) - q_kinetic(i, 1), q_debug(i+1, 2) - q_kinetic(i, 2)); 
        end
    end
    MSE_L = MSE_L / data_len_l
    MSE_LS = MSE_LS / data_len_ls
    MSE_LSL = MSE_LSL / data_len_lsl
    MSE_N = MSE_N / data_len_n
    MSE_K = MSE_K / data_len_k
    RMSE_L = sqrt(RMSE_L / data_len_l)
    RMSE_LS = sqrt(RMSE_LS / data_len_ls)
    RMSE_LSL = sqrt(RMSE_LSL / data_len_lsl)
    RMSE_N = sqrt(RMSE_N / data_len_n)
    RMSE_K = sqrt(RMSE_K / data_len_k)
end

%% mode == 2
if mode == 2
    figure(1)
    plot(time(1:end-1), qd_lagrange_multi(1:end-1, 1), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_multi(1:end-1, 1), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles_multi(1:end-1, 1), 'm--');
    hold on
    plot(time(1:end-1), qd_newton_multi(1:end-1, 1), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic_multi(1:end-1, 1), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 1));
    axis equal
    
    title('Vehicle x speed compare')
    xlabel('t (s)')
    ylabel('Y (m/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(2)
    plot(time(1:end-1), qd_lagrange_multi(1:end-1, 2), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_multi(1:end-1, 2), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles_multi(1:end-1, 2), 'm--');
    hold on
    plot(time(1:end-1), qd_newton_multi(1:end-1, 2), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic_multi(1:end-1, 2), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 2));
    axis equal
    
    title('Vehicle y speed compare')
    xlabel('t (s)')
    ylabel('Y (m/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(3)
    plot(time(1:end-1), qd_lagrange_multi(1:end-1, 3), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_multi(1:end-1, 3), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles_multi(1:end-1, 3), 'm--');
    hold on
    plot(time(1:end-1), qd_newton_multi(1:end-1, 3), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic_multi(1:end-1, 3), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 3));
    axis equal
    
    title('Vehicle yawrate1 compare')
    xlabel('t (s)')
    ylabel('yawrate1 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(4)
    plot(time(1:end-1), qd_lagrange_multi(1:end-1, 4), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_multi(1:end-1, 4), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles_multi(1:end-1, 4), 'm--');
    hold on
    plot(time(1:end-1), qd_newton_multi(1:end-1, 4), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic_multi(1:end-1, 4), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 4));
    axis equal
    
    title('Vehicle yawrate2 compare')
    xlabel('t (s)')
    ylabel('yawrate2 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(5)
    plot(time(1:end-1), qd_lagrange_multi(1:end-1, 5), 'g--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_multi(1:end-1, 5), 'y--');
    hold on
    plot(time(1:end-1), qd_lagrange_single_littles_multi(1:end-1, 5), 'm--');
    hold on
    plot(time(1:end-1), qd_newton_multi(1:end-1, 5), 'b--');
    hold on
    plot(time(1:end-1), qd_kinetic_multi(1:end-1, 5), 'r--');
    hold on
    plot(time(1:end-1), qd_debug(2:data_len, 5));
    axis equal
    
    title('Vehicle yawrate3 compare')
    xlabel('t (s)')
    ylabel('yawrate3 (rad/s)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    figure(6)
    plot(q_lagrange_multi(1:end-1, 1), q_lagrange_multi(1:end-1, 2), 'g--');
    hold on
    plot(q_lagrange_single_multi(1:end-1, 1), q_lagrange_single_multi(1:end-1, 2), 'y--');
    hold on
    plot(q_lagrange_single_littles_multi(1:end-1, 1), q_lagrange_single_littles_multi(1:end-1, 2), 'm--');
    hold on
    plot(q_newton_multi(1:end-1, 1), q_newton_multi(1:end-1, 2), 'b--');
    hold on
    plot(q_kinetic_multi(1:end-1, 1), q_kinetic_multi(1:end-1, 2), 'r--');
    hold on
    plot(q_debug(2:data_len, 1), q_debug(2:data_len, 2));
    axis equal
    
    title('Vehicle trajectory compare')
    xlabel('X(m)')
    ylabel('Y(m)')
    legend("lagrange","lagrange single","lagrange single little steer","newton","kinetic","real");

    MSE_L = 0;
    MSE_LS = 0;
    MSE_LSL = 0;
    MSE_N = 0;
    MSE_K = 0;
    RMSE_L = 0;
    RMSE_LS = 0;
    RMSE_LSL = 0;
    RMSE_N = 0;
    RMSE_K = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data_len_l = data_len-1;
    data_len_ls = data_len-1;
    data_len_lsl = data_len-1;
    data_len_n = data_len-1;
    data_len_k = data_len-1;
    for i = 1:data_len-1
        if isnan(q_lagrange_multi(i, 1)) || isnan(q_lagrange_multi(i, 2))
            data_len_l = data_len_l - 1;
        else
            RMSE_L = RMSE_L + hypot(q_debug(i+1, 1) - q_lagrange_multi(i, 1), q_debug(i+1, 2) - q_lagrange_multi(i, 2))^2;
            MSE_L = MSE_L + hypot(q_debug(i+1, 1) - q_lagrange_multi(i, 1), q_debug(i+1, 2) - q_lagrange_multi(i, 2)); 
        end
        if isnan(q_lagrange_single_multi(i, 1)) || isnan(q_lagrange_single_multi(i, 2))
            data_len_ls = data_len_ls - 1;
        else
            RMSE_LS = RMSE_LS + hypot(q_debug(i+1, 1) - q_lagrange_single_multi(i, 1), q_debug(i+1, 2) - q_lagrange_single_multi(i, 2))^2;
            MSE_LS = MSE_LS + hypot(q_debug(i+1, 1) - q_lagrange_single_multi(i, 1), q_debug(i+1, 2) - q_lagrange_single_multi(i, 2)); 
        end
        if isnan(q_lagrange_single_littles_multi(i, 1)) || isnan(q_lagrange_single_littles_multi(i, 2))
            data_len_lsl = data_len_lsl - 1;
        else
            RMSE_LSL = RMSE_LSL + hypot(q_debug(i+1, 1) - q_lagrange_single_littles_multi(i, 1), q_debug(i+1, 2) - q_lagrange_single_littles_multi(i, 2))^2;
            MSE_LSL = MSE_LSL + hypot(q_debug(i+1, 1) - q_lagrange_single_littles_multi(i, 1), q_debug(i+1, 2) - q_lagrange_single_littles_multi(i, 2)); 
        end
        if isnan(q_newton_multi(i, 1)) || isnan(q_newton_multi(i, 2))
            data_len_n = data_len_n - 1;
        else
            RMSE_N = RMSE_N + hypot(q_debug(i+1, 1) - q_newton_multi(i, 1), q_debug(i+1, 2) - q_newton_multi(i, 2))^2;
            MSE_N = MSE_N + hypot(q_debug(i+1, 1) - q_newton_multi(i, 1), q_debug(i+1, 2) - q_newton_multi(i, 2)); 
        end
        if isnan(q_kinetic_multi(i, 1)) || isnan(q_kinetic_multi(i, 2))
            data_len_k = data_len_k - 1;
        else
            RMSE_K = RMSE_K + hypot(q_debug(i+1, 1) - q_kinetic_multi(i, 1), q_debug(i+1, 2) - q_kinetic_multi(i, 2))^2;
            MSE_K = MSE_K + hypot(q_debug(i+1, 1) - q_kinetic_multi(i, 1), q_debug(i+1, 2) - q_kinetic_multi(i, 2)); 
        end
    end
    MSE_L = MSE_L / data_len_l
    MSE_LS = MSE_LS / data_len_ls
    MSE_LSL = MSE_LSL / data_len_lsl
    MSE_N = MSE_N / data_len_n
    MSE_K = MSE_K / data_len_k
    RMSE_L = sqrt(RMSE_L / data_len_l)
    RMSE_LS = sqrt(RMSE_LS / data_len_ls)
    RMSE_LSL = sqrt(RMSE_LSL / data_len_lsl)
    RMSE_N = sqrt(RMSE_N / data_len_n)
    RMSE_K = sqrt(RMSE_K / data_len_k)
end