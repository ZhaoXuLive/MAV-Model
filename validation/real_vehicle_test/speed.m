data_mode = 4;
deg = pi / 180;
dt = 0.01;
if data_mode == 1
    load('quarter_turn.mat')
    load('quarter_turn_steer.mat')
    startpoint = 1;
    endpoint = 10000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 2
    load('circle.mat')
    load('circle_steer.mat')
    startpoint = 1;
    endpoint = 38000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 3
    load('double_change.mat')
    load('double_change_steer.mat')
    startpoint = 1;
    endpoint = 7000;
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 4
    load('2318.mat')
    startpoint = 1;
    endpoint = 17000;  % 20 - 100
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 5
    load('2122.mat')
    startpoint = 1;
    endpoint = 8500;  % 6000-8000
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 6
    load('2146.mat')
    startpoint = 1;
    endpoint = 13000; % 80-100
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end
if data_mode == 7
    load('2148.mat')
    startpoint = 1;
    endpoint = 7500;  % 10-20 55-70
    q_debug = qqd_est_out(startpoint:endpoint, 1:5);
    qd_debug = qqd_est_out(startpoint:endpoint, 6:10);
    steer = steer_m(startpoint:endpoint,:) * deg;
    data_len = size(qqd_est_out(startpoint:endpoint,:),1);
    time = linspace(0,(data_len-1)*0.01,data_len)';
end

speed_a = zeros(data_len, 1);

for i = 1:data_len
    speed_a(i) = hypot(qd_debug(i, 1), qd_debug(i, 2));
end

plot(time, speed_a);
hold on
plot(time, qd_debug(:, 1));
hold on
plot(time, qd_debug(:, 2));
legend("all", "x", "y")
