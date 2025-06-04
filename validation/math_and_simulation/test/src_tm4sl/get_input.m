function [tstr_in, str_in] = get_input(k)

%% explain
% k defines the type of test
%   1 = step-steer
%   2 = ramp (slowly increasing steer)
%   3 = random-steer
%   4 = sinus
%   5 = square
% swa_in at the steering wheel, in radians

% k = 5;

deg = pi/180;
T = 100; 
Ts = 10;

%% choose situation
switch k
    % step
    case 1 
        amp = 20 * deg; 
        dt = 0.2;
        pts = [0 0
               Ts 0
               Ts + dt,amp
               T, amp];
    % ramp
    case 2  
        amp = 180 * deg;
        tf = 50;
        pts = [0 0
               Ts 0
               tf amp
               T amp];
    % random steer    
    case 3 
        dt = 0.1;
        tf = Ts + 30;
        stdev = 10 * deg;
        n = round((tf - Ts) / dt);
        tvec = (1:n)' * dt;
        strvec = stdev * randn(size(tvec));
        pts = [0 0
               Ts 0
               Ts + tvec, strvec
               tf + dt 0
               T 0];
    % series of sinusoids   
    case 4 
        dt = 0.01;
        % amplitudes and frequencies set below 
        % (T = 2s, A = 5deg, length = 40s, 即20个T)
        f = 0.5;
        tvec1 = (0:dt:20/f)';
        strvec1 = 5*deg*sin(2*pi*f*tvec1);
        pts = [0 0
               Ts 0
               Ts + dt + tvec1, strvec1];
        tf = pts(end, 1);
        pts = [pts
               tf + dt, 0
               T 0];
    % square wave  
    case 5 
        amp = 45*deg; 
        dt = 5;
        ddt = 0.2;
        pts = [0 0
               Ts 0
               Ts + ddt, amp
               Ts + dt, amp
               Ts + dt + ddt, -amp
               Ts + 2*dt, -amp
               Ts + 2*dt + ddt,amp
               Ts + 3*dt, amp
               Ts + 3*dt + ddt, -amp
               Ts + 4*dt, -amp
               Ts + 4*dt + ddt, amp
               Ts + 5*dt, amp
               Ts + 5*dt + ddt, -amp
               Ts + 6*dt, -amp
               Ts + 6*dt + ddt, 0
               T, 0];
end

t_in = pts(:, 1)';
swa_in = pts(:, 2)';
            
%% test/plot the function

% plot(t_in, swa_in*180/pi);

tstr_in = (0:0.01:100)';
str_in = interp1(t_in, swa_in, tstr_in, 'linear', 'extrap');
% plot(tstr_in, str_in*180/pi);

end
