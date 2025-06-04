function [t_in,swa_in]=open_loop_steer(k)
% open loop steering inputs for 'repeating sequence' block in simulink
% k defines the type of test
%   1 = step-steer
%   2 = ramp (slowly increasing steer)
%   3 = random-steer
%   4 = sinus
%
% swa_in at the steering wheel, in radians

% further adjustable parameters below

deg=pi/180;
T=100; %default repeat time

switch k
    
    case 1 %step
        t0=5; amp=20*deg; dt=0.2;
        pts=[
            0 0
            t0 0
            t0+dt,amp
            T, amp];
        
    case 2 %ramp
        t0=2; amp=180*deg;tf=50;
        pts=[
            0 0
            t0 0
            tf amp
            T amp];
        
    case 3 %random steer
        t0=2;dt=0.1;
        tf=t0+30;
        stdev=10*deg;
        n=round((tf-t0)/dt);
        tvec=(1:n)'*dt;
        strvec=stdev*randn(size(tvec));
        pts=[
            0 0
            t0 0
            t0+tvec,strvec
            tf+dt 0
            T 0];
        
    case 4 %series of sinusoids (can comment out those not needed)
        t0=4;dt=0.01;
        %amplitudes and frequencies set below
        f=0.5;tvec1=(0:dt:20/f)';strvec1=5*deg*sin(2*pi*f*tvec1);
        pts=[
            0 0
            t0 0
            t0+dt+tvec1, strvec1];
        tf=pts(end,1);
        pts=[pts
            tf+dt,0
            T 0];
    case 5 %series of sinusoids (can comment out those not needed)
        t0=4;dt=0.01;
        %amplitudes and frequencies set below
        f=1.0;tvec2=(0:dt:40/f)';strvec2=5*deg*sin(2*pi*f*tvec2);
        pts=[
            0 0
            t0 0
            t0+dt+tvec2, strvec2];
        tf=pts(end,1);
        pts=[pts
            tf+dt,0
            T 0];
    case 6 %series of sinusoids (can comment out those not needed)
        t0=4;dt=0.01;
        %amplitudes and frequencies set below
        ts=8; %settling time between cycles
        f=0.5;tvec1=(0:dt:20/f)';strvec1=10*deg*sin(2*pi*f*tvec1);
        pts=[
            0 0
            t0 0
            t0+dt+tvec1, strvec1];
        tf=pts(end,1);
        pts=[pts
            tf+dt,0
            T 0];
     case 7 %series of sinusoids (can comment out those not needed)
        t0=4;dt=0.01;
        %amplitudes and frequencies set below
        f=2;tvec2=(0:dt:80/f)';strvec2=5*deg*sin(2*pi*f*tvec2);
        pts=[
            0 0
            t0 0
            t0+dt+tvec2, strvec2];
        tf=pts(end,1);
        pts=[pts
            tf+dt,0
            T 0];
    
        
     case 8 %lane change sinus
        t0=2;dt=0.01;
        %amplitudes and frequencies set below
        
        f=0.3;tvec1=(0:dt:1/f)';strvec1=70*deg*sin(2*pi*f*tvec1);

        pts=[
            0 0
            t0 0
            t0+dt+tvec1, strvec1];
        tf=pts(end,1);
        pts=[pts
            tf+dt,0
            T 0];
        
      case 9 %square wave
        t0=10; amp=45*deg; dt=5;ddt=0.2;
        pts=[
            0 0
            t0 0
            t0+ddt,amp
            t0+dt, amp
            t0+dt+ddt,-amp
            t0+2*dt,-amp
            t0+2*dt+ddt,amp
            t0+3*dt, amp
            t0+3*dt+ddt,-amp
            t0+4*dt,-amp
            t0+4*dt+ddt,amp
            t0+5*dt, amp
            t0+5*dt+ddt,-amp
            t0+6*dt,-amp
            t0+6*dt+ddt,0
            T,0];
end

t_in=pts(:,1)';
swa_in=pts(:,2)';


            
return

%% test/plot the function
close all
[t_in,swa_in]=open_loop_steer(1);
figure, plot(t_in,swa_in*180/pi)
            
            
            
            
