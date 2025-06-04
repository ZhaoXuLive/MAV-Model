function [t_in,swa_in]=open_loop_steer(f, a)

% further adjustable parameters below

deg=pi/180;
t0=5; dt=0.01;
t1 = (0:dt:t0)';
str1 = zeros(t0/dt + 1, 1);

% series of sinusoids 

tvec=(0:dt:10/f)';
strvec=a*deg*sin(2*pi*f*tvec);

pts=[ t1, str1
     t0+dt+tvec, strvec];
% for i = 1:19
%     tf=pts(end,1);
%     pts=[pts
%         tf+dt,0
%         tf+dt+tvec, strvec];
% end
         
t_in=pts(:,1)';
swa_in=pts(:,2)';
         
return

%% test/plot the function

% [t_in,swa_in]=open_loop_steer(1);
% figure, plot(t_in,swa_in*180/pi)
            