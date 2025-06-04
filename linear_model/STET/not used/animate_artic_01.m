function []=animate_artic_01(q,nstep,PARSu,TRACK)
% derived from ART1c, but with track removed
% use for open-loop testing
% similar to animate_car1.m etc.
% q=generalized coords: 
% q1=x of front point (front joint)
% q2=y
% q3=yaw angle of unit 1
% q4=yaw angle of unit 2
% q5=yaw angle of unit 3
% Xa and Xj are added points, format Xa(2 x Na x N) etc.

if exist('nstep','var')==0, nstep=1;end
if exist('TRACK','var')==0, TRACK=[];end
if ~isempty(TRACK)
    [~,marking]=map(TRACK,[-2,2]);
    close
end
small=1e-3;
%find the number of active joints and axles
% test=squeeze(sum(abs(Xj(1,:,:)),3)); Nj=find(test>small,1,'last');
% test=squeeze(sum(abs(Xa(1,:,:)),3)); Na=find(test>small,1,'last');

figure('DoubleBuffer','on');
% axis('equal')
set(gca,'XTickMode','manual','XTick',[-1000:10:1000])
set(gca,'YTickMode','manual','YTick',[-1000:10:1000])
grid on
frame=[-10 10 -10 10]*4;
cam=[q(:,1),q(:,2)];
cam=cam-10*[cos(q(:,3)),sin(q(:,3))];
cam=cam-5*[cos(q(:,4)),sin(q(:,4))]; %so approx at centre of middle unit

for i=1:nstep:length(q)
    frameshift=[cam(i,1),cam(i,1),cam(i,2),cam(i,2)];
    cla
    artic_shape_01(q(i,:),PARSu)
    if ~isempty(TRACK)
        plot(marking(:,1),marking(:,2),'g-');
        plot(marking(:,3),marking(:,4),'g-');
    end
    axis('equal')
    grid on
    
    axis(frame+frameshift)
    drawnow
    if i==1
        pause
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



