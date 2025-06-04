% basic post-processing

close all

time = ans.time;
q_out = ans.q_out;
qd_out = ans.qd_out;
STRa_out = ans.STRa_out;
Vgg_out = ans.Vgg_out;

%plotting

% plotting styles
sty={'r.-','g-','b--','c:','m-.','b-','b--','b:','b-.','r-','r--','r:','r-.'};
axle_str={'axle 1','axle 2','axle 4','axle 5','axle 5','axle 6','axle 7','axle 8'};
unit_str={'unit 1','unit 2','unit 3'};

% cell(Na,1); for i=1:Na, axle_str{i}=['axle ',num2str(i)];end
% unit_str={}; for i=1:Nu, unit_str{i}=['unit ',num2str(i)];end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~exist('Rg','var') 
%data postprocess
    
dt=time(2)-time(1);
Nt=length(time);
Nu=size(PARSu,1);
Na=size(PARSa,1);

L=PARSu(:,3);%unit lengths
Xg=PARSu(:,4); %CG location relative to rear joint

%joint and cg positions and velocities (globals)
Rj=zeros(Nt,2,Nu+1);Vj=zeros(Nt,2,Nu+1);
Rg=zeros(Nt,2,Nu);Vg=zeros(Nt,2,Nu);
Rj(:,:,1)=q_out(:,1:2);
Vj(:,:,1)=qd_out(:,1:2);
for i=1:Nu
    psi=q_out(:,2+i);psidot=qd_out(:,2+i);
    u=[cos(psi),sin(psi)];n=[-sin(psi),cos(psi)];
    r0=squeeze(Rj(:,:,i));v0=squeeze(Vj(:,:,i));%front joint
    r1=r0-L(i)*u;Rj(:,:,i+1)=r1;%rear joint
    rg=r1+Xg(i)*u;Rg(:,:,i)=rg;%cg location
    %velocities
    v1=v0-L(i)*(psidot*[1 1]).*n;Vj(:,:,i+1)=v1; %rear joint velocities
    vg=v1+Xg(i)*(psidot*[1 1]).*n;Vg(:,:,i)=vg; %cg velocities
end

%joint trajectories, with final point marked as circle for joint and dot
%for cg
figure
for i=1:Nu+1
    plot(squeeze(Rj(:,1,i)),squeeze(Rj(:,2,i))), hold on
end

% for i=1:Nu
%     plot(squeeze(Rj(end,1,i)),squeeze(Rj(end,2,i)),'ko')
%     plot(squeeze(Rg(end,1,i)),squeeze(Rg(end,2,i)),'r.')
% end

title('Articulated joint trajectory')
xlabel('X(m)')
ylabel('Y(m)')
legend('joint 1','joint 2','joint 3','joint 4')

%compare velocities with differentiated positions (as check)
if 0
    i=3;r0=squeeze(Rg(:,:,i));v0=squeeze(Vg(:,:,i)); %choose point to compare
    [~,dr]=zpempa2(r0,5);v=dr/dt;
    figure,nn=50;
    subplot(211), plot(time,v0(:,1)),hold on,plot(time(1:nn:end),v(1:nn:end,1),'r.')
    subplot(212), plot(time,v0(:,2)),hold on,plot(time(1:nn:end),v(1:nn:end,2),'r.')
end

Ag=zeros(Nt,2,Nu); %globals
Agx=zeros(Nt,Nu);Agy=zeros(Nt,Nu); %body coords
for i=1:Nu
    v0=squeeze(Vg(:,:,i)); %choose point to compare
    [~,dv]=zpempa2(v0,5);a=dv/dt;
    Ag(:,:,i)=a;
    psi=q_out(:,2+i);u=[cos(psi),sin(psi)];n=[-sin(psi),cos(psi)];
    Agx(:,i)=dot(a,u,2);
    Agy(:,i)=dot(a,n,2);
end

%speed from cg velocity of unit 1
vg1=squeeze(Vg(:,:,1));speed=hypot(vg1(:,1),vg1(:,2));

% Articulation angles
art_ang=q_out(:,3:Nu+1) - q_out(:,4:Nu+2);
art_ang_vel=qd_out(:,3:Nu+1)-qd_out(:,4:Nu+2);
    
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOTTING

% axle steer angles
figure('Name','steer_angle')
for i=1:Na
    plot(time,STRa_out(:,i)*180/pi,sty{i},'LineWidth', 1);hold on
end
grid on
xlabel('Time [s]'),ylabel('Axle steer angle [deg]','FontSize',14)
legend(axle_str)
title('Steering Angle of Axles','FontSize',14)
%%
% articulation angles and velocities
figure('Name','articulation');

subplot(2,1,1),
for i=1:Nu-1
    plot(time,art_ang(:,i)*180/pi,sty{i},'LineWidth', 1);hold on
end
ylabel('Angle[deg]','FontSize',14);
xlabel('Time [s]','FontSize',14)
title('Articulation Angle of Joints','FontSize',14);
grid on
hold on

subplot(2,1,2),
for i=1:Nu-1
    plot(time,art_ang_vel(:,i)*180/pi,sty{i},'LineWidth', 1);hold on
end
grid on
ylabel('Velocity[dps]','FontSize',14);
title('Angular Velocity of Joints','FontSize',14);
xlabel('Time [s]','FontSize',14)
legend({'joint 1','joint 2'})


%% speeds at each mass centre
figure('Name','vel_long')
for i=1:Nu
    y=squeeze(Vgg_out(1,i,:));
    plot(time,y/kph,sty{i}),hold on
end
grid on
xlabel('Time [s]','FontSize',14)
ylabel('Velocity [kph]','FontSize',14)
legend(unit_str)
title('Longitudinal Velocity Unit CG','FontSize',14);

%% lateral velocities at each mass centre
figure('Name','vel_lat')
for i=1:Nu
    y=squeeze(Vgg_out(2,i,:));
    plot(time,y/kph,sty{i}),hold on
end
grid on
xlabel('Time [s]','FontSize',14)
ylabel('Velocity [kph]','FontSize',14)
legend(unit_str)
title('Lateral Velocity Unit CG','FontSize',14);

%% side slip angle

figure('Name','vel_lat')
for i=1:Nu
    vx=squeeze(Vgg_out(1,i,:));
    vy=squeeze(Vgg_out(2,i,:));
    beta=-atan2(vy,vx);
    plot(time,beta*180/pi,sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14)
ylabel('Sideslip angle [deg]','FontSize',14)
legend(unit_str)
title('CG Sideslip Angle','FontSize',14)

%% yaw rates
% close all
figure('Name','yawrate')
for i=1:Nu
    plot(time,qd_out(:,i+2)*180/pi,sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14),ylabel('Yaw rate [deg/s]','FontSize',14)
legend(unit_str)
title('Yaw Rate','FontSize',14);

%% lateral accelerations
figure('Name','acc_lat')
for i=1:Nu
    plot(time,Agy(:,i),sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14),ylabel('CG lateral acceleration [ms^{-2}]','FontSize',14)
legend(unit_str)
title('Lateral Acceleration','FontSize',14);

%% longitudinal accelerations
figure('Name','acc_long')
for i=1:Nu
    plot(time,Agx(:,i),sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14),ylabel('CG longitudinal acceleration [ms^{-2}]','FontSize',14)
legend('unit 1','unit 2','unit 3','Location','Best')
title('Longitudinal Acceleration','FontSize',14);

return


%% lat acc vs speed*yaw rate
figure,nn=50;
for i=1:Nu
    plot(time,Agy(:,i)),hold on
    plot(time(1:nn:end), speed(1:nn:end).*qd_out(1:nn:end,2+i),'r.')
    if i<Nu, pause, cla, end
end