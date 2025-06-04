%% date load

close all

time = 0:0.001:100;
q_out = q_out;
qd_out = qd_out;
STRa_out = steer8_out;
Vgg_out = Vgg_out;
sty={'r.-','g-','b--','c:','m-.','b-','b--','b:','b-.','r-','r--','r:','r-.'};
axle_str={'axle 1','axle 2','axle 4','axle 5','axle 5','axle 6','axle 7','axle 8'};
unit_str={'unit 1','unit 2','unit 3'};

%% calculate

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

% calculate acc/yawratw 
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
    
%% PLOTTING

% final unit cg trajectories
figure
plot(squeeze(Rg(:,1,3)),squeeze(Rg(:,2,3)))
title('final unit cg trajectory')
xlabel('X(m)')
ylabel('Y(m)')

% axle steer angles
figure('Name','steer_angle')
for i=1:Na
    plot(time,STRa_out(:,i)*180/pi,sty{i},'LineWidth', 1);hold on
end
grid on
xlabel('Time [s]'),ylabel('Axle steer angle [deg]','FontSize',14)
legend(axle_str)
title('Steering Angle of Axles','FontSize',14)

% yaw rates
figure('Name','yawrate')
for i=1:Nu
    plot(time,qd_out(:,i+2)*180/pi,sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14),ylabel('Yaw rate [deg/s]','FontSize',14)
legend(unit_str)
title('Yaw Rate','FontSize',14);

% lateral accelerations
figure('Name','acc_lat')
for i=1:Nu
    plot(time,Agy(:,i),sty{i},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontSize',14),ylabel('CG lateral acceleration [ms^{-2}]','FontSize',14)
legend(unit_str)
title('Lateral Acceleration','FontSize',14);

% %% lat acc vs speed*yaw rate
% figure,nn=50;
% for i=1:Nu
%     plot(time,Agy(:,i)),hold on
%     plot(time(1:nn:end), speed(1:nn:end).*qd_out(1:nn:end,2+i),'r.')
%     if i<Nu, pause, cla, end
% end