%% date load

close all

global PARSu
global PARSa

Q_ART_par

% time_tm = 0:0.001:100;
% q_out_tm = q_out;
% qd_out_tm = qd_out;
% STRa_out_tm = steer8_out;
% Vgg_out_tm = Vgg_out;

sty={'r:','r-','r--','b:','b-','b--'};
lines = ["r-"; ":"; "-."; "--"; "--."; ":."];
axle_str={'Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8'};
unit_str={'Unit 1 (math model)','Unit 2 (math model)','Unit 3 (math model)', 'Unit 1 (truckmaker)', 'Unit 2 (truckmaker)', 'Unit 3 (truckmaker)'};

%% PLOTTING

[Agy_tm, Rg_tm] = calculate_data(time_tm, q_out_tm, qd_out_tm);
[Agy, Rg] = calculate_data(time, q_out, qd_out);

% case1 : 2091 20000
% case4 : 3800 20000

% yaw rates
% figure('Name','yaw rate')
subplot(1,2,1)
for i=1:Nu
    plot(time(1:2091),qd_out(1:2091,i+2)*180/pi,sty{i},'LineWidth', 1),hold on
end
for i=1:Nu
    plot(time_tm(1:20000),qd_out_tm(1:20000,i+2)*180/pi,sty{i+3},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontName','Times New Roman','FontSize',12);
ylabel('Yaw Rate [deg/s]','FontName','Times New Roman','FontSize',12);
legend(unit_str,'FontName','Times New Roman','FontSize',9);
title('Yaw Rate of Three Unit','FontName','Times New Roman','FontSize',12);
axis([0 20 -1 3.5]);
% axis([0 20 -0.3 0.5])

% lateral accelerations
subplot(1,2,2)
for i=1:Nu
    plot(time(1:2091),Agy(1:2091,i),sty{i},'LineWidth', 1),hold on
end
for i=1:Nu
    plot(time_tm(1:20000),Agy_tm(1:20000,i),sty{i+3},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]','FontName','Times New Roman','FontSize',12);
ylabel('Lateral Acceleration [ms^{-2}]','FontName','Times New Roman','FontSize',12);
legend(unit_str,'FontName','Times New Roman','FontSize',9);
title('Lateral Acceleration of Three Unit','FontName','Times New Roman','FontSize',12);
axis([0 20 -0.1 0.3]);
% axis([0 20 -0.06 0.12])

%% data calculate
function [Agy, Rg] = calculate_data(time, q_out, qd_out)

%% calculate

global PARSu
global PARSa
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
    
end