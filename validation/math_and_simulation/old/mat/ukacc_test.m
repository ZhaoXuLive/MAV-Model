%%
clear; clc; close all
 
Q_ART_par

% 设定图片中的统一字体和字号。单幅图片可以另行设定。设定20号字主要考虑图片缩小后放在论文中的效果。
set(0,'DefaultAxesFontname','Times New Roman');
set(0,'DefaultAxesFontsize',20);

% 建立一个专用的文件夹，用于存放数据和图片。可根据实际情况设定。
FigurePath = ['./SimuData/figures/', datestr(now,'yyyymmdd')];
mkdir(FigurePath)

%% date load

global PARSu
global PARSa

mode = 1;
if mode == 1
    load('truckmaker_case1.mat')
    time_tm = 0:0.001:100;
    q_out_tm = q_out;
    qd_out_tm = qd_out;
    STRa_out_tm = steer8_out;
    Vgg_out_tm = Vgg_out;
    load('lagrange_step.mat')
else
    load('truckmaker_case4.mat')
    time_tm = 0:0.001:100;
    q_out_tm = q_out;
    qd_out_tm = qd_out;
    STRa_out_tm = steer8_out;
    Vgg_out_tm = Vgg_out;
    load('lagrange_sin.mat')
end

sty={'r-','g-','b-','r--','g--','b--'};
axle_str={'Axle 1','Axle 2','Axle 3','Axle 4','Axle 5','Axle 6','Axle 7','Axle 8'};
unit_str={'Unit 1 (Math)','Unit 2 (Math)','Unit 3 (Math)', 'Unit 1 (IPG)', 'Unit 2 (IPG)', 'Unit 3 (IPG)'};

%% PLOTTING

[Agy_tm, Rg_tm] = calculate_data(time_tm, q_out_tm, qd_out_tm);
[Agy, Rg] = calculate_data(time, q_out, qd_out);

% yaw rates
h1 = figure;
for i=1:Nu
    plot(time(1:2091),qd_out(1:2091,i+2)*180/pi,sty{i},'LineWidth', 1),hold on
end
for i=1:Nu
    plot(time_tm(1:20000),qd_out_tm(1:20000,i+2)*180/pi,sty{i+3},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]');
ylabel('Yaw rate [deg/s]');
set(gcf, 'Position', [100, 100, 800, 300]); % 设置Figure窗口的位置和大小
legend(unit_str, 'Location','EastOutside');
if mode == 1
    axis([0 20 -0.5 2]);
    print(h1,'-depsc2','-loose',[FigurePath,'yawrate_ms_step.eps']); 
    savefig(h1, [FigurePath,'/yawrate_ms.fig']);
else
    axis([0 11 -0.3 0.3]);
    print(h1,'-depsc2','-loose',[FigurePath,'yawrate_ms_sin.eps']); 
    savefig(h1, [FigurePath,'/yawrate_ms.fig']);
end

h2 = figure;
for i=1:Nu
    plot(time(1:2091),Agy(1:2091,i),sty{i},'LineWidth', 1),hold on
end
for i=1:Nu
    plot(time_tm(1:20000),Agy_tm(1:20000,i),sty{i+3},'LineWidth', 1),hold on
end
grid on
xlabel('Time [s]');
ylabel('Lateral acc. [m/s^{2}]');
set(gcf, 'Position', [100, 100, 800, 300]); % 设置Figure窗口的位置和大小
legend(unit_str, 'Location','EastOutside');
if mode == 1
    axis([0 20 -0.1 0.2]);
    print(h2,'-depsc2','-loose',[FigurePath,'latAcc_ms_step.eps']); 
    savefig(h2, [FigurePath,'/yawrate_ms.fig']);
else
    axis([0 11 -0.08 0.08]);
    print(h2,'-depsc2','-loose',[FigurePath,'latAcc_ms_sin.eps']); 
    savefig(h2, [FigurePath,'/yawrate_ms.fig']);
end

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