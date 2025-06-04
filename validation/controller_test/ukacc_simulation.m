%%
clear; clc; close all
 
call_AFG_s1v1

% 设定图片中的统一字体和字号。单幅图片可以另行设定。设定20号字主要考虑图片缩小后放在论文中的效果。
set(0,'DefaultAxesFontname','Times New Roman');
set(0,'DefaultAxesFontsize',20);

% 建立一个专用的文件夹，用于存放数据和图片。可根据实际情况设定。
FigurePath = ['./SimuData/figures/', datestr(now,'yyyymmdd')];
mkdir(FigurePath)

%%
sampleTime = 0.1;
sim_model = 'model2_19b'; 
sim(sim_model, 200)

%%
q = ans.q_out;
qd = ans.qd_out;
qdd = ans.qdd_out;

latAcc1 = qdd(:,2) - 2.5*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J1 = qdd(:,2) - 11.4*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J2 = latAcc_J1 - 9.4*(qdd(:,4).*cos(q(:,4)) - qd(:,4).*sin(q(:,4)));
latAcc3 = latAcc_J2 - 8.9*(qdd(:,5).*cos(q(:,5)) - qd(:,5).*sin(q(:,5)));

time = ans.tout;

load('track_circle50.mat');
TRACK = temptrack;

l = [8.9 7.4 1.5 7.9 1.5 7.4 1.5 0];
ld = diag(l,0);
l_1 = 11.4;
l_2 = 9.4; 
l_3 = 8.9; 

yaw1 = q(:,3);
yaw2 = q(:,4);
yaw3 = q(:,5);

x_J1 = q(:,1) - l_1.*cos(yaw1);
y_J1 = q(:,2) - l_1.*sin(yaw1);

x_J2 = x_J1 - l_2.*cos(yaw2);
y_J2 = y_J1 - l_2.*sin(yaw2);

x_A8 = x_J2 - l_3.*cos(yaw3);
y_A8 = y_J2 - l_3.*sin(yaw3);

x0 = [x_J1 x_J1 x_J1 x_J2 x_J2 x_A8 x_A8 x_A8];
y0 = [y_J1 y_J1 y_J1 y_J2 y_J2 y_A8 y_A8 y_A8];
yaw = [yaw1 yaw1 yaw1 yaw2 yaw2 yaw3 yaw3 yaw3 ];
x = x0 + cos(yaw)*ld;
y = y0 + sin(yaw)*ld;

s_ref = TRACK(:,1);
x_ref = TRACK(:,2);
y_ref = TRACK(:,3);

s = zeros(size(x));
offtrack = zeros(size(y));

[s_J1,offtrack_J1] = segment_project(x_ref,y_ref,s_ref,x_J1,y_J1);
[s_J2,offtrack_J2] = segment_project(x_ref,y_ref,s_ref,x_J2,y_J2);
for j = 1:8
    [s(:,j),offtrack(:,j)] = segment_project(x_ref,y_ref,s_ref,x(:,j),y(:,j));
end  

%%

colors = ["#EDB120","#4DBEEE","#77AC30","#0072BD","#7E2F8E","#6f6f6f","#bb9727","#c76da2","#D95319","#A2142F","#9ac9db"];
lines = ["-"; ":"; "-."; "--"; "--."; ":."];

h1 = figure;
% set(get(gca, 'Xlabel'));
% set(get(gca, 'Ylabel'));
plot(x_ref,y_ref,'b','Linewidth', 1);
hold on
j = 1;
for i = 1:8
    if i~=2&&i~=7
        plot(x_ref,y_ref,lines(j),'Linewidth', 1);
        hold on
        j = j+1;
    end
end
xlabel('X (East) [m]');
ylabel('Y (North) [m]');
legend('Ref','Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8');
axis([0 450 -10 200])
print(h1,'-depsc2','-loose',[FigurePath,'traj.eps']); 
savefig(h1, [FigurePath,'/traj.fig']);

h2 = figure;
% set(get(gca, 'Xlabel'));
% set(get(gca, 'Ylabel'));
j = 1;
for i = 1:8
    if i~=2&&i~=7
        plot(s(:,i),offtrack(:,i),lines(j),'Linewidth', 1);
        hold on
        j = j+1;
    end
end
xlabel('Travel distance [m]');
ylabel('Offtrack [m]');
legend('Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8');
print(h2,'-depsc2','-loose',[FigurePath,'offtrack.eps']); 
savefig(h2, [FigurePath,'/offtrack.fig']);

h3 = figure;
RWA_yawrate = max(abs(qd(:,5)))/max(abs(qd(:,3)));
% set(get(gca, 'Xlabel'));
% set(get(gca, 'Ylabel'));
plot(time,qd(:,3),lines(1),'Linewidth', 1);
hold on
plot(time,qd(:,4),lines(2),'Linewidth', 1);
hold on
plot(time,qd(:,5),lines(3),'Linewidth', 1);
xlabel('Time [s]');
ylabel('Yaw rate  [rad/s]');
legend('Unit 1','Unit 2','Unit 3');
print(h3,'-depsc2','-loose',[FigurePath,'RWA_yawrate.eps']); 
savefig(h3, [FigurePath,'/RWA_yawrate.fig']);

h4 = figure;
RWA_latAcc = max(abs(latAcc3))/max(abs(latAcc1));
% set(get(gca, 'Xlabel'));
% set(get(gca, 'Ylabel'));
title(['Lateral acc. of 1&3 Unit, RWA = ',num2str(RWA_latAcc,6) ]);
plot(time,latAcc1,lines(1),'Linewidth', 1);
hold on
plot(time,latAcc3,lines(2),'Linewidth', 1);
hold on
plot(time,latAcc_J1,lines(3),'Linewidth', 1);
hold on
plot(time,latAcc_J2,lines(4),'Linewidth', 1);
hold on
xlabel('Time [s]');
ylabel('Lateral acc. [m/s^2]');
legend('First Unit','Third Unit','Joint 1','Joint 2');
print(h4,'-depsc2','-loose',[FigurePath,'RWA_latAcc.eps']); 
savefig(h4, [FigurePath,'/RWA_latAcc.fig']);

function [s,offtrack] = segment_project(x_ref,y_ref,s_traj,x,y)
    s = zeros(size(x));
    offtrack = zeros(size(y));
    for i = 1:size(x,1)
        index_start = 1;
        index_end = size(s_traj,1);
        if(i~= 1)
           [~, index_start] = min(abs(s(i-1)-50 - s_traj));
           [~, index_end] = min(abs(s(i-1)+50 - s_traj));
        end
        
        if(index_end-index_start < 10)
            index_end = index_end+10;
            index_start = index_start-10;
        end
        if(index_end-index_start == 1)
            disp('fxxxk');
        end
        if(index_start<1)
            index_start = 1;
        end
        if(index_end>size(s_traj,1))
            index_end = size(s_traj,1);
        end
        [s(i),offtrack(i),index] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));
    end
end

function [s,d,index,err_all] = cart2frenet(x_ref,y_ref,s_traj,x_point,y_point)

    traj = [x_ref y_ref];
    section = diff(traj); 
    vector = [x_point - x_ref y_point - y_ref];
    vector = vector(1:end-1,:);
    sectionLength = sum(section.*section,2);
    t = (vector(:,1).*section(:,1)+vector(:,2).*section(:,2)) ./ sectionLength;
    point_proj = zeros(size(section));

    for i = 1:size(t,1)
        if ((t(i) > 0.0 && t(i) < 1.0) || (t(i) >= 1.0 && i == size(t,1)) || (t(i) <= 0.0 && i <= 1))
            point_proj(i,:) = traj(i,:) + t(i)*section(i,:);
        elseif(t(i) <= 0.0)
            point_proj(i,:) = traj(i,:);
        else%(t(i) >= 1.0)
            point_proj(i,:) = traj(i+1,:);
        end
    end
    err_all = hypot(point_proj(:,1) - x_point, point_proj(:,2) - y_point);
    [err,index] = min(err_all);
    point_proj_out = [point_proj(index,1),point_proj(index,2)];
    s = s_traj(index) + (vector(index,1).*section(index,1)+vector(index,2).*section(index,2))/sqrt(sectionLength(index)); 
    d = hypot(x_point - point_proj_out(1), y_point - point_proj_out(2));
    dir = (x_point - point_proj_out(1))*section(index,2) - (y_point - point_proj_out(2))*section(index,1);
    d = d*sign(dir);
end
