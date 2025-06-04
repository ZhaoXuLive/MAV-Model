% clear;clc;
% close all;
%% 车身相关参数
l = [8.9 7.4 1.5 7.9 1.5 7.4 1.5 0];
ld = diag(l,0);
l_1 = 11.4;
l_2 = 9.4; %铰接点间距
l_3 = 8.9; %后铰接点到第8轴距离

sampleTime = 0.1;
%% 读取数据

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('lagrange.mat');
% load('matlab0814.mat');
% load('circle50.mat');
% q = out.q_out;
% qd = out.qd_out;
% TRACK = pts2trk_v2(roadRecord);

q = out.q_out;
qd = out.qd_out;
qdd = out.qdd_out;

latAcc1 = qdd(:,2) - 2.5*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J1 = qdd(:,2) - 11.4*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J2 = latAcc_J1 - 9.4*(qdd(:,4).*cos(q(:,4)) - qd(:,4).*sin(q(:,4)));
latAcc3 = latAcc_J2 - 8.9*(qdd(:,5).*cos(q(:,5)) - qd(:,5).*sin(q(:,5)));

time = out.tout;

load('track_circle50.mat');
TRACK = temptrack;
% TRACK = trackgen4(2,50);

yaw1 = q(:,3);
yaw2 = q(:,4);
yaw3 = q(:,5);

%% 铰接点坐标计算
x_J1 = q(:,1) - l_1.*cos(yaw1);
y_J1 = q(:,2) - l_1.*sin(yaw1);

x_J2 = x_J1 - l_2.*cos(yaw2);
y_J2 = y_J1 - l_2.*sin(yaw2);

%% 8轴坐标 T车横摆角 横摆角速度计算
x_A8 = x_J2 - l_3.*cos(yaw3);
y_A8 = y_J2 - l_3.*sin(yaw3);

%% 各轴坐标计算
x0 = [x_J1 x_J1 x_J1 x_J2 x_J2 x_A8 x_A8 x_A8];
y0 = [y_J1 y_J1 y_J1 y_J2 y_J2 y_A8 y_A8 y_A8];
yaw = [yaw1 yaw1 yaw1 yaw2 yaw2 yaw3 yaw3 yaw3 ];
x = x0 + cos(yaw)*ld;
y = y0 + sin(yaw)*ld;

s_ref = TRACK(:,1);
x_ref = TRACK(:,2);
y_ref = TRACK(:,3);

%% 对于各段轨迹分别投影 计算offtrack

s = zeros(size(x));
offtrack = zeros(size(y));
s_J1 = zeros(size(x_J1));
s_J2 = zeros(size(x_J2));
offtrack_J1 = zeros(size(y_J1));
offtrack_J2 = zeros(size(y_J2));

[s_J1,offtrack_J1] = segment_project(x_ref,y_ref,s_ref,x_J1,y_J1);
[s_J2,offtrack_J2] = segment_project(x_ref,y_ref,s_ref,x_J2,y_J2);
for j = 1:8
    [s(:,j),offtrack(:,j)] = segment_project(x_ref,y_ref,s_ref,x(:,j),y(:,j));
end  

%% 画图
close all

colors = ["#EDB120","#4DBEEE","#77AC30","#0072BD","#7E2F8E","#6f6f6f","#bb9727","#c76da2","#D95319","#A2142F","#9ac9db"];

% figure(1);
% hold on;
% title('Trajctory of 8 Axle， Dynamic');
% xlabel('x (East) [m]');
% ylabel('y (North) [m]');
% axis equal;
% xlim([-50 400]);
% ylim([-20 120]);
% 
% plot(x_ref,y_ref,'color',colors(11),'Linewidth', 2);
% v = hypot(qd(:,1),qd(:,2));
% for i = 1:floor(size(q,1)/10)-1
%     plot(x_J1(1+(i-1)*10:(i-1)*10+10),y_J1(1+(i-1)*10:(i-1)*10+10),'color',colors(9),'Linewidth', 2);
%     plot(x_J2(1+(i-1)*10:(i-1)*10+10),y_J2(1+(i-1)*10:(i-1)*10+10),'color',colors(10),'Linewidth', 2);
% %     plot(x(1+(i-1)*10:(i-1)*10+10,:),y(1+(i-1)*10:(i-1)*10+10,:),'b','Linewidth', 2);
%     for i = 1:8
%         plot(x(1+(i-1)*10:(i-1)*10+10,i),y(1+(i-1)*10:(i-1)*10+10,i),'color',colors(i),'Linewidth', 2);
%     end
% %     pause(0.001/v(1+(i-1)*10));
% %     pause(0.001);
% end

figure(1);
hold on;
title('Trajctory of 2 Joint and Axle 1&8');
xlabel('x (East) [m]');
ylabel('y (North) [m]');
axis equal;
xlim([0 400]);
ylim([-10 40]);
% xlim([0 350]);
% ylim([-10 120]);
% xlim([0 120]);
% ylim([-10 10]);
plot(x_ref,y_ref,'color',colors(11),'Linewidth', 2);
plot(x_J1,y_J1,'color',colors(9),'Linewidth', 2);
plot(x_J2,y_J2,'color',colors(10),'Linewidth', 2);
plot(x(:,1),y(:,1),'color',colors(1),'Linewidth', 2);
plot(x(:,8),y(:,8),'color',colors(8),'Linewidth', 2);
legend('Ref','Joint 1','Joint 2','Axle 1','Axle 8');%'Joint1','Joint2',
hold off;

figure(2);
hold on;
title('Trajctory of 8 Axle');
xlabel('x (East) [m]');
ylabel('y (North) [m]');
axis equal;
xlim([0 400]);
ylim([-10 40]);
% xlim([0 250]);
% ylim([-150 10]);
% xlim([0 350]);
% ylim([-10 120]);
% xlim([0 120]);
% ylim([-10 10]);
plot(x_ref,y_ref,'color',colors(11),'Linewidth', 2);
for i = 1:8
    if i~=2&&i~=7
        plot(x(:,i),y(:,i),'color',colors(i),'Linewidth', 2);
    end
end
legend('Ref','Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8');%'Joint1','Joint2',
hold off;

figure(3);
% plot(s,offtrack,'Linewidth', 2);
title('Offtrack of 8 Axle');
xlabel('s [m]');
ylabel('offtrack  [m]');
% xlim([0 120]);
% ylim([-0.05 0.05]);
hold on;
for i = 1:8
    if i~=2&&i~=7
        plot(s(:,i),offtrack(:,i),'color',colors(i),'Linewidth', 2);
    end
end
legend('Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8');

figure(4);
% plot(s,offtrack,'Linewidth', 2);
title('Offtrack of 2 Joint and Axle 1&8');
xlabel('s [m]');
ylabel('offtrack  [m]');
% xlim([0 120]);
% ylim([-0.05 0.05]);
hold on;
plot(s_J1,offtrack_J1,'color',colors(10),'Linewidth', 2);
plot(s_J2,offtrack_J2,'color',colors(11),'Linewidth', 2);
plot(s(:,1),offtrack(:,1),'color',colors(1),'Linewidth', 2);
plot(s(:,8),offtrack(:,8),'color',colors(8),'Linewidth', 2);

legend('Joint 1','Joint 2','Axle 1','Axle 8');

figure(5);
hold on;
RWA_yawrate = max(abs(qd(:,5)))/max(abs(qd(:,3)));

title(['YawRate of 3 Unit, RWA = ',num2str(RWA_yawrate,6) ]);
xlabel('Time [s]');
ylabel('YawRate  [rad/s]');

plot(time,qd(:,3),'color',colors(1),'Linewidth', 2);
plot(time,qd(:,4),'color',colors(2),'Linewidth', 2);
plot(time,qd(:,5),'color',colors(3),'Linewidth', 2);

legend('YawRate 1','YawRate 2','YawRate 3');

figure(6);
hold on;
RWA_latAcc = max(abs(latAcc3))/max(abs(latAcc1));

title(['Lateral Acceleration of 1&3 Unit, RWA = ',num2str(RWA_latAcc,6) ]);
xlabel('Time [s]');
ylabel('LatAcc  [m/(s^2)]');

plot(time,latAcc1,'color',colors(1),'Linewidth', 2);
plot(time,latAcc3,'color',colors(3),'Linewidth', 2);

plot(time,latAcc_J1,'color',colors(4),'Linewidth', 2);
plot(time,latAcc_J2,'color',colors(5),'Linewidth', 2);

legend('LatAcc 1','LatAcc 3','LatAcc J1','LatAcc J2');

%% 投影函数 整段
function [s,offtrack] = segment_project(x_ref,y_ref,s_traj,x,y)
    s = zeros(size(x));
    offtrack = zeros(size(y));
%     index = 0;
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
%         if(index ~= 0)
%             index_start = index - 150;
%             index_end = index + 150;
%         end
        if(index_start<1)
            index_start = 1;
        end
        if(index_end>size(s_traj,1))
            index_end = size(s_traj,1);
        end

%         [s(i),offtrack(i)] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));

        [s(i),offtrack(i),index] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));
    end
end

%% 投影函数 单点

function [s,d,index,err_all] = cart2frenet(x_ref,y_ref,s_traj,x_point,y_point)
% function [s,d] = cart2frenet(traj,section,sectionLength,x_point,y_point)
% x_ref = x_ref(index_start:index_end);
% y_ref = y_ref(index_start:index_end);
traj = [x_ref y_ref];
% point = [x_point y_point];
section = diff(traj); % 各路径段向量
vector = [x_point - x_ref y_point - y_ref];% 被投影点到各路径点的距离
vector = vector(1:end-1,:);
sectionLength = sum(section.*section,2);% 各路径段长度平方
t = (vector(:,1).*section(:,1)+vector(:,2).*section(:,2)) ./ sectionLength;% 被投影点在各路径的投影与相应路径段长度的比值
point_proj = zeros(size(section));

% t_temp = t(end);
% t(t < 0.0) = 0;
% t(t > 1.0) = ;

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
% s = s_traj(index_start+index) + vector(index,1).*section(index,1)+vector(index,2).*section(index,2);
s = s_traj(index) + (vector(index,1).*section(index,1)+vector(index,2).*section(index,2))/sqrt(sectionLength(index)); %%%%%%
d = hypot(x_point - point_proj_out(1), y_point - point_proj_out(2));
dir = (x_point - point_proj_out(1))*section(index,2) - (y_point - point_proj_out(2))*section(index,1);
d = d*sign(dir);

% vector = point - section;

end
