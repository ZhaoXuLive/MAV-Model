
%% 车身相关参数
l = [8.9 7.4 1.5 7.9 1.5 7.4 1.5 0];
ld = diag(l,0);
l_1 = 11.4;
l_2 = 9.4; %铰接点间距
l_3 = 8.9; %后铰接点到第8轴距离

sampleTime = 0.1;
%% 读取数据

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = out.q_out;
qd = out.qd_out;
% TRACK = trackgen4(2,50);
% TRACK = pts2trk_v2(roadRecord);

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

figure(1);
hold on;
xlim([-50 350]);
ylim([-20 120]);
plot(x_ref,y_ref,'Linewidth', 2);
% v = hypot(qd(:,1),qd(:,2));
for i = 1:floor(size(q,1)/10)-1
    plot(x_J1(1+(i-1)*10:(i-1)*10+10),y_J1(1+(i-1)*10:(i-1)*10+10),'r','Linewidth', 2);
%     pause(0.001/v(1+(i-1)*10));
    pause(0.001);
end

figure(2);
hold on;
plot(x_J1,y_J1);
plot(x_J2,y_J2);
for i = 1:8
%     for j = 1:size()
    plot(x(:,i),y(:,i));
end
hold off;

figure(3);
plot(s,offtrack);

%% 投影函数 整段
function [s, offtrack] = segment_project(x_ref,y_ref,s_traj,x,y)
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
%         if(index ~= 0)
%             index_start = index - 150;
%             index_end = index + 150;
%         end
%         if(index_start<1)
%             index_start = 1;
%         end
%         if(index_end>size(s_traj,1))
%             index_end = size(s_traj,1);
%         end

%         [s(i),offtrack(i)] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));

        [s(i),offtrack(i), ~] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));
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
    if ((t(i) > 0.0 && t(i) < 1.0) || (t(i) >= 1.0 && i == size(t,1)) || (t(i) <= 0.0 && i <= 20))
        point_proj(i,:) = traj(i,:) + t(i)*section(i,:);
    elseif(t(i) <= 0.0)
        point_proj(i,:) = traj(i,:);
    else%(t(i) >= 1.0)
        point_proj(i,:) = traj(i+1,:);
    end
end

err_all = hypot(point_proj(:,1) - x_point, point_proj(:,2) - y_point);
[~, index] = min(err_all);
point_proj_out = [point_proj(index,1),point_proj(index,2)];
% s = s_traj(index_start+index) + vector(index,1).*section(index,1)+vector(index,2).*section(index,2);
s = s_traj(index) + vector(index,1).*section(index,1)+vector(index,2).*section(index,2);
d = hypot(x_point - point_proj_out(1), y_point - point_proj_out(2));
dir = (x_point - point_proj_out(1))*section(index,2) - (y_point - point_proj_out(2))*section(index,1);
d = d*sign(dir);

% vector = point - section;

end
