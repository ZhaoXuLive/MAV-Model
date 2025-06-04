
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
qdd = out.qdd_out;

latAcc1 = qdd(:,2) - 2.5*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J1 = qdd(:,2) - 11.4*(qdd(:,3).*cos(q(:,3)) - qd(:,3).*sin(q(:,3)));
latAcc_J2 = latAcc_J1 - 9.4*(qdd(:,4).*cos(q(:,4)) - qd(:,4).*sin(q(:,4)));
latAcc3 = latAcc_J2 - 8.9*(qdd(:,5).*cos(q(:,5)) - qd(:,5).*sin(q(:,5)));

time = out.tout;

yaw1 = q(:,3);
yaw2 = q(:,4);
yaw3 = q(:,5);

%% 画图
close all

colors = ["#EDB120","#4DBEEE","#77AC30","#0072BD","#7E2F8E","#6f6f6f","#bb9727","#c76da2","#D95319","#A2142F","#9ac9db"];

figure(1);
hold on;
RWA_yawrate = max(abs(qd(:,5)))/max(abs(qd(:,3)));

title(['YawRate of 3 Unit, RWA = ',num2str(RWA_yawrate,6) ]);
xlabel('Time [s]');
ylabel('YawRate  [rad/s]');

plot(time,qd(:,3),'color',colors(1),'Linewidth', 2);
plot(time,qd(:,4),'color',colors(2),'Linewidth', 2);
plot(time,qd(:,5),'color',colors(3),'Linewidth', 2);

legend('YawRate 1','YawRate 2','YawRate 3');

figure(2);
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
