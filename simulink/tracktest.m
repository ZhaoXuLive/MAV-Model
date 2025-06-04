%% 参数
l_j1j2 = 11.4;
l_j2j3 = 9.4;
l_j1g1 = 11.4 - 5.538; 
l_j2g2 = 9.4 - 4.757;
l_j3g3 = 5.538; 
l_j1a1 = 2.5; 
l_j1a2 = 4;
l_j1a3 = 9.9;
l_j2a4 = 1.5;
l_j2a5 = 7.9;
l_j3a6 = 1.5;
l_j3a7 = 7.4;
l_j3a8 = 8.9;

%% 输入
trackErr1 = 0.1;
trackErr2 = 0.1;
trackErr3 = 0.1;
trackErr4 = 0.1;

trackErrd1 = 0.1;
trackErrd2 = 0.1;
trackErrd3 = 0.1;
trackErrd4 = 0.1;

load('track_circle50.mat');
TRACK = temptrack;

S0 = 0;
Srange = 20;

q = [40, 0, 0, 0, 0];
qd = [3, 0, 0, 0, 0];

qs = [q(1); trackErr1; trackErr2; trackErr3; trackErr4];
qsd = [qd(1); trackErrd1; trackErrd2; trackErrd3; trackErrd4];

%% 根据当前车辆状态获得各控制点投影点的角度 w_psi_psi

% 第一个虚拟铰接点的投影点
[Scur_j1, psi_pj1] = projectPsi(q(1), q(2), S0, Srange, TRACK);

% 四个控制点的投影点

% 首先求出四个控制点的实际位置
s1x = q(1) - cos(q(3)) * l_j1a1;
sly = q(2) - sin(q(3)) * l_j1a1;
s2x = q(1) - cos(q(3)) * l_j1j2;
s2y = q(2) - sin(q(3)) * l_j1j2;
s3x = q(1) - cos(q(3)) * l_j1j2 - cos(q(4)) * l_j2j3;
s3y = q(2) - sin(q(3)) * l_j1j2 - sin(q(4)) * l_j2j3;
s4x = q(1) - cos(q(3)) * l_j1j2 - cos(q(4)) * l_j2j3 - cos(q(5)) * l_j3a8;
s4y = q(2) - sin(q(3)) * l_j1j2 - sin(q(4)) * l_j2j3 - sin(q(5)) * l_j3a8;

% 求投影点
[Scur_s1, psi_ps1, t] = projectPsi(s1x, sly, S0, Srange, TRACK);
[Scur_s2, psi_ps2, ~] = projectPsi(s1x, sly, S0, Srange, TRACK);
[Scur_s3, psi_ps3, ~] = projectPsi(s1x, sly, S0, Srange, TRACK);
[Scur_s4, psi_ps4, ~] = projectPsi(s1x, sly, S0, Srange, TRACK);

%% 求解各单元横摆角的期望部分

ts = zeros(4, 1);
% 把车辆放置在轨迹上（第一个控制点正交投影，其余各点按距离摆放）
Scur_s1_des = Scur_s1;
psi_ps1_des = psi_ps1;
ts(1) = t;
[Scur_s2_des, psi_ps2_des, t] = desprojectPsi(Scur_s1_des, t, L12, Srange, TRACK);
ts(2) = t;
[Scur_s3_des, psi_ps3_des, t] = desprojectPsi(Scur_s2_des, t, L23, Srange, TRACK);
ts(3) = t;
[Scur_s4_des, psi_ps4_des, t] = desprojectPsi(Scur_s3_des, t, L34, Srange, TRACK);
ts(4) = t;

[psi1_des, psi2_des, psi3_des] = psi2ides(Scur_s1_des, Scur_s2_des, Scur_s3_des, Scur_s4_des, ts, TRACK);

%% 根据跟踪误差得到各单元横摆角的误差部分

trackErr_T_psiErr = [cos(psi_pj1 - psi_ps1), l_j1a1 * cos(psi1_des - psi_ps1), 0, 0;
                     cos(psi_pj1 - psi_ps2), l_j1j2 * cos(psi1_des - psi_ps2), 0, 0;
                     cos(psi_pj1 - psi_ps3), l_j1j2 * cos(psi1_des - psi_ps3), l_j2j3 * cos(psi2_des - psi_ps3), 0;
                     cos(psi_pj1 - psi_ps4), l_j1j2 * cos(psi1_des - psi_ps4), l_j2j3 * cos(psi2_des - psi_ps4), l_j3a8 * cos(psi3_des - psi_ps4);];

psi_err = inv(trackErr_T_psiErr) * [qs(2); qs(3); qs(4); qs(5)];

psi1_err = psi_err(2);
psi2_err = psi_err(3);
psi3_err = psi_err(4);

psi1 = psi1_des + psi1_err;
psi2 = psi2_des + psi2_err;
psi3 = psi3_des + psi3_err;

%% 计算得到q、qd、qdd(用误差计算得到)

% 根据psi_ps1、psi_ps2、psi_ps3、psi_ps4、psi1、psi2、psi3得出Js
Js = [1, 0, 0, 0, 0;
      sin(psi1 - psi_ps1), cos(psi1 - psi_ps1), -ls(1) * cos(psi1 - psi_ps1), -ls(2) * cos(psi2 - psi_ps1), -ls(3) * cos(psi3 - psi_ps1);
      sin(psi1 - psi_ps2), cos(psi1 - psi_ps2), -ls(1) * cos(psi1 - psi_ps2), -ls(2) * cos(psi2 - psi_ps2), -ls(3) * cos(psi3 - psi_ps2);
      sin(psi1 - psi_ps3), cos(psi1 - psi_ps3), -ls(1) * cos(psi1 - psi_ps3), -ls(2) * cos(psi2 - psi_ps3), -ls(3) * cos(psi3 - psi_ps3);
      sin(psi1 - psi_ps4), cos(psi1 - psi_ps4), -ls(1) * cos(psi1 - psi_ps4), -ls(2) * cos(psi2 - psi_ps4), -ls(3) * cos(psi3 - psi_ps4);];

% 计算得到q、qd
q_sc = [0; 0; psi1; psi2; psi3];
qd_sc = inv(Js) * qsd;

% 求控制点速度
xd = qd(1);
yd = qd(2);
xd2 = xd + l_j1j2*sin(q(3))*qd(3);
yd2 = yd - l_j1j2*cos(q(3))*qd(3);
xd3 = xd + l_j1j2*sin(q(3))*qd(3) + l_j2j3*sin(q(4))*qd(4);
yd3 = yd - l_j1j2*cos(q(3))*qd(3) - l_j2j3*cos(q(4))*qd(4);

v_s1 = hypot(xd, yd);
v_s2 = hypot(xd, yd);
v_s3 = hypot(xd2, yd2);
v_s4 = hypot(xd3, yd3);

% 求轨道坐标系速度
cigema_ps1 = v_s1 * cos(psi1 - psi_ps1);
cigema_ps2 = v_s2 * cos(psi1 - psi_ps2);
cigema_ps3 = v_s3 * cos(psi2 - psi_ps3);
cigema_ps4 = v_s4 * cos(psi3 - psi_ps4);

psi_1j = qd_sc(3);
psi_2j = qd_sc(4);
psi_3j = qd_sc(5);

% 计算Js的导数
Jsd = [0, 0, 0, 0, 0;
      cos(psi1 - psi_ps1) * (psi_1j - cigema_ps1), -sin(psi1 - psi_ps1) * (psi_1j - cigema_ps1), ls(1) * sin(psi1 - psi_ps1) * (psi_1j - cigema_ps1), ls(2) * sin(psi2 - psi_ps1) * (psi_2j - cigema_ps1), ls(3) * sin(psi3 - psi_ps1) * (psi_3j - cigema_ps1);
      cos(psi1 - psi_ps2) * (psi_1j - cigema_ps2), -sin(psi1 - psi_ps2) * (psi_1j - cigema_ps2), ls(1) * sin(psi1 - psi_ps2) * (psi_1j - cigema_ps2), ls(2) * sin(psi2 - psi_ps2) * (psi_2j - cigema_ps2), ls(3) * sin(psi3 - psi_ps2) * (psi_3j - cigema_ps2);
      cos(psi1 - psi_ps3) * (psi_1j - cigema_ps3), -sin(psi1 - psi_ps3) * (psi_1j - cigema_ps3), ls(1) * sin(psi1 - psi_ps3) * (psi_1j - cigema_ps3), ls(2) * sin(psi2 - psi_ps3) * (psi_2j - cigema_ps3), ls(3) * sin(psi3 - psi_ps3) * (psi_3j - cigema_ps3);
      cos(psi1 - psi_ps4) * (psi_1j - cigema_ps4), -sin(psi1 - psi_ps4) * (psi_1j - cigema_ps4), ls(1) * sin(psi1 - psi_ps4) * (psi_1j - cigema_ps4), ls(2) * sin(psi2 - psi_ps4) * (psi_2j - cigema_ps4), ls(3) * sin(psi3 - psi_ps4) * (psi_3j - cigema_ps4);];

qdd_sc = inv(Js) * qsdd - inv(Js) * Jsd * inv(Js) * qsd; 
