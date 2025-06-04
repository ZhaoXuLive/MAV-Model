clear; clc; close all
 
% 设定图片中的统一字体和字号。单幅图片可以另行设定。设定20号字主要考虑图片缩小后放在论文中的效果。
set(0,'DefaultAxesFontname','Times New Roman');
set(0,'DefaultAxesFontsize',20);

% 建立一个专用的文件夹，用于存放数据和图片。可根据实际情况设定。
FigurePath = ['./SimuData/figures/', datestr(now,'yyyymmdd')];
mkdir(FigurePath)

%% fixed parameter

m1 = 11142;
m3 = 11142;
m2 = 7651; 
l_j1j2 = 11.4;
l_j2j3 = 9.4;
l_j1a1 = 2.5;
l_j1a2 = 4;
l_j1a3 = 9.9;
l_j2a4 = 1.5;
l_j2a5 = 7.9;
l_j3a6 = 1.5;
l_j3a7 = 7.4;
l_j3a8 = 8.9;
l_j1g1 = 5.862; 
l_j2g2 = 4.643; 
l_j3g3 = 5.538; 

%% steady state

R = 50;

delta_dy_v = zeros(70, 8);
delta_ki_v = zeros(70, 8);
for i = 1:70
    v = i / 3.6;
    x = 0;
    y = 0;
    psi1 = asin(4.45 / R);
    psi2 = -asin(4.7 / R);
    psi3 = -(2 * asin(4.7 / R) + asin(4.45 / R));
    fix_q = [x; y; psi1; psi2; psi3];
    xd = v;
    yd = 0;
    psi1d = v/R;
    psi2d = v/R;
    psi3d = v/R;
    fix_qd = [xd; yd; psi1d; psi2d; psi3d];
    delta_d = dynamic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d);
    delta_k = kinetic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d);
    delta_dy_v(i, :) = delta_d' * 180 / pi;
    delta_ki_v(i, :) = delta_k' * 180 / pi;
end

delta_dy_r = zeros(15, 8);
delta_ki_r = zeros(15, 8);
for i = 25:5:100
    v = 50 / 3.6;
    x = 0;
    y = 0;
    psi1 = asin(4.45 / i);
    psi2 = -asin(4.7 / i);
    psi3 = -(2 * asin(4.7 / i) + asin(4.45 / i));
    fix_q = [x; y; psi1; psi2; psi3];
    xd = v;
    yd = 0;
    psi1d = v/i;
    psi2d = v/i;
    psi3d = v/i;
    fix_qd = [xd; yd; psi1d; psi2d; psi3d];
    delta_d = dynamic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d);
    delta_k = kinetic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d);
    delta_dy_r(i / 5 - 4, :) = delta_d' * 180 / pi;
    delta_ki_r(i / 5 - 4, :) = delta_k' * 180 / pi;
end

%% plot

% v change
h1 = figure;
v_list = 1:1:70;
plot(v_list, delta_dy_v(:, 1), '--');
hold on
plot(v_list, delta_ki_v(:, 1));
xlabel('Longitudinal speed [km/h]');
ylabel('Steer angle [deg]');
legend('dynamic', 'kinetic');
axis([0 70 -4 -2])
print(h1,'-depsc2','-loose',[FigurePath,'speed_a.eps']); 
savefig(h1, [FigurePath,'/figure1.fig']);


h2 = figure;
plot(v_list, delta_dy_v(:, 4), '--');
hold on
plot(v_list, delta_ki_v(:, 4));
xlabel('Longitudinal speed [km/h]');
ylabel('Steer angle [deg]');
legend('dynamic', 'kinetic');
axis([0 70 -5 -3])
print(h2,'-depsc2','-loose',[FigurePath,'speed_b.eps']); 
savefig(h2, [FigurePath,'/figure2.fig']);

h3 = figure;
plot(v_list, delta_dy_v(:, 5), '--');
hold on
plot(v_list, delta_ki_v(:, 5));
xlabel('Longitudinal speed [km/h]');
ylabel('Steer angle [deg]');
legend('dynamic', 'kinetic');
axis([0 70 -12 -10])
print(h3,'-depsc2','-loose',[FigurePath,'speed_c.eps']); 
savefig(h3, [FigurePath,'/figure3.fig']);

h4 = figure;
plot(v_list, delta_dy_v(:, 8), '--');
hold on
plot(v_list, delta_ki_v(:, 8));
xlabel('Longitudinal speed [km/h]');
ylabel('Steer angle [deg]');
legend('dynamic', 'kinetic');
axis([0 70 -12 -10])
print(h4,'-depsc2','-loose',[FigurePath,'speed_d.eps']); 
savefig(h4, [FigurePath,'/figure4.fig']);

% r change
h5 = figure;
r_list = 25:5:100;
plot(r_list, delta_dy_r(:, 1), '--');
hold on
plot(r_list, delta_ki_r(:, 1));
xlabel('Turning Radius [m]');
ylabel('Steer Angle [deg]');
legend('dynamic', 'kinetic');
axis([25 100 -8 2])
print(h5,'-depsc2','-loose',[FigurePath,'radius_a.eps']); 
savefig(h5, [FigurePath,'/figure5.fig']);

h6 = figure;
plot(r_list, delta_dy_r(:, 4), '--');
hold on
plot(r_list, delta_ki_r(:, 4));
xlabel('Turning Radius [m]');
ylabel('Steer Angle [deg]');
legend('dynamic', 'kinetic');
axis([25 100 -8 2])
print(h6,'-depsc2','-loose',[FigurePath,'radius_b.eps']); 
savefig(h6, [FigurePath,'/figure6.fig']);

h7 = figure;
plot(r_list, delta_dy_r(:, 5), '--');
hold on
plot(r_list, delta_ki_r(:, 5));
xlabel('Turning Radius [km/h]');
ylabel('Steer Angle [deg]');
legend('dynamic', 'kinetic');
axis([25 100 -18 -2])
print(h7,'-depsc2','-loose',[FigurePath,'radius_c.eps']); 
savefig(h7, [FigurePath,'/figure7.fig']);

h8 = figure;
plot(r_list, delta_dy_r(:, 8), '--');
hold on
plot(r_list, delta_ki_r(:, 8));
xlabel('Turning Radius [m]');
ylabel('Steer Angle [deg]');
legend('dynamic', 'kinetic');
axis([25 100 -18 -2])
print(h8,'-depsc2','-loose',[FigurePath,'radius_d.eps']); 
savefig(h8, [FigurePath,'/figure8.fig']);


%% dynamic model
function delta_d = dynamic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d)

m1 = 11142;
m3 = 11142;
m2 = 7651; 
l_j1j2 = 11.4;
l_j2j3 = 9.4;
l_j1a1 = 2.5;
l_j1a2 = 4;
l_j1a3 = 9.9;
l_j2a4 = 1.5;
l_j2a5 = 7.9;
l_j3a6 = 1.5;
l_j3a7 = 7.4;
l_j3a8 = 8.9;
l_j1g1 = 5.862; 
l_j2g2 = 4.643; 
l_j3g3 = 5.538; 

% trans
xd_v = xd;
yd_v = yd;
xd_w = xd_v * cos(psi1) - yd_v * sin(psi1);
yd_w = xd_v * sin(psi1) + yd_v * cos(psi1);
xd = xd_w;
yd = yd_w;

global L_long Ca L_lat beta C D;

C = [0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, l_j3g3*m3*cos(psi3)*psi3d;
     0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, l_j3g3*m3*sin(psi3)*psi3d;
     1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, 1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*(xd*cos(psi1)+yd*sin(psi1))+1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi2d+1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi3d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, 1/2*(l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(l_j2g2*m2+l_j2j3*m3)*(xd*cos(psi2)+yd*sin(psi2))-1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi1d+1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi3d, l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*l_j3g3*m3*cos(psi3)*psi3d, 1/2*l_j3g3*m3*sin(psi3)*psi3d, l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-psi1d+1/2*psi3d), l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*l_j3g3*m3*(xd*cos(psi3)+yd*sin(psi3))-1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi1d-1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi2d;];

d1 = 10000;
d2 = 10000;

delta_kj = eye(3);
dj = [0, d1, d2, 0]';
D = zeros(5);

for k = 1:5
    for l = 1:5
        if k > 2 && l > 2
            if  l == k
                D(k, l) = (1 - delta_kj(1, k-2))*dj(k - 2) + (1 - delta_kj(3, k-2))*dj(k - 1);
            elseif l == k-1
                D(k, l) = -(1 - delta_kj(1, k-2))*dj(k - 2);
            elseif l == k+1
                D(k, l) = -(1 - delta_kj(3, k-2))*dj(k - 1);
            else
                D(k, l) = 0;
            end
        end
    end
end

psij = [psi1, psi2, psi3]';
Rot_w2v = zeros(2, 2, 3);
for i = 1:3
    Rot_w2v(:, :, i) = [cos(psij(i)), -sin(psij(i));
                        sin(psij(i)), cos(psij(i));];
end

u = zeros(2, 3);
n = zeros(2, 3);
for i = 1:3
    u(:, i)=[cos(psij(i)), sin(psij(i))]';
    n(:, i)=[-sin(psij(i)), cos(psij(i))]';
end

IU = [1, 1, 1, 2, 2, 3, 3, 3]';
l_ja = [l_j1a1, l_j1a2, l_j1a3, l_j2a4, l_j2a5, l_j3a6, l_j3a7, l_j3a8]';

L_long = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           0, 0, 0, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, 0, 0, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, 0, 0, 0;];

L_lat = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;];

vxj = zeros(3, 1);
vyj = zeros(3, 1);
vxj(1) = xd;
vyj(1) = yd;
vxj(2) = xd + l_j1j2*sin(psi1)*psi1d;
vyj(2) = yd - l_j1j2*cos(psi1)*psi1d;
vxj(3) = xd + l_j1j2*sin(psi1)*psi1d + l_j2j3*sin(psi2)*psi2d;
vyj(3) = yd - l_j1j2*cos(psi1)*psi1d - l_j2j3*cos(psi2)*psi2d;

va = zeros(2, 10);
psidj = [psi1d, psi2d, psi3d]';
for i = 1:8
    va(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i));
end

va_v = zeros(2, 10);
for i = 1:8
    va_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va(:, i);
end

beta = zeros(8, 1);
for i = 1:8
    beta(i) = atan2(va_v(2, i), va_v(1, i));
end

ca = [400e3, 400e3, 400e3, 400e3, 400e3, 400e3, 400e3, 400e3]';
Ca = diag(ca);

global dynamic_qd;
dynamic_qd = [xd; yd; psi1d; psi2d; psi3d];

options = optimoptions('fsolve');
options.Algorithm = 'levenberg-marquardt';
initial_delta = zeros(8, 1);
delta_d = fsolve(@dynamic_func, initial_delta, options);

end

function F = dynamic_func(delta)
    global L_long Ca L_lat beta C D; 
    global dynamic_qd;
    F = L_long * (- diag(sin(delta)) * Ca * (delta - beta)) + L_lat * (diag(cos(delta)) * Ca * (delta - beta)) - (C + D) * dynamic_qd;
end

%% kinetic model
function delta_k = kinetic_model(psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d)

l_j1j2 = 11.4;
l_j2j3 = 9.4;
l_j1a1 = 2.5;
l_j1a2 = 4;
l_j1a3 = 9.9;
l_j2a4 = 1.5;
l_j2a5 = 7.9;
l_j3a6 = 1.5;
l_j3a7 = 7.4;
l_j3a8 = 8.9;

% trans
xd_v = xd;
yd_v = yd;
xd_w = xd_v * cos(psi1) - yd_v * sin(psi1);
yd_w = xd_v * sin(psi1) + yd_v * cos(psi1);
xd = xd_w;
yd = yd_w;

psij = [psi1, psi2, psi3]';
Rot_w2v = zeros(2, 2, 3);
for i = 1:3
    Rot_w2v(:, :, i) = [cos(psij(i)), -sin(psij(i));
                        sin(psij(i)), cos(psij(i));];
end

u = zeros(2, 3);
n = zeros(2, 3);
for i = 1:3
    u(:, i)=[cos(psij(i)), sin(psij(i))]';
    n(:, i)=[-sin(psij(i)), cos(psij(i))]';
end

IU = [1, 1, 1, 2, 2, 3, 3, 3]';
l_ja = [l_j1a1, l_j1a2, l_j1a3, l_j2a4, l_j2a5, l_j3a6, l_j3a7, l_j3a8]';

vxj = zeros(3, 1);
vyj = zeros(3, 1);
vxj(1) = xd;
vyj(1) = yd;
vxj(2) = xd + l_j1j2*sin(psi1)*psi1d;
vyj(2) = yd - l_j1j2*cos(psi1)*psi1d;
vxj(3) = xd + l_j1j2*sin(psi1)*psi1d + l_j2j3*sin(psi2)*psi2d;
vyj(3) = yd - l_j1j2*cos(psi1)*psi1d - l_j2j3*cos(psi2)*psi2d;

va = zeros(2, 10);
psidj = [psi1d, psi2d, psi3d]';
for i = 1:8
    va(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i));
end

va_v = zeros(2, 10);
for i = 1:8
    va_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va(:, i);
end

beta = zeros(8, 1);
for i = 1:8
    beta(i) = atan2(va_v(2, i), va_v(1, i));
end

delta_k = beta;

end

%% newton model
function delta_n = newton_model(q_input, qd_input, V)

% the number of vehicle units
unit_n = 3;
% the number of vehicle axles
axle_m = 8;
% the number of axles of every vehicle unit 
getNumAxle = [0, 3, 5, 8];
% the unit of the respective axle
getNumUnit = [1, 1, 1, 2, 2, 3, 3, 3];

% the mass of every vehicle unit
m = [11142, 7651, 11142];
% the rotational inertia of every vehicle unit
Iz = [143570, 37800, 143570];
% the tire cornering stiffness
Ca = 3065 * 180 / pi;
% the articulated angle damping
art_d = [10000, 10000];
% art_d = [0, 0];

% the length of every vehicle unit
lu = [11.4, 9.4, 11.4];
% the length of every vehicle unit from the front joint to the mass center
lg = [5.862, 4.643, 5.538];
% the length of every vehicle unit from the front joint to the axle center
la = [2.5, 4, 9.9, 1.5, 7.9, 1.5, 7.4, 8.9];

% calculate U \beta r \theta
psi = q_input(3:unit_n+2);
dpsi = qd_input(3:unit_n+2);

% trans
xd_v = qd_input(1);
yd_v = qd_input(2);
xd_w = xd_v * cos(psi(1)) - yd_v * sin(psi(1));
yd_w = xd_v * sin(psi(1)) + yd_v * cos(psi(1));
xd = xd_w;
yd = yd_w;

qd_input(1) = xd;
qd_input(2) = yd;

r = dpsi;
theta = zeros(unit_n-1, 1);
for i = 1:unit_n-1
    theta(i) = psi(i) - psi(i+1);
end

v_gw = zeros(2, unit_n);
v_gv = zeros(2, unit_n);
beta = zeros(unit_n, 1);

v_jw(:, 1) = qd_input(1:2);
for i = 1:unit_n-1
    v_jw(1, i+1) = v_jw(1, i) + sin(psi(i)) * lu(i) * dpsi(i);
    v_jw(2, i+1) = v_jw(2, i) - cos(psi(i)) * lu(i) * dpsi(i);
end

for i = 1:unit_n
    v_gw(1, i) = v_jw(1, i) + sin(psi(i)) * lg(i) * dpsi(i);
    v_gw(2, i) = v_jw(2, i) - cos(psi(i)) * lg(i) * dpsi(i);
    v_gv(1, i) = v_gw(1, i) * cos(psi(i)) + v_gw(2, i) * sin(psi(i));
    v_gv(2, i) = -v_gw(1, i) * sin(psi(i)) + v_gw(2, i) * cos(psi(i));
    U = v_gv(1, 1);
    beta(i) = v_gv(2, i) / U;
end

L = cell(1, unit_n);
for k = 1:unit_n
    % L1: 鍚勮酱涓績鍒板崟鍏冭川蹇冪殑璺濈锛屽墠涓烘锛屽悗涓鸿礋
    L1 = lg(k) - la(getNumUnit == k);
    L{k} = L1;
end

nx = 3*unit_n-1; % number of states
nv = 4*unit_n-2; % number of variables (including joint forces)

% S * (\dot{x} R)' + P x = By y 鏂圭▼
% x = {\beta r \theta}  n n n-1 鐘舵?侀噺
% r = \dot{\psi} 妯憜瑙掗?熷害
% \theta = \psi_j - \psi_j+1  閾版帴瑙?

% dimensions of primary matrices
S = zeros(nv, nv);
P = zeros(nv, nx);
BY = zeros(nv, axle_m);
K1 = zeros(axle_m, nx);
K2 = zeros(axle_m, axle_m);

% joint damp
% input matrix for external yaw moments
BM = zeros(nv, unit_n); 

% populate S,P,B equation-by-equation

% first define end-index of sub-matrix blocks
% 3n-1涓姸鎬侊紝鍓峮涓槸璐ㄥ績渚у亸瑙掞紝涓棿n涓槸妯憜瑙掗?熷害锛屾渶鍚巒-1涓槸閾版帴瑙?
b1 = unit_n;         % last row of beta_dot
b2 = 2 * unit_n;     % last row of r_dot
b3 = 3 * unit_n - 1; % last row of theta_dot

% also need number of 'preceding' axles at each unit
% 姣忎釜鍗曞厓鍓嶉潰鎵?鏈夌殑鍗曞厓鍏辨湁澶氬皯涓酱
pax = zeros(1, unit_n); % preceding axle numbers
for i = 2:unit_n
    pax(i) = pax(i - 1) + length(L{i - 1});
end

% equation (1)
for i = 1:unit_n 
    S(i, i) = m(i) * U;
    P(i, b1 + i) = m(i) * U;
    if i > 1
        % R_(i - 1) term in equ(1)
        S(i, b3 + i - 1) = -1;
    end 
    if i < unit_n
        % R_i term 
        S(i, b3 + i) = 1;
    end 
    for j = 1:length(L{i})
        % Yij terms
       BY(i, pax(i) + j) = 1; 
    end
end

% equation (2)
for i = 1:unit_n 
    S(b1 + i, b1 + i) = Iz(i);
    if i > 1
        % R_(i-1) term in equ(1)
        S(b1 + i, b3 + i - 1) = -lg(i);
    end 
    if i < unit_n
        % R_i term 
        S(b1 + i, b3 + i) = lg(i) - lu(i);
    end 
    for j = 1:length(L{i})
        % Mij terms
       BY(b1 + i, pax(i) + j) = L{i}(j); 
    end
    BM(b1+i) = 1; % direct yaw moment term
end

% equation (3)
for i = 1:unit_n-1
    % theta-dot term
    S(b2 + i, b2 + i) = 1; 
    % r_i in notes = unit in front of the joint
    P(b2 + i, b1 + i) = -1; 
    % r_(i+1)
    P(b2 + i, b1 + i + 1) = 1; 
end

% equation (4) 
for i = 1:unit_n-1 
    S(b3 + i, i + 1) = U;
    S(b3 + i, i) = -U; 
    S(b3 + i, b1 + i + 1) = lg(i + 1);
    S(b3 + i, b1 + i) = -(lg(i) - lu(i));
    S(b3 + i, b2 + i) = -U; 
end

% K1 matrix, Y = K1 * X for each axle force
for i = 1:unit_n
    for j = 1:length(L{i})
        % 'flat' axle index
        ia = pax(i) + j; 
        K1(ia, i) = U;
        K1(ia, b1 + i) = L{i}(j);
        K1(ia, :) = -K1(ia, :) * Ca / U;
    end
end

% K2 matrix, go axle by axle
for i = 1:unit_n
    for j = 1:length(L{i})
        % 'flat' axle index
        ia = pax(i) + j; 
        K2(ia, ia) = Ca;
    end
end

AA = S \ (BY * K1 - P);
A = AA(1:nx, :);

BB = S \ BY * K2;
B = BB(1:nx, :);

% /dot x = AX +B /delta
x = [beta; r; theta];

dx = zeros(8, 1);

% dr -- ddpsi
dx(4:6) = zeros(3, 1);
% dtheta -- dpsi - dpsi+1
dx(7) = qd_input(3) - qd_input(4);
dx(8) = qd_input(4) - qd_input(5);
% dbeta -- 
% new qdd to dx
qdd_v = zeros(2, 1);
dx(1) = (qdd_v(2) - lg(1) * dx(4)) / V;
dx(2) = (qdd_v(2) - lu(1) * dx(4) - lg(2) * dx(5)) / V;

delta = B \ (dx - A * x);

end