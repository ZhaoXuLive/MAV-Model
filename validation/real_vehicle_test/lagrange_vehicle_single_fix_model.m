function qdd_output = lagrange_vehicle_single_fix_model(q, qd, delta, V, sat)

% fixed parameter
m1 = 11142;
m3 = 11142; 
m2 = 7651; 

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

h11 = 0;
h12 = 0;
h13 = 0;
h21 = 0; 
h22 = 0;
h31 = 0;
h32 = 0;
h33 = 0;

ca_T = 3065 * 180 / pi;
ca1 = ca_T;
ca2 = ca_T;
ca3 = ca_T;
ca4 = ca_T;
ca5 = ca_T;
ca6 = ca_T;
ca7 = ca_T;
ca8 = ca_T;

Iz1 = 143570;
Iz3 = 143570;
Iz2 = 37800;

d1 = 10000;
d2 = 10000;

% trans
x = q(1);
y = q(2);
psi1 = q(3);
psi2 = q(4);
psi3 = q(5);
xd = qd(1);
yd = qd(2);
psi1d = qd(3);
psi2d = qd(4);
psi3d = qd(5);

% disp(qd);
% change V
xd_v = V;
yd_v = - xd * sin(psi1) + yd * cos(psi1);
xd_w = xd_v * cos(psi1) - yd_v * sin(psi1);
yd_w = xd_v * sin(psi1) + yd_v * cos(psi1);
xd = xd_w;
yd = yd_w;
qd(1) = xd;
qd(2) = yd;
% disp(qd);

%% matrix

% M
l_kj = [l_j1g1, l_j1j2, l_j1j2;
             0, l_j2g2, l_j2j3;
             0,      0, l_j3g3;];
delta_kj = eye(3);
Jc = zeros(3, 5, 3);
for j = 1:3
    Jc(:, :, j) = [1, 0,  l_kj(1, j) * sin(psi1),  l_kj(2, j) * sin(psi2),  l_kj(3, j) * sin(psi3);
          0, 1, -l_kj(1, j) * cos(psi1), -l_kj(2, j) * cos(psi2), -l_kj(3, j) * cos(psi3);
          0, 0,          delta_kj(1, j),          delta_kj(2, j),          delta_kj(3, j);];
end

mj = [m1, m2, m3]';
izj = [Iz1, Iz2, Iz3]';
Mj = zeros(3, 3, 3);
for j = 1:3
    Mj(:, :, j) = diag([mj(j), mj(j), izj(j)]);
end

M = zeros(5, 5);
for j = 1:3
    M = M + Jc(:, :, j)' * Mj(:, :, j) * Jc(:, :, j);
end

% C
C = [0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, l_j3g3*m3*cos(psi3)*psi3d;
     0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, l_j3g3*m3*sin(psi3)*psi3d;
     1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, 1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*(xd*cos(psi1)+yd*sin(psi1))+1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi2d+1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi3d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, 1/2*(l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(l_j2g2*m2+l_j2j3*m3)*(xd*cos(psi2)+yd*sin(psi2))-1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi1d+1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi3d, l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*l_j3g3*m3*cos(psi3)*psi3d, 1/2*l_j3g3*m3*sin(psi3)*psi3d, l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-psi1d+1/2*psi3d), l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*l_j3g3*m3*(xd*cos(psi3)+yd*sin(psi3))-1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi1d-1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi2d;];

% D
D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, d1, -d1, 0;
     0, 0, -d1, d1+d2, -d2;
     0, 0, 0, -d2, d2;];

% Q
psij = [psi1, psi2, psi3]';
Rot_w2v = zeros(2, 2, 3);
for i = 1:3
    Rot_w2v(:, :, i) = [cos(psij(i)), -sin(psij(i));
                        sin(psij(i)), cos(psij(i));];
end

Rot_v2t = zeros(2, 2, 8);
for i = 1:8
    Rot_v2t(:, :, i) = [cos(delta(i)), -sin(delta(i));
                        sin(delta(i)), cos(delta(i));];
end

u = zeros(2, 3);
n = zeros(2, 3);
for i = 1:3
    u(:, i)=[cos(q(i + 2)), sin(q(i + 2))]';
    n(:, i)=[-sin(q(i + 2)), cos(q(i + 2))]';
end

IU = [1, 1, 1, 2, 2, 3, 3, 3]';
l_ja = [l_j1a1, l_j1a2, l_j1a3, l_j2a4, l_j2a5, l_j3a6, l_j3a7, l_j3a8]';
ha = [h11, h12, h13, h21, h22, h31, h32, h33]';

L_long1 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           -h11, -h12, -h13, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, -h21, -h22, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, -h31, -h32, -h33;];

L_long2 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           h11, h12, h13, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, h21, h22, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, h31, h32, h33;];

L_lat1 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;];

L_lat2 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;];

Delta_c = zeros(8, 8);
Delta_s = zeros(8, 8);
for i = 1:8
    Delta_c(i, i) = cos(delta(i));
    Delta_s(i, i) = sin(delta(i));
end

vxj = zeros(3, 1);
vyj = zeros(3, 1);
vxj(1) = xd;
vyj(1) = yd;
vxj(2) = xd + l_j1j2*sin(psi1)*psi1d;
vyj(2) = yd - l_j1j2*cos(psi1)*psi1d;
vxj(3) = xd + l_j1j2*sin(psi1)*psi1d + l_j2j3*sin(psi2)*psi2d;
vyj(3) = yd - l_j1j2*cos(psi1)*psi1d - l_j2j3*cos(psi2)*psi2d;
va1 = zeros(2, 8);
va2 = zeros(2, 8);
psidj = [psi1d, psi2d, psi3d]';
for i = 1:8
    va1(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i)) + (-1)*u(:, IU(i))*ha(i)*psidj(IU(i));
    va2(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i)) + u(:, IU(i))*ha(i)*psidj(IU(i));
end

va1_v = zeros(2, 8);
va2_v = zeros(2, 8);
for i = 1:8
    va1_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va1(:, i);
    va2_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va2(:, i);
end

beta1 = zeros(8, 1);
beta2 = zeros(8, 1);
for i = 1:8
    beta1(i) = atan2(va1_v(2, i), va1_v(1, i));
    beta2(i) = atan2(va2_v(2, i), va2_v(1, i));
end

ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8]';
Ca = diag(ca);

alpha1 = delta' - beta1;
alpha2 = delta' - beta2;
for i = 1:8
    if abs(alpha1(i)) > sat * pi / 180
        alpha1(i) = sign(alpha1(i)) * sat * pi / 180;
    end
    if abs(alpha2(i)) > sat * pi / 180
        alpha2(i) = sign(alpha2(i)) * sat * pi / 180;
    end
end

Tau_tire = L_long1 * (- Delta_s * Ca * alpha1) ...
         + L_long2 * (- Delta_s * Ca * alpha2) ...
         + L_lat1 * (Delta_c * Ca * alpha1) ...
         + L_lat2 * (Delta_c * Ca * alpha2);

Q = Tau_tire;

qdd_output = inv(M)*(Q - C * qd - D * qd);
