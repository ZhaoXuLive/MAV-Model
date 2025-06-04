function qdd_v = linear_variable(q_v, qd_v, delta_v)

%% parameter
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

h11 = 0.9625;
h12 = 0.9625;
h13 = 0.98;
h21 = 0.98; 
h22 = 0.98;
h31 = 0.98;
h32 = 0.9625;
h33 = 0.9625;

ca1 = 3065 * 180 / pi;
ca2 = 3065 * 180 / pi;
ca3 = 3065 * 180 / pi;
ca4 = 3065 * 180 / pi;
ca5 = 3065 * 180 / pi;
ca6 = 3065 * 180 / pi;
ca7 = 3065 * 180 / pi;
ca8 = 3065 * 180 / pi;

Iz1 = 143570;
Iz3 = 143570;
Iz2 = 37800;

d1 = 10000;
d2 = 10000;

% input
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real

% state variable
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [xd; yd; psi1d; psi2d; psi3d;];
delta = [delta1; delta2; delta3; delta4; delta5; delta6; delta7; delta8;];

%% matrix

% M
l_kj = [l_j1g1, l_j1j2, l_j1j2;
             0, l_j2g2, l_j2j3;
             0,      0, l_j3g3;];
delta_kj = eye(3);
Jc = sym(zeros(3, 5, 3));
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

c33 = -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3)) * ( (xd*cos(psi1)-yd*sin(psi1)-xd) * cos(psi1) + (xd*sin(psi1)+yd*cos(psi1)-yd) * sin(psi1) );
c44 = -1/2*(l_j2g2*m2 + l_j2j3*m3) * ( (xd*cos(psi1)-yd*sin(psi1)-xd) * cos(psi2) + (xd*sin(psi1)+yd*cos(psi1)-yd) * sin(psi2) );
c55 = -1/2*l_j3g3*m3 * ( (xd*cos(psi1)-yd*sin(psi1)-xd) * cos(psi3) + (xd*sin(psi1)+yd*cos(psi1)-yd) * sin(psi3) );
C_wv = [0, 0, 0, 0, 0;
        0, 0 ,0, 0, 0;
        0, 0, c33, 0, 0;
        0, 0, 0, c44, 0;
        0, 0, 0, 0, c55;];

C = C + C_wv;

% D
D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, d1, -d1, 0;
     0, 0, -d1, d1+d2, -d2;
     0, 0, 0, -d2, d2;];

% Q
psij = [psi1, psi2, psi3]';
Rot_w2v = sym(zeros(2, 2, 3));
for i = 1:3
    Rot_w2v(:, :, i) = [cos(psij(i)), -sin(psij(i));
                        sin(psij(i)), cos(psij(i));];
end

Rot_v2t = sym(zeros(2, 2, 8));
for i = 1:8
    Rot_v2t(:, :, i) = [cos(delta(i)), -sin(delta(i));
                        sin(delta(i)), cos(delta(i));];
end

xj = sym(zeros(3, 1));
yj = sym(zeros(3, 1));
xj(1) = x;
yj(1) = y;
xj(2) = x - l_j1j2*cos(psi1);
yj(2) = y - l_j1j2*sin(psi1);
xj(3) = x - l_j1j2*cos(psi1) - l_j2j3*cos(psi2);
yj(3) = y - l_j1j2*sin(psi1) - l_j2j3*sin(psi2);

u = sym(zeros(2, 3));
n = sym(zeros(2, 3));
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

Delta_c = sym(zeros(8, 8));
Delta_s = sym(zeros(8, 8));
for i = 1:8
    Delta_c(i, i) = cos(delta(i));
    Delta_s(i, i) = sin(delta(i));
end

vxj = sym(zeros(3, 1));
vyj = sym(zeros(3, 1));
vxj(1) = xd;
vyj(1) = yd;
vxj(2) = xd + l_j1j2*sin(psi1)*psi1d;
vyj(2) = yd - l_j1j2*cos(psi1)*psi1d;
vxj(3) = xd + l_j1j2*sin(psi1)*psi1d + l_j2j3*sin(psi2)*psi2d;
vyj(3) = yd - l_j1j2*cos(psi1)*psi1d - l_j2j3*cos(psi2)*psi2d;
va1 = sym(zeros(2, 8));
va2 = sym(zeros(2, 8));
psidj = [psi1d, psi2d, psi3d]';
for i = 1:8
    va1(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i)) + (-1)*u(:, IU(i))*ha(i)*psidj(IU(i));
    va2(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i)) + u(:, IU(i))*ha(i)*psidj(IU(i));
end

va1_v = sym(zeros(2, 8));
va2_v = sym(zeros(2, 8));
for i = 1:8
    va1_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va1(:, i);
    va2_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va2(:, i);
end

beta_wv_y = sym(zeros(3, 1));
beta_wv_y(1) = (-(xd*cos(psi1) - yd*sin(psi1))*sin(psi1) + (xd*sin(psi1) + yd*cos(psi1))*cos(psi1)) - (-xd * sin(psi1) + yd * cos(psi1));
beta_wv_y(2) = (-(xd*cos(psi1) - yd*sin(psi1))*sin(psi2) + (xd*sin(psi1) + yd*cos(psi1))*cos(psi2)) - (-xd * sin(psi2) + yd * cos(psi2));
beta_wv_y(3) = (-(xd*cos(psi1) - yd*sin(psi1))*sin(psi3) + (xd*sin(psi1) + yd*cos(psi1))*cos(psi3)) - (-xd * sin(psi3) + yd * cos(psi3));
beta_wv_x = sym(zeros(3, 1));
beta_wv_x(1) = ((xd*cos(psi1) - yd*sin(psi1))*cos(psi1) + (xd*sin(psi1) + yd*cos(psi1))*sin(psi1)) - (xd * cos(psi1) + yd * sin(psi1));
beta_wv_x(2) = ((xd*cos(psi1) - yd*sin(psi1))*cos(psi2) + (xd*sin(psi1) + yd*cos(psi1))*sin(psi2)) - (xd * cos(psi2) + yd * sin(psi2));
beta_wv_x(3) = ((xd*cos(psi1) - yd*sin(psi1))*cos(psi3) + (xd*sin(psi1) + yd*cos(psi1))*sin(psi3)) - (xd * cos(psi3) + yd * sin(psi3));

beta1 = sym(zeros(8, 1));
beta2 = sym(zeros(8, 1));
for i = 1:8
    beta1(i) = atan2(va1_v(2, i) + beta_wv_y(IU(i)), va1_v(1, i) + beta_wv_x(IU(i)));
    beta2(i) = atan2(va2_v(2, i) + beta_wv_y(IU(i)), va2_v(1, i) + beta_wv_x(IU(i)));
end

ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8]';
Ca = diag(ca);

Tau_tire = L_long1 * (- Delta_s * Ca * (delta - beta1)) ...
         + L_long2 * (- Delta_s * Ca * (delta - beta2)) ...
         + L_lat1 * (Delta_c * Ca * (delta - beta1)) ...
         + L_lat2 * (Delta_c * Ca * (delta - beta2));

%% trans
Twv = [cos(psi1), -sin(psi1), 0, 0, 0;
     sin(psi1), cos(psi1), 0, 0, 0;
     0, 0, 1, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1;];

Tdwv = [-sin(psi1)*psi1d, -cos(psi1)*psi1d, 0, 0, 0;
        cos(psi1)*psi1d, -sin(psi1)*psi1d, 0, 0, 0;
        0, 0, 0, 0, 0;
        0, 0, 0, 0, 0;
        0, 0, 0, 0, 0;];

Mv = Twv' * M * Twv;
Cv = Twv' * C * Twv + Twv' * M * Tdwv;
Dv = Twv' * D * Twv;
Qv = Twv' * Tau_tire;

% main equation
qdd = inv(Mv)*(Qv - Cv * qd - Dv * qd);
qdd(1) = 0;

A1 = jacobian(qdd, q);
B1 = jacobian(qdd, qd);
C1 = jacobian(qdd, delta);

%% sub

A1l = subs(A1,{x, y, psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d, ...
    delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8} ...
    , {q_v(1), q_v(2), q_v(3), q_v(4), q_v(5), qd_v(1), qd_v(2), qd_v(3), qd_v(4), qd_v(5), ...
    delta_v(1), delta_v(2), delta_v(3), delta_v(4), delta_v(5), delta_v(6), delta_v(7), delta_v(8)});
A2l = vpa(A1l, 3);

B1l = subs(B1,{x, y, psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d, ...
    delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8} ...
    , {q_v(1), q_v(2), q_v(3), q_v(4), q_v(5), qd_v(1), qd_v(2), qd_v(3), qd_v(4), qd_v(5), ...
    delta_v(1), delta_v(2), delta_v(3), delta_v(4), delta_v(5), delta_v(6), delta_v(7), delta_v(8)});
B2l = vpa(B1l, 3);

C1l = subs(C1,{x, y, psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d, ...
    delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8} ...
    , {q_v(1), q_v(2), q_v(3), q_v(4), q_v(5), qd_v(1), qd_v(2), qd_v(3), qd_v(4), qd_v(5), ...
    delta_v(1), delta_v(2), delta_v(3), delta_v(4), delta_v(5), delta_v(6), delta_v(7), delta_v(8)});
C2l = vpa(C1l, 3);

qdd_v = simplify(A2l * q_v + B2l * qd_v + C2l * delta_v);
qdd_v = double(vpa(qdd_v, 6));

end
