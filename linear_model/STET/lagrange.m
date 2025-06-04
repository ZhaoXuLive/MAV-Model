% q = [1;2;3;4;5];
% qd = [6;7;8;9;10;];
% delta = [pi / 6;pi / 6;pi / 6;pi / 6;pi / 6;pi / 6;pi / 6;pi / 6];
% 
% qdd = linear_variable(q, qd, delta)

% fixed parameter
syms m1 m2 m3 real
% m1 = 11142;
% m3 = 11142; 
% m2 = 7651; 

syms l_j1j2 l_j2j3 l_j1g1 l_j2g2 l_j3g3 l_j1a1 l_j1a2 l_j1a3 l_j2a4 l_j2a5 l_j3a6 l_j3a7 l_j3a8 real
% l_j1j2 = 11.4;
% l_j2j3 = 9.4;
% l_j1g1 = 11.4 - 5.538; 
% l_j2g2 = 9.4 - 4.757;
% l_j3g3 = 5.538; 
% l_j1a1 = 2.5; 
% l_j1a2 = 4;
% l_j1a3 = 9.9;
% l_j2a4 = 1.5;
% l_j2a5 = 7.9;
% l_j3a6 = 1.5;
% l_j3a7 = 7.4;
% l_j3a8 = 8.9;

syms h11 h12 h13 h21 h22 h31 h32 h33 real
% h11 = 0.9625;
% h12 = 0.9625;
% h13 = 0.98;
% h21 = 0.98; 
% h22 = 0.98;
% h31 = 0.98;
% h32 = 0.9625;
% h33 = 0.9625;

syms ca1 ca2 ca3 ca4 ca5 ca6 ca7 ca8 real
% ca1 = 3065 * 180 / pi;
% ca2 = 3065 * 180 / pi;
% ca3 = 3065 * 180 / pi;
% ca4 = 3065 * 180 / pi;
% ca5 = 3065 * 180 / pi;
% ca6 = 3065 * 180 / pi;
% ca7 = 3065 * 180 / pi;
% ca8 = 3065 * 180 / pi;

syms Iz1 Iz3 Iz2 d1 d2 real
% Iz1 = 143570;
% Iz3 = 143570;
% Iz2 = 37800;
% 
% d1 = 10000;
% d2 = 10000;

% input
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real

% state variable
syms x y psi1 psi2 psi3 real
syms U yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [U; yd; psi1d; psi2d; psi3d;];
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
Mj = sym(zeros(3, 3, 3));
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
     1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, 1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*(U*cos(psi1)+yd*sin(psi1))+1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi2d+1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi3d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, 1/2*(l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(l_j2g2*m2+l_j2j3*m3)*(U*cos(psi2)+yd*sin(psi2))-1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi1d+1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi3d, l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*l_j3g3*m3*cos(psi3)*psi3d, 1/2*l_j3g3*m3*sin(psi3)*psi3d, l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-psi1d+1/2*psi3d), l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*l_j3g3*m3*(U*cos(psi3)+yd*sin(psi3))-1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi1d-1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi2d;];

c33 = -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3)) * ( (U*cos(psi1)-yd*sin(psi1)-U) * cos(psi1) + (U*sin(psi1)+yd*cos(psi1)-yd) * sin(psi1) );
c44 = -1/2*(l_j2g2*m2 + l_j2j3*m3) * ( (U*cos(psi1)-yd*sin(psi1)-U) * cos(psi2) + (U*sin(psi1)+yd*cos(psi1)-yd) * sin(psi2) );
c55 = -1/2*l_j3g3*m3 * ( (U*cos(psi1)-yd*sin(psi1)-U) * cos(psi3) + (U*sin(psi1)+yd*cos(psi1)-yd) * sin(psi3) );
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

L_lat = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;];

vxj = sym(zeros(3, 1));
vyj = sym(zeros(3, 1));
vxj(1) = U;
vyj(1) = yd;
vxj(2) = U + l_j1j2*sin(psi1)*psi1d;
vyj(2) = yd - l_j1j2*cos(psi1)*psi1d;
vxj(3) = U + l_j1j2*sin(psi1)*psi1d + l_j2j3*sin(psi2)*psi2d;
vyj(3) = yd - l_j1j2*cos(psi1)*psi1d - l_j2j3*cos(psi2)*psi2d;
va = sym(zeros(2, 8));
psidj = [psi1d, psi2d, psi3d]';
for i = 1:8
    va(:, i) = [vxj(IU(i)); vyj(IU(i))] - n(:, IU(i))*l_ja(i)*psidj(IU(i));
end

va_v = sym(zeros(2, 8));
for i = 1:8
    va_v(:, i) = inv(Rot_w2v(:, :, IU(i)))*va(:, i);
end

beta_wv_y = sym(zeros(3, 1));
beta_wv_y(1) = (-(U*cos(psi1) - yd*sin(psi1))*sin(psi1) + (U*sin(psi1) + yd*cos(psi1))*cos(psi1)) - (-U * sin(psi1) + yd * cos(psi1));
beta_wv_y(2) = (-(U*cos(psi1) - yd*sin(psi1))*sin(psi2) + (U*sin(psi1) + yd*cos(psi1))*cos(psi2)) - (-U * sin(psi2) + yd * cos(psi2));
beta_wv_y(3) = (-(U*cos(psi1) - yd*sin(psi1))*sin(psi3) + (U*sin(psi1) + yd*cos(psi1))*cos(psi3)) - (-U * sin(psi3) + yd * cos(psi3));
beta_wv_x = sym(zeros(3, 1));
beta_wv_x(1) = ((U*cos(psi1) - yd*sin(psi1))*cos(psi1) + (U*sin(psi1) + yd*cos(psi1))*sin(psi1)) - (U * cos(psi1) + yd * sin(psi1));
beta_wv_x(2) = ((U*cos(psi1) - yd*sin(psi1))*cos(psi2) + (U*sin(psi1) + yd*cos(psi1))*sin(psi2)) - (U * cos(psi2) + yd * sin(psi2));
beta_wv_x(3) = ((U*cos(psi1) - yd*sin(psi1))*cos(psi3) + (U*sin(psi1) + yd*cos(psi1))*sin(psi3)) - (U * cos(psi3) + yd * sin(psi3));

beta = sym(zeros(8, 1));
for i = 1:8
    beta(i) = atan2(va_v(2, i) + beta_wv_y(IU(i)), va_v(1, i) + beta_wv_x(IU(i)));
end

ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8]';
Ca = diag(ca);

% Tau_tire = L_long1 * (- Delta_s * Ca * (delta - beta1)) ...
%          + L_long2 * (- Delta_s * Ca * (delta - beta2)) ...
%          + L_lat1 * (Delta_c * Ca * (delta - beta1)) ...
%          + L_lat2 * (Delta_c * Ca * (delta - beta2));
Tau_tire = 2 * L_lat * Ca * (delta - beta);

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

Mv = simplify(Twv' * M * Twv);
Cv = simplify(Twv' * C * Twv + Twv' * M * Tdwv);
Dv = simplify(Twv' * D * Twv);
L_lat_V = simplify(Twv' * L_lat);
% disp(simplify(beta));
Qv = Twv' * Tau_tire;

% constant vehicle forward speed
qd_u = [yd; psi1d; psi2d; psi3d;];

% Mu * qdd_u + Cu * qd_u + Du * qd_u = 2 * L_lat_u * Ca * (delta - beta)
Mu = Mv(2:5, 2:5);
Cu1 = Cv(2:5, 2:5);
Cu2 = [0, (m1 + m2 + m3) * U, 0, 0;
       0, -(l_j1g1 * m1 + l_j1j2 * m2 + l_j1j2 * m3) * U / 2, 0, 0;
       0, -(l_j2g2 * m2 + l_j2j3 * m3) * cos(psi1 - psi2) * U, (l_j2g2 * m2 + l_j2j3 * m3) * cos(psi1 - psi2) * U / 2, 0;
       0, -l_j3g3 * m3 * cos(psi1 - psi3) * U, 0, l_j3g3 * m3 * cos(psi1 - psi3) * U / 2];
Cu = Cu1 + Cu2;
Du = Dv(2:5, 2:5);
L_lat_u = L_lat_V(2:5, :);
% disp(simplify(beta));

inv_Mu = simplify(inv(Mu))
