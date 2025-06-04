clear

%%
% coefficient
% m_1 = m_3 = 11142, m_2 = 7651 
% l_j1_j2 = 11.4, l_j2_j3 = 9.4, l_j1_g1 = 5.862, l_j2_g2 = 4.643, l_j3_g3 = 5.538 
% l_j1_a1 = 2.5, l_j1_a2 = 4, l_j1_a3 = 9.9, l_j2_a4 = 1.5, l_j2_a5 = 7.9, l_j3_a6 = 1.5, l_j3_a7 = 7.4, l_j3_a8 = 8.9
% h_11 = h_12 = 0.9625, h_13 = 0.98, h_21 = h_22 = 0.98, h_31 = 0.98, h_32 = h_33 = 0.9625
% ca1 = ca2 = ca3 = ca4 = ca5 = ca6 = ca7 = ca8 = 400e3
% Izz1 = Izz3 = 143570, Izz2 = 37800
% d_1 = d_2 = 0
% r_1 = r_2 = r_3 = r_4 = r_5 = r_6 = r_7 = r_8 = 0.4918

% fixed parameter
syms m1 m2 m3 real
syms iz1 iz2 iz3 real
syms l_j1j2 l_j2j3 l_j1g1 l_j2g2 l_j3g3 real
syms l_j1a1 l_j1a2 l_j1a3 l_j2a4 l_j2a5 l_j3a6 l_j3a7 l_j3a8 real
syms h11 h12 h13 h21 h22 h31 h32 h33 real
syms Ca [8, 1] 
syms Rt [8, 1]
syms d1 d2 real

% input
syms delta [8, 1] 
syms Tdr [8, 1] 
% state variable
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [xd; yd; psi1d; psi2d; psi3d;];

%%
% Parameter given
% M J L xg = mass, moment of inertia, length, cg location measured from rear joint

MM = [m1, m2, m3]';
JJ = [iz1, iz2, iz3]';
LL = [l_j1j2, l_j2j3, l_j1j2]';
XX = [l_j1g1, l_j2g2, l_j3g3]';
XG = LL - XX;

% Rotation matrix

for i = 1:3
    uu(:, i)=[cos(q(i + 2)), sin(q(i + 2))]';
    nn(:, i)=[-sin(q(i + 2)), cos(q(i + 2))]';
end

% Matrix generate
% Specify matrix dimensions in advance

D1 = diag(kron(MM', [1 1]));
D2 = diag([0 0 JJ']);

P1 = kron(ones(3, 1), eye(2));
syms P2 [6, 3]
for i = 1:3
    for j = 1:3
        if i > j
            P2(2 * i - 1 : 2 * i, j) = -LL(j) * nn(:, j);
        elseif i == j
            P2(2 * i - 1 : 2 * i, j) = -XX(j) * nn(:, j);
        else
            P2(2 * i - 1 : 2 * i, j) = 0;
        end
    end
end
P = [P1, P2];

M = simplify(P' * D1 * P + D2);

R1 = zeros(6, 2);
syms R2 [6, 3]
for i = 1:3
    for j = 1:3
        if i > j
            R2(2 * i - 1 : 2 * i, j) = LL(j) * uu(:, j);
        elseif i==j
            R2(2 * i - 1 : 2 * i, j) = XX(j) * uu(:, j);
        else
            R2(2 * i - 1 : 2 * i, j) = 0;
        end
    end
end
R = [R1, R2];

dTdq = simplify((qd' * P' * D1 * R)' .* qd);

Pd1 = zeros(6, 2);
syms Pd2 [6, 3]
for i = 1:3
    for j = 1:3
        if i > j
            Pd2(2 * i - 1 : 2 * i, j) = LL(j) * uu(:, j) * qd(j + 2);
        elseif i==j
            Pd2(2 * i - 1 : 2 * i, j) = XX(j) * uu(:, j) * qd(j + 2);
        else
            Pd2(2 * i - 1 : 2 * i, j) = 0;
        end
    end
end
Pd = [Pd1, Pd2];
Md = simplify(Pd'*D1*P + P'*D1*Pd)

C = simplify(Md*qd - dTdq)

%%
JA = [l_j1a1, l_j1a2, l_j1a3, l_j2a4, l_j2a5, l_j3a6, l_j3a7, l_j3a8]';
XA = [l_j1j2, l_j1j2, l_j1j2, l_j2j3, l_j2j3, l_j1j2, l_j1j2, l_j1j2]' - JA;
IU = [1, 1, 1, 2, 2, 3, 3, 3]';

syms Va [2, 8]      % axle velocities in globals
syms Vj [2, 4]      % joint velocities in globals
Va_v = Va;          % axle velocities in locals

Vj(:, 1) = [xd, yd]';
for i = 1:3
    Vj(:, i+1) = Vj(:, i) - qd(i + 2) * LL(i) * nn(:, i);
end

for i = 1:8
    Va(:, i) = Vj(:, IU(i) + 1) + qd(IU(i) + 2) * XA(i) * nn(:, IU(i));
    Vax_v = dot(Va(:, i), uu(:, IU(i)));
    Vay_v = dot(Va(:, i), nn(:, IU(i)));
    Va_v(:, i) = [Vax_v; Vay_v];
end
% disp(simplify(Va_v));

deg = pi/180;
syms Fy [8, 1]; 
syms Fx [8, 1]; 

% for i = 1:8
%     vx = Va_v(1, i) * cos(delta(i)) + Va_v(2, i) * sin(delta(i));
%     vy = -Va_v(1, i) * sin(delta(i)) + Va_v(2, i) * cos(delta(i));
%     
%     alpha = atan2(vy, vx);
%     Fy(i) = Ca(i) * alpha;
% %     Fx(i) = Tdr(i) / Rt(i);
% end

syms J [2, 5]
Q = zeros(5,1);

for i = 1:8
    ua(:, i)=[cos(q(IU(i) + 2) + delta(i)), sin(q(IU(i) + 2) + delta(i))]';
    na(:, i)=[-sin(q(IU(i) + 2) + delta(i)), cos(q(IU(i) + 2) + delta(i))]';
end

for i = 1:8
    Fa = [Fx(i); Fy(i)];
    iu = IU(i); 
    u = ua(:, i);
    n = na(:, i); 

    for j = 1:5
        J(:, j) = [0; 0];
    end
    J(1, 1) = 1; 
    J(2, 2) = 1;

    for j = 1:iu
        if j < iu
            J(:, j + 2) = -LL(j) * nn(:, j);
        else
            J(:, j + 2) = -(LL(j) - XA(i)) * nn(:, j);
        end
    end
    Q = Q + J' * [u n] * Fa;
end
disp(simplify(Q));

