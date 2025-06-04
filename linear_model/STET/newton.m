% articulated vehicle linear model
% Linear Model of n units with N axles
% n and N determined from vehicle data
% states: body slip angle for unit 1, yaw rates, articulation angles (2n)
%   = [beta_1,r1,r2,...,rn,theta_1,theta_2,...theta_n-1]
% inputs: steer angles at each axle (N inputs) .. input matrix B1
%     AND yaw moments at each unit (n inputs)  .. input matrix B2
% 
% outputs: all states
%
% M= n-vector of masses (row or col, n elements, mass in kg)
% J= n-vector of CG yaw moments of inertia (ditto, kg m^2)
% L= n-element cell array of axle locations {L1,L2,...}
%   Li=[x1,x2,...] x coords relative to the mass centre (m) on each unit
%  (so the units may have different numbers of axles)
% Ca=cell array of axle cornering stiffnesses - same structure as L
%   units: N/rad
% Xf= n-vector of front joint x coords rel to CG (front of unit 1 is not used)
% Xr= n-vector of rear joints (similar)
% U=vehicle speed (m/s) ... assumed constant

% PARSu: [M J L xg Bj]
% PARSa: [xa u Ca DRIVE BRAKE STEER AUTOSTEER]
% xg: 各单元质心到后面的铰接点的距离 均为正 xr = -xg 均为负
% xa: 各轴中心到后面的铰接点的距离  
% xf: 各单元质心到前面铰接点的距离 均为正 xf = L - xg

n=3; % number of units
N=8; % number of axles

syms m1 m2 m3 real
syms Iz1 Iz3 Iz2 real
syms l_j1j2 l_j2j3 l_j1g1 l_j2g2 l_j3g3 l_j1a1 l_j1a2 l_j1a3 l_j2a4 l_j2a5 l_j3a6 l_j3a7 l_j3a8 real
syms ca1 ca2 ca3 ca4 ca5 ca6 ca7 ca8 real
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real
syms U real

M = [m1, m2, m3];
J = [Iz1, Iz2, Iz3];

l_g1j2 = l_j1j2 - l_j1g1;
l_g2j3 = l_j2j3 - l_j2g2;
l_g3j4 = l_j1j2 - l_j3g3;
Xr = -[l_g1j2, l_g2j3, l_g3j4];

Xf = [l_j1g1, l_j2g2, l_j3g3];

L = cell(1,n);
Ca = cell(1,n);
ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8];

l_a1j2 = l_j1j2 - l_j1a1;
l_a2j2 = l_j1j2 - l_j1a2;
l_a3j2 = l_j1j2 - l_j1a3;
l_a4j3 = l_j2j3 - l_j2a4;
l_a5j3 = l_j2j3 - l_j2a5;
l_a6j4 = l_j1j2 - l_j3a6;
l_a7j4 = l_j1j2 - l_j3a7;
l_a8j4 = l_j1j2 - l_j3a8;
Xa = [l_a1j2, l_a2j2, l_a3j2, l_a4j3, l_a5j3, l_a6j4, l_a7j4, l_a8j4];

axle4unit = [1, 1, 1, 2, 2, 3, 3, 3];
for k = 1:n
    % ind: 第k单元对应的轴，10*1 logical
    ind = (axle4unit == k);
    % L1: 各轴中心到单元质心的距离，前为正
    L1 = Xa(ind) + Xr(k);
    % 对应轴侧偏刚度
    Ca1 = ca(ind);
    L{k} = L1;
    Ca{k} = Ca1;
end

nx = 3*n-1; % number of states
nv = 4*n-2; % number of variables (including joint forces)

% dimensions of primary matrices
S = sym(zeros(nv, nv));
P = sym(zeros(nv, nx));
BY = sym(zeros(nv, N));
K1 = sym(zeros(N, nx));
K2 = sym(zeros(N, N));

% input matrix for external yaw moments
BM = zeros(nv,n); 

% populate S,P,B equation-by-equation

% first define end-index of sub-matrix blocks
% 3n-1个状态，前n个是质心侧偏角，中间n个是横摆角速度，最后n-1个是铰接角
b1 = n;         % last row of beta_dot
b2 = 2 * n;     % last row of r_dot
b3 = 3 * n - 1; % last row of theta_dot

% also need number of 'preceding' axles at each unit
% 每个单元前面所有的单元共有多少个轴
pax = sym(zeros(1,n)); % preceding axle numbers
for i = 2:n
    pax(i) = pax(i - 1) + length(L{i - 1});
end

% equation (1)
for i = 1:n 
    S(i, i) = M(i) * U;
    P(i, b1 + i) = M(i) * U;
    if i > 1
        % R_(i - 1) term in equ(1)
        S(i, b3 + i - 1) = -1;
    end 
    if i < n
        % R_i term 
        S(i, b3 + i) = 1;
    end 
    for j = 1:length(L{i})
        % Yij terms
       BY(i, pax(i) + j) = 1; 
    end
end

% equation (2)
for i = 1:n 
    S(b1 + i, b1 + i) = J(i);
    if i > 1
        % R_(i-1) term in equ(1)
        S(b1 + i, b3 + i - 1) = -Xf(i);
    end 
    if i < n
        % R_i term 
        S(b1 + i, b3 + i) = Xr(i);
    end 
    for j = 1:length(L{i})
        % Mij terms
       BY(b1 + i, pax(i) + j) = L{i}(j); 
    end
    BM(b1+i) = 1; % direct yaw moment term
end

% equation (3)
for i = 1:n-1
    % theta-dot term
    S(b2 + i, b2 + i) = 1; 
    % r_i in notes = unit in front of the joint
    P(b2 + i, b1 + i) = -1; 
    % r_(i+1)
    P(b2 + i, b1 + i + 1) = 1; 
end

% equation (4) 
for i = 1:n-1 
    S(b3 + i, i + 1) = U;
    S(b3 + i, i) = -U; 
    S(b3 + i, b1 + i + 1) = Xf(i + 1);
    S(b3 + i, b1 + i) = -Xr(i);
    S(b3 + i, b2 + i) = -U; 
end

% K1 matrix, Y = K1 * X for each axle force
for i = 1:n
    for j = 1:length(L{i})
        % 'flat' axle index
        ia = pax(i) + j; 
        K1(ia, i) = U;
        K1(ia, b1 + i) = L{i}(j);
        K1(ia, :) = -K1(ia, :) * Ca{i}(j) / U;
    end
end

% K2 matrix, go axle by axle
for i = 1:n
    for j = 1:length(L{i})
        % 'flat' axle index
        ia = pax(i) + j; 
        K2(ia, ia) = Ca{i}(j);
    end
end

% disp(L);
% disp(pax);

AA = S \ (BY * K1 - P);
A = simplify(AA(1:nx, :));

disp(BM);

BB = S \ BY * K2;
B = simplify(BB(1:nx, :));

% n_state: beta1 beta2 beta3 psi1d psi2d psi3d theta1d theta2d
% l_state: yj1 psi1 psi2 psi3

% beta1 * U = v_yj1_d; U = v_xj1_d;
% transform
% 第一步，建立真正意义上的恒定速度、车辆坐标系下的动力学方程，state以l_state为例；
% 第二步，建立两种状态之间的关系，通过变换得到牛顿方程下的l_state方程，统一二者；

