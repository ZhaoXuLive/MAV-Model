function qdd_output = newton_model(unit_n, axle_m, ..., 
    getNumUnit, m, Iz, Ca, lu, lg, la, q_input, qd_input, delta_input, sat)

% calculate U \beta r \theta
psi = q_input(3:unit_n+2);
dpsi = qd_input(3:unit_n+2);
delta = delta_input';

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

% judge the slip angle 

v_aw = zeros(2, axle_m);
v_av = zeros(2, axle_m);
beta_a = zeros(axle_m, 1);
for i = 1:axle_m
    v_aw(1, i) = v_jw(1, getNumUnit(i)) + sin(psi(getNumUnit(i))) * la(getNumUnit(i)) * dpsi(getNumUnit(i));
    v_aw(2, i) = v_jw(2, getNumUnit(i)) - cos(psi(getNumUnit(i))) * la(getNumUnit(i)) * dpsi(getNumUnit(i));
    v_av(1, i) = v_aw(1, i) * cos(psi(getNumUnit(i))) + v_aw(2, i) * sin(psi(getNumUnit(i)));
    v_av(2, i) = -v_aw(1, i) * sin(psi(getNumUnit(i))) + v_aw(2, i) * cos(psi(getNumUnit(i)));
    beta_a(i) = v_av(2, i) / U;
end

alpha = delta - beta_a;
for i = 1:8
    if abs(alpha(i)) > sat * pi / 180
        alpha(i) = sign(alpha(i)) * sat * pi / 180;
    end
end

yij = Ca * alpha;

% newton 
L = cell(1, unit_n);
for k = 1:unit_n
    % L1: 各轴中心到单元质心的距离，前为正，后为负
    L1 = lg(k) - la(getNumUnit == k);
    L{k} = L1;
end

nx = 3*unit_n-1; % number of states
nv = 4*unit_n-2; % number of variables (including joint forces)

% S * (\dot{x} R)' + P x = By y 方程
% x = {\beta r \theta}  n n n-1 状态量
% r = \dot{\psi} 横摆角速度
% \theta = \psi_j - \psi_j+1  铰接角

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
% 3n-1个状态，前n个是质心侧偏角，中间n个是横摆角速度，最后n-1个是铰接角
b1 = unit_n;         % last row of beta_dot
b2 = 2 * unit_n;     % last row of r_dot
b3 = 3 * unit_n - 1; % last row of theta_dot

% also need number of 'preceding' axles at each unit
% 每个单元前面所有的单元共有多少个轴
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

% AA = S \ (BY * K1 - P);
% A = AA(1:nx, :);
% 
% BB = S \ BY * K2;
% B = BB(1:nx, :);
% 
% % /dot x = AX +B /delta
% x = [beta; r; theta];
% dx = A * x + B * delta;
% dx= dx(1:nx, :);

x = [beta; r; theta];
dxx = S \ (BY * yij - P * x);
dx= dxx(1:nx, :);

% dx to new qdd
dbeta_1 = dx(1);
qdd_v = zeros(2, 1);
qdd_v(1) = 0;
qdd_v(2) = dbeta_1 * U + lg(1) * dx(unit_n+1) + dx(unit_n+1) * U;

qd_v = zeros(2, 1);
qd_v(1) = qd_input(1) * cos(psi(1)) + qd_input(2) * sin(psi(1));
qd_v(2) = -qd_input(1) * sin(psi(1)) + qd_input(2) * cos(psi(1));
Twv = [cos(psi(1)), -sin(psi(1));
     sin(psi(1)), cos(psi(1))];
Tdwv = [-sin(psi(1))*dpsi(1), -cos(psi(1))*dpsi(1);
        cos(psi(1))*dpsi(1), -sin(psi(1))*dpsi(1)];

qdd_output = zeros(unit_n+2, 1);
qdd_output(1:2) = Twv * qdd_v + Tdwv * qd_v;
qdd_output(3:unit_n+2) = dx(unit_n+1:2*unit_n);

end

% BB = S \ BY * K2;
% B1 = BB(1:nx, :);
% BB = S \ BM;
% B2 = BB(1:nx, :)
% 
% Bj = [100000, 100000, 100000]'; %joint damping n dim
% 
% KM = zeros(unit_n, 2*unit_n);
% %structure of the individual damping sub-matrix
% S = [-1, 1; 1, -1]; 
% for i = 1:unit_n-1
%     KM(i:i+1, i+1:i+2) = KM(i:i+1, i+2:i+2) + Bj(i)*S;
% end
% 
% A=A+BM*KM; %close the loop via joint damping