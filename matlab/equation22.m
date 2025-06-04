clear

%%
% fixed parameter
syms m1 m2 m3 real
syms iz1 iz2 iz3 real
syms l_j1j2 l_j2j3 l_j1g1 l_j2g2 l_j3g3 real
syms l_j1a1 l_j1a2 l_j1a3 l_j2a4 l_j2a5 l_j3a6 l_j3a7 l_j3a8 real
syms h11 h12 h13 h21 h22 h31 h32 h33 real
syms ca1 ca2 ca3 ca4 ca5 ca6 ca7 ca8 real
syms d1 d2 real
syms r1 r2 r3 r4 r5 r6 r7 r8 real

% input
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real
syms Tdr11 Tdr12 Tdr13 Tdr14 Tdr15 Tdr16 Tdr17 Tdr18 real
syms Tdr21 Tdr22 Tdr23 Tdr24 Tdr25 Tdr26 Tdr27 Tdr28 real
% state variable
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [xd; yd; psi1d; psi2d; psi3d;];

%%
syms l_kj [3, 3]
l_kj = [l_j1g1, l_j1j2, l_j1j2;
             0, l_j2g2, l_j2j3;
             0,      0, l_j3g3;];
delta_kj = eye(3);

for j = 1:3
    Jc(:, :, j) = [1, 0,  l_kj(1, j) * sin(psi1),  l_kj(2, j) * sin(psi2),  l_kj(3, j) * sin(psi3);
          0, 1, -l_kj(1, j) * cos(psi1), -l_kj(2, j) * cos(psi2), -l_kj(3, j) * cos(psi3);
          0, 0,          delta_kj(1, j),          delta_kj(2, j),          delta_kj(3, j);];
end

mj = [m1, m2, m3]';
izj = [iz1, iz2, iz3]';
syms Mj [3, 3, 3]
for j = 1:3
    Mj(:, :, j) = diag([mj(j), mj(j), izj(j)]);
end

M = zeros(5, 5);
for j = 1:3
    M = M + Jc(:, :, j)' * Mj(:, :, j) * Jc(:, :, j);
end
disp(simplify(M));

C = sym(zeros(5));
for i = 1:5
    for l = 1:5
        for k = 1:5
            C(i, l) = C(i, l) + (diff(M(i, l), q(k)) - 1/2 * diff(M(k, l), q(i))) * qd(k);
        end
    end
end
disp(simplify(C*qd));    

dj = [0, d1, d2, 0]';
D = sym(zeros(5));
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
% disp(D);

psij = [psi1, psi2, psi3]';
syms Rot_w2v [2, 2, 3]
for i = 1:3
    Rot_w2v(:, :, i) = [cos(psij(i)), -sin(psij(i));
                        sin(psij(i)), cos(psij(i));];
end

deltaj = [delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8]';
syms Rot_v2t [2, 2, 8]
for i = 1:8
    Rot_v2t(:, :, i) = [cos(deltaj(i)), -sin(deltaj(i));
                        sin(deltaj(i)), cos(deltaj(i));];
end

xj(1) = x;
yj(1) = y;
xj(2) = x - l_j1j2*cos(psi1);
yj(2) = y - l_j1j2*sin(psi1);
xj(3) = x - l_j1j2*cos(psi1) - l_j2j3*cos(psi2);
yj(3) = y - l_j1j2*sin(psi1) - l_j2j3*sin(psi2);

for i = 1:3
    u(:, i)=[cos(q(i + 2)), sin(q(i + 2))]';
    n(:, i)=[-sin(q(i + 2)), cos(q(i + 2))]';
end

IU = [1, 1, 1, 2, 2, 3, 3, 3]';
l_ja = [l_j1a1, l_j1a2, l_j1a3, l_j2a4, l_j2a5, l_j3a6, l_j3a7, l_j3a8]';
% ha = [h11, h12, h13, h21, h22, h31, h32, h33]';
ha = zeros(8, 1);
r_a1 = sym(zeros(2, 8));
r_a2 = sym(zeros(2, 8));
for i = 1:8
    r_a1(:, i) = [xj(IU(i)); yj(IU(i))] - u(:, IU(i))*l_ja(i) + n(:, IU(i))*ha(i);
    r_a2(:, i) = [xj(IU(i)); yj(IU(i))] - u(:, IU(i))*l_ja(i) + (-1)*n(:, IU(i))*ha(i);
end
% disp(r_a1);

L_long1 = sym(zeros(5, 8));
L_long2 = sym(zeros(5, 8));
L_lat1 = sym(zeros(5, 8));
L_lat2 = sym(zeros(5, 8));

% for i = 1:8
%     for j = 1:5
%         L_long1(j, i) = (diff(r_a1(:, i), q(j)))' * Rot_w2v(:, 1, IU(i));
%         L_long2(j, i) = (diff(r_a2(:, i), q(j)))' * Rot_w2v(:, 1, IU(i));
%         L_lat1(j, i) = (diff(r_a1(:, i), q(j)))' * Rot_w2v(:, 2, IU(i));
%         L_lat2(j, i) = (diff(r_a2(:, i), q(j)))' * Rot_w2v(:, 2, IU(i));
%     end
% end
% disp(simplify(L_long1));
% disp(simplify(L_lat1));
for i = 1:8
    L_long1(:, i) = (jacobian(r_a1(:, i), q))' * Rot_w2v(:, 1, IU(i));
    L_long2(:, i) = (jacobian(r_a2(:, i), q))' * Rot_w2v(:, 1, IU(i));
    L_lat1(:, i) = (jacobian(r_a1(:, i), q))' * Rot_w2v(:, 2, IU(i));
    L_lat2(:, i) = (jacobian(r_a2(:, i), q))' * Rot_w2v(:, 2, IU(i));
end
% disp(simplify(L_long1));
% disp(simplify(L_lat1));

Delta_c = sym(zeros(8, 8));
Delta_s = sym(zeros(8, 8));
for i = 1:8
    Delta_c(i,i) = cos(deltaj(i));
    Delta_s(i,i) = sin(deltaj(i));
end

Tdrj1 = [Tdr11, Tdr12, Tdr13, Tdr14, Tdr15, Tdr16, Tdr17, Tdr18]';
Tdrj2 = [Tdr21, Tdr22, Tdr23, Tdr24, Tdr25, Tdr26, Tdr27, Tdr28]';
rj = [r1, r2, r3, r4, r5, r6, r7, r8]';
Rj = diag(rj);
% F_dr1 = inv(Rj)*Tdrj1

F_cor1 = sym(zeros(8, 1));
F_cor2 = sym(zeros(8, 1));

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
disp(simplify(va1_v));
disp(simplify(va2_v));

% beta1 = sym(zeros(8, 1));
% beta2 = sym(zeros(8, 1));
% for i = 1:8
%     beta1(i) = atan2(va1_v(2, i), va1_v(1, i));
%     beta2(i) = atan2(va2_v(2, i), va2_v(1, i));
% end
% 
% ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8]';
% Ca = diag(ca);
% % tau_tire = L_long1 * (Delta_c * inv(Rj) * Tdrj1 - Delta_s * Ca * (deltaj - beta1)) ...
% %          + L_long2 * (Delta_c * inv(Rj) * Tdrj2 - Delta_s * Ca * (deltaj - beta2)) ...
% %          + L_lat1 * (Delta_s * inv(Rj) * Tdrj1 + Delta_c * Ca * (deltaj - beta1)) ...
% %          + L_lat2 * (Delta_s * inv(Rj) * Tdrj2 + Delta_c * Ca * (deltaj - beta2));
% 
syms Fx [8, 1]; 
syms Fy [8, 1]; 
tau_tire = L_long1 * (Delta_c * Fx - Delta_s * Fy) + L_lat1 * (Delta_s * Fx + Delta_c * Fy);
disp(simplify(tau_tire));
