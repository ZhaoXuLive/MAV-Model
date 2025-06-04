clear

% m_1 = m_3 = 11142, m_2 = 7651 
% l_j1_j2 = 11.4, l_j2_j3 = 9.4, l_j1_g1 = 5.862, l_j2_g2 = 4.643, l_j3_g3 = 5.538 
% l_j1_a1 = 2.5, l_j1_a2 = 4, l_j1_a3 = 9.9, l_j2_a4 = 1.5, l_j2_a5 = 7.9, l_j3_a6 = 1.5, l_j3_a7 = 7.4, l_j3_a8 = 8.9
% h_11 = h_12 = 0.9625, h_13 = 0.98, h_21 = h_22 = 0.98, h_31 = 0.98, h_32 = h_33 = 0.9625
% ca1 = ca2 = ca3 = ca4 = ca5 = ca6 = ca7 = ca8 = 400e3
% Izz1 = Izz3 = 143570, Izz2 = 37800
% d_1 = d_2 = 0

% fixed parameter
% input
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real
% state variable
syms y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real
syms vx1 real

qd = [xd; yd; psi1d; psi2d; psi3d];

M = [11142 + 7651 + 11142, 0, (5.862*11142 + 11.4*(7651 + 11142))*sin(psi1), (4.643*7651 + 9.4*11142)*sin(psi2), 5.538*11142*sin(psi3);
     0, 11142 + 7651 + 11142, -(5.862*11142 + 11.4*(7651 + 11142))*cos(psi1), -(4.643*7651 + 9.4*11142)*cos(psi2), -5.538*11142*cos(psi3);
     (5.862*11142 + 11.4*(7651 + 11142))*sin(psi1), -(5.862*11142 + 11.4*(7651 + 11142))*cos(psi1), 143570+11142*(5.862^2)+(7651+11142)*(11.4^2), (4.643*11.4*7651+11.4*9.4*11142)*cos(psi1-psi2), 5.538*11.4*11142*cos(psi1-psi3);
     (4.643*7651 + 9.4*11142)*sin(psi2), -(4.643*7651 + 9.4*11142)*cos(psi2), (4.643*11.4*7651+11.4*9.4*11142)*cos(psi1-psi2), 37800+7651*(4.643^2)+11142*(9.4^2), 5.538*9.4*11142*cos(psi2-psi3);
     5.538*11142*sin(psi3), -5.538*11142*cos(psi3), 5.538*11.4*11142*cos(psi1-psi3), 5.538*9.4*11142*cos(psi2-psi3), 143570+11142*(5.538^2);]

C = [0, 0, (5.862*11142 + 11.4*(7651 + 11142))*cos(psi1)*psi1d, (4.643*7651 + 9.4*11142)*cos(psi2)*psi2d, 5.538*11142*cos(psi3)*psi3d;
     0, 0, (5.862*11142 + 11.4*(7651 + 11142))*sin(psi1)*psi1d, (4.643*7651 + 9.4*11142)*sin(psi2)*psi2d, 5.538*11142*sin(psi3)*psi3d;
     1/2*(5.862*11142 + 11.4*(7651 + 11142))*cos(psi1)*psi1d, 1/2*(5.862*11142 + 11.4*(7651 + 11142))*sin(psi1)*psi1d, -1/2*(5.862*11142 + 11.4*(7651 + 11142))*(xd*cos(psi1)+yd*sin(psi1))+1/2*(4.643*11.4*7651+11.4*9.4*11142)*sin(psi1-psi2)*psi2d+1/2*5.538*11.4*11142*sin(psi1-psi3)*psi3d, (4.643*11.4*7651+11.4*9.4*11142)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), 5.538*11.4*11142*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(4.643*7651 + 9.4*11142)*cos(psi2)*psi2d, 1/2*(4.643*7651 + 9.4*11142)*sin(psi2)*psi2d, (4.643*11.4*7651+11.4*9.4*11142)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(4.643*7651+9.4*11142)*(xd*cos(psi2)+yd*sin(psi2))-1/2*(4.643*11.4*7651+11.4*9.4*11142)*sin(psi1-psi2)*psi1d+1/2*5.538*9.4*11142*sin(psi2-psi3)*psi3d, 5.538*9.4*11142*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*5.538*11142*cos(psi3)*psi3d, 1/2*5.538*11142*sin(psi3)*psi3d, 5.538*11.4*11142*sin(psi1-psi3)*(-psi1d+1/2*psi3d), 5.538*9.4*11142*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*5.538*11142*(xd*cos(psi3)+yd*sin(psi3))-1/2*5.538*11.4*11142*sin(psi1-psi3)*psi1d-1/2*5.538*9.4*11142*sin(psi2-psi3)*psi2d;]

D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;];

L_long1 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           -0.9625, -0.9625, -0.98, 11.4*sin(psi1-psi2), 11.4*sin(psi1-psi2), 11.4*sin(psi1-psi3), 11.4*sin(psi1-psi3), 11.4*sin(psi1-psi3);
           0, 0, 0, -0.98, -0.98, 9.4*sin(psi2-psi3), 9.4*sin(psi2-psi3), 9.4*sin(psi2-psi3);
           0, 0, 0, 0, 0, -0.98, -0.9625, -0.9625;]

L_long2 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           0.9625, 0.9625, 0.98, 11.4*sin(psi1-psi2), 11.4*sin(psi1-psi2), 11.4*sin(psi1-psi3), 11.4*sin(psi1-psi3), 11.4*sin(psi1-psi3);
           0, 0, 0, 0.98, 0.98, 9.4*sin(psi2-psi3), 9.4*sin(psi2-psi3), 9.4*sin(psi2-psi3);
           0, 0, 0, 0, 0, 0.98, 0.9625, 0.9625;]

L_lat1 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -2.5, -4, -9.9, -11.4*cos(psi1-psi2), -11.4*cos(psi1-psi2), -11.4*cos(psi1-psi3), -11.4*cos(psi1-psi3), -11.4*cos(psi1-psi3);
          0, 0, 0, -1.5, -7.9, -9.4*cos(psi2-psi3), -9.4*cos(psi2-psi3), -9.4*cos(psi2-psi3);
          0, 0, 0, 0, 0, -1.5, -7.4, -8.9;]

L_lat2 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -2.5, -4, -9.9, -11.4*cos(psi1-psi2), -11.4*cos(psi1-psi2), -11.4*cos(psi1-psi3), -11.4*cos(psi1-psi3), -11.4*cos(psi1-psi3);
          0, 0, 0, -1.5, -7.9, -9.4*cos(psi2-psi3), -9.4*cos(psi2-psi3), -9.4*cos(psi2-psi3);
          0, 0, 0, 0, 0, -1.5, -7.4, -8.9;]

ca = [400e3; 400e3; 400e3; 400e3; 400e3; 400e3; 400e3; 400e3;];
Ca = diag(ca)

delta = [delta1; delta2; delta3; delta4; delta5; delta6; delta7; delta8;];
delta_s = [sin(delta1); sin(delta2); sin(delta3); sin(delta4); sin(delta5); sin(delta6); sin(delta7); sin(delta8);];
delta_c = [cos(delta1); cos(delta2); cos(delta3); cos(delta4); cos(delta5); cos(delta6); cos(delta7); cos(delta8);];
Delta_s = diag(delta_s);
Delta_c = diag(delta_c);

beta11 = atan((-xd*sin(psi1) + yd*cos(psi1) - 2.5*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - 0.9625*psi1d)); 
beta12 = atan((-xd*sin(psi1) + yd*cos(psi1) - 2.5*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + 0.9625*psi1d)); 
beta21 = atan((-xd*sin(psi1) + yd*cos(psi1) - 4*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - 0.9625*psi1d)); 
beta22 = atan((-xd*sin(psi1) + yd*cos(psi1) - 4*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + 0.9625*psi1d)); 
beta31 = atan((-xd*sin(psi1) + yd*cos(psi1) - 9.9*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - 0.98*psi1d)); 
beta32 = atan((-xd*sin(psi1) + yd*cos(psi1) - 9.9*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + 0.98*psi1d)); 

beta41 = atan((-xd*sin(psi2) + yd*cos(psi2) - 11.4*cos(psi1-psi2)*psi1d - 1.5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + 11.4*sin(psi1-psi2)*psi1d - 0.98*psi2d)); 
beta42 = atan((-xd*sin(psi2) + yd*cos(psi2) - 11.4*cos(psi1-psi2)*psi1d - 1.5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + 11.4*sin(psi1-psi2)*psi1d + 0.98*psi2d));
beta51 = atan((-xd*sin(psi2) + yd*cos(psi2) - 11.4*cos(psi1-psi2)*psi1d - 7.9*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + 11.4*sin(psi1-psi2)*psi1d - 0.98*psi2d)); 
beta52 = atan((-xd*sin(psi2) + yd*cos(psi2) - 11.4*cos(psi1-psi2)*psi1d - 7.9*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + 11.4*sin(psi1-psi2)*psi1d + 0.98*psi2d));

beta61 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 1.5*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d - 0.98*psi3d)); 
beta62 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 1.5*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d + 0.98*psi3d));
beta71 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 7.4*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d - 0.9625*psi3d)); 
beta72 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 7.4*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d + 0.9625*psi3d));
beta81 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 8.9*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d - 0.9625*psi3d)); 
beta82 = atan((-xd*sin(psi3) + yd*cos(psi3) - 11.4*cos(psi1-psi3)*psi1d - 9.4*cos(psi2-psi3)*psi2d - 8.9*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + 11.4*sin(psi1-psi3)*psi1d + 9.4*sin(psi2-psi3)*psi2d + 0.9625*psi3d));

beta1 = [beta11; beta21; beta31; beta41; beta51; beta61; beta71; beta81]
beta2 = [beta12; beta22; beta32; beta42; beta52; beta62; beta72; beta82]

% main equation
% M*qdd + C*qd + D*qd = L_long1*((Ru^-1)*Tdr - Delta*Ca*(Delta-Beta1)) + L_lat1*((Ru^-1)*Tdr_D*delta + Ca*(delta-beta1)) + L_long2*((Ru^-1)*Tdr - Delta*Ca*(Delta-Beta2)) + L_lat2*((Ru^-1)*Tdr_D*delta + Ca*(delta-beta2));

Twv = [cos(psi1), -sin(psi1), 0, 0, 0;
     sin(psi1), cos(psi1), 0, 0, 0;
     0, 0, 1, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1;]
% need to check
Tdwv = diff(Twv) * psi1d

Mv = Twv'*M*Twv
Cv = Twv'*C*Twv + Twv'*M*Tdwv
Dv = Twv'*D
L_long1_v = Twv'*L_long1
L_long2_v = Twv'*L_long2
L_lat1_v = Twv'*L_lat1
L_lat2_v = Twv'*L_lat2

qdd = inv(Mv)*(-L_long1_v*Delta_s*Ca*(delta-beta1) + L_lat1_v*Delta_c*Ca*(delta-beta1) + L_long2_v*(-Delta_s*Ca*(delta-beta2)) + L_lat2_v*Delta_c*Ca*(delta-beta2) - C_v*qd + D_v*qd);

xqdd = subs(qdd, {xd, yd, psi1d, psi2d, psi3d}, {vx1, vx1*yd, vx1*psi1d, vx1*psi2d, vx1*psi3d})

xi = [y; psi1; psi2; psi3; yd; psi1d; psi2d; psi3d;];
xid = [1/vx1*yd;
       1/vx1*psi1d;
       1/vx1*psi2d;
       1/vx1*psi3d;
       1/(vx1^2)*xqdd(2) + 1/vx1*yd;
       1/(vx1^2)*xqdd(3) + 1/vx1*psi1d;
       1/(vx1^2)*xqdd(4) + 1/vx1*psi2d;
       1/(vx1^2)*xqdd(5) + 1/vx1*psi3d;];

A1t = jacobian(xid, xi);
B1t = jacobian(xid, delta);

% t0 = cputime;
tic;
A2lt = subs(A1t,{psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d, ...
    delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8} ...
    , {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
A2lt = vpa(A2lt, 3)
B2lt = subs(B1t,{psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d, ...
    delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8} ...
    , {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
B2lt = vpa(B2lt, 3)
% elapsed_time = cputime - t0
elapsed_time = toc

