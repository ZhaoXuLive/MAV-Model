clear

% m_1 = m_3 = 11142, m_2 = 7651 
syms m1 m2 m3 real
% l_j1_j2 = 11.4, l_j2_j3 = 9.4, l_j1_g1 = 11.4 - 5.538, l_j2_g2 = 9.4 - 4.757, l_j3_g3 = 5.538 
syms l_j1j2 l_j2j3 l_j1g1 l_j2g2 l_j3g3 real
% l_j1_a1 = 2.5, l_j1_a2 = 4, l_j1_a3 = 9.9, l_j2_a4 = 1.5, l_j2_a5 = 7.9, l_j3_a6 = 1.5, l_j3_a7 = 7.4, l_j3_a8 = 8.9
syms l_j1a1 l_j1a2 l_j1a3 l_j2a4 l_j2a5 l_j3a6 l_j3a7 l_j3a8 real
% h_11 = h_12 = 0.9625, h_13 = 0.98, h_21 = h_22 = 0.98, h_31 = 0.98, h_32 = h_33 = 0.9625
syms h1 h2 h3 real
% ca1 = ca2 = ca3 = ca4 = ca5 = ca6 = ca7 = ca8 = 400e3
syms ca_1 ca_2 ca_3 ca_4 ca_5 ca_6 ca_7 ca_8 real
% Izz1 = Izz3 = 143570, Izz2 = 37800
syms Iz1 Iz2 Iz3 real
% d_1 = d_2 = 0
syms d1 d2 real
% r_1 = r_2 = r_3 = r_4 = r_5 = r_6 = r_7 = r_8 = 0.4918
syms r1 r2 r3 r4 r5 r6 r7 r8 real
% fixed parameter
syms Tdr1 Tdr2 Tdr3 Tdr4 Tdr5 Tdr6 Tdr7 Tdr8 real
% input
syms delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 real
% state variable
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [xd; yd; psi1d; psi2d; psi3d;];

M = [m1 + m2 + m3, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1), (l_j2g2*m2 + l_j2j3*m3)*sin(psi2), l_j3g3*m3*sin(psi3);
     0, m1 + m2 + m3, -(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1), -(l_j2g2*m2 + l_j2j3*m3)*cos(psi2), -l_j3g3*m3*cos(psi3);
     (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1), -(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1), Iz1+m1*(l_j1g1^2)+(m2+m3)*(l_j1j2^2), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*cos(psi1-psi2), l_j3g3*l_j1j2*m3*cos(psi1-psi3);
     (l_j2g2*m2 + l_j2j3*m3)*sin(psi2), -(l_j2g2*m2 + l_j2j3*m3)*cos(psi2), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*cos(psi1-psi2), Iz2+m2*(l_j2g2^2)+m3*(l_j2j3^2), l_j3g3*l_j2j3*m3*cos(psi2-psi3);
     l_j3g3*m3*sin(psi3), -l_j3g3*m3*cos(psi3), l_j3g3*l_j1j2*m3*cos(psi1-psi3), l_j3g3*l_j2j3*m3*cos(psi2-psi3), Iz3+m3*(l_j3g3^2);]

C = [0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, l_j3g3*m3*cos(psi3)*psi3d;
     0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, l_j3g3*m3*sin(psi3)*psi3d;
     1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, 1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*(xd*cos(psi1)+yd*sin(psi1))+1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi2d+1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi3d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, 1/2*(l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(l_j2g2*m2+l_j2j3*m3)*(xd*cos(psi2)+yd*sin(psi2))-1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi1d+1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi3d, l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*l_j3g3*m3*cos(psi3)*psi3d, 1/2*l_j3g3*m3*sin(psi3)*psi3d, l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-psi1d+1/2*psi3d), l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*l_j3g3*m3*(xd*cos(psi3)+yd*sin(psi3))-1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi1d-1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi2d;]

D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, d1, -d1, 0;
     0, 0, -d1, d1+d2, -d2;
     0, 0, 0, -d2, d2;]

L_long1 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           -h1, -h1, -h1, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, -h2, -h2, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, -h3, -h3, -h3;]

L_long2 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           h1, h1, h1, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, h2, h2, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, h3, h3, h3;]

L_lat1 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;]

L_lat2 = [-sin(psi1), -sin(psi1), -sin(psi1), -sin(psi2), -sin(psi2), -sin(psi3), -sin(psi3), -sin(psi3);
          cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi2), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3), -l_j1j2*cos(psi1-psi3);
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3), -l_j2j3*cos(psi2-psi3);
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;]

R = [r1, r2, r3, r4, r5, r6, r7, r8];
Ru = diag(R)

Tdr = [Tdr1; Tdr2; Tdr3; Tdr4; Tdr5; Tdr6; Tdr7; Tdr8;]
Tdr_D = diag(Tdr)

ca = [ca_1; ca_2; ca_3; ca_4; ca_5; ca_6; ca_7; ca_8;];
Ca = diag(ca)

delta = [delta1; delta2; delta3; delta4; delta5; delta6; delta7; delta8;]
Delta = diag(delta)

beta11 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a1*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h1*psi1d)); 
beta12 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a1*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h1*psi1d)); 
beta21 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a2*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h1*psi1d)); 
beta22 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a2*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h1*psi1d)); 
beta31 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a3*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h1*psi1d)); 
beta32 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a3*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h1*psi1d)); 

beta41 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a4*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d - h2*psi2d)); 
beta42 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a4*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d + h2*psi2d));
beta51 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d - h2*psi2d)); 
beta52 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d + h2*psi2d));

beta61 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a6*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h3*psi3d)); 
beta62 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a6*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h3*psi3d));
beta71 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a7*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h3*psi3d)); 
beta72 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a7*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h3*psi3d));
beta81 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a8*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h3*psi3d)); 
beta82 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a8*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h3*psi3d));

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

qdd = inv(M)*(L_long1*(inv(Ru)*Tdr - Delta*Ca*(delta-beta1)) + L_lat1*(inv(Ru)*Tdr_D*delta + Ca*(delta-beta1)) + L_long2*(inv(Ru)*Tdr - Delta*Ca*(delta-beta2)) + L_lat2*(inv(Ru)*Tdr_D*delta + Ca*(delta-beta2)) - C*qd + D*qd);

A1 = jacobian(qdd, q);
B1 = jacobian(qdd, qd);

% A1l = subs(A1,{x, y, psi1, psi2, psi3, xd, yd, psi1d, psi2d, psi3d},{10, 5, 0, 0, 0, 1, 1, 0, 0, 0})


