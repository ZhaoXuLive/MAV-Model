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
syms Tdr1 Tdr2 Tdr3 Tdr4 Tdr5 Tdr6 Tdr7 Tdr8 real
% state variable
syms x y psi1 psi2 psi3 real
syms xd yd psi1d psi2d psi3d real

q = [x; y; psi1; psi2; psi3;];
qd = [xd; yd; psi1d; psi2d; psi3d;];

%%
M = [m1 + m2 + m3, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1), (l_j2g2*m2 + l_j2j3*m3)*sin(psi2), l_j3g3*m3*sin(psi3);
     0, m1 + m2 + m3, -(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1), -(l_j2g2*m2 + l_j2j3*m3)*cos(psi2), -l_j3g3*m3*cos(psi3);
     (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1), -(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1), iz1+m1*(l_j1g1^2)+(m2+m3)*(l_j1j2^2), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*cos(psi1-psi2), l_j3g3*l_j1j2*m3*cos(psi1-psi3);
     (l_j2g2*m2 + l_j2j3*m3)*sin(psi2), -(l_j2g2*m2 + l_j2j3*m3)*cos(psi2), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*cos(psi1-psi2), iz2+m2*(l_j2g2^2)+m3*(l_j2j3^2), l_j3g3*l_j2j3*m3*cos(psi2-psi3);
     l_j3g3*m3*sin(psi3), -l_j3g3*m3*cos(psi3), l_j3g3*l_j1j2*m3*cos(psi1-psi3), l_j3g3*l_j2j3*m3*cos(psi2-psi3), iz3+m3*(l_j3g3^2);]

C = [0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, l_j3g3*m3*cos(psi3)*psi3d;
     0, 0, (l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, (l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, l_j3g3*m3*sin(psi3)*psi3d;
     1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*cos(psi1)*psi1d, 1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*sin(psi1)*psi1d, -1/2*(l_j1g1*m1 + l_j1j2*(m2 + m3))*(xd*cos(psi1)+yd*sin(psi1))+1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi2d+1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi3d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-1/2*psi1d+psi2d), l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-1/2*psi1d+psi3d);
     1/2*(l_j2g2*m2 + l_j2j3*m3)*cos(psi2)*psi2d, 1/2*(l_j2g2*m2 + l_j2j3*m3)*sin(psi2)*psi2d, (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*(-psi1d+1/2*psi2d), -1/2*(l_j2g2*m2+l_j2j3*m3)*(xd*cos(psi2)+yd*sin(psi2))-1/2*(l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3)*sin(psi1-psi2)*psi1d+1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi3d, l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-1/2*psi2d+psi3d);
     1/2*l_j3g3*m3*cos(psi3)*psi3d, 1/2*l_j3g3*m3*sin(psi3)*psi3d, l_j3g3*l_j1j2*m3*sin(psi1-psi3)*(-psi1d+1/2*psi3d), l_j3g3*l_j2j3*m3*sin(psi2-psi3)*(-psi2d+1/2*psi3d), -1/2*l_j3g3*m3*(xd*cos(psi3)+yd*sin(psi3))-1/2*l_j3g3*l_j1j2*m3*sin(psi1-psi3)*psi1d-1/2*l_j3g3*l_j2j3*m3*sin(psi2-psi3)*psi2d;]

disp(simplify(C*qd));

D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, d1, -d1, 0;
     0, 0, -d1, d1+d2, -d2;
     0, 0, 0, -d2, d2;];

%%
L_long1 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           -h11, -h12, -h13, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, -h21, -h22, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, -h31, -h32, -h33;]

L_long2 = [cos(psi1), cos(psi1), cos(psi1), cos(psi2), cos(psi2), cos(psi3), cos(psi3), cos(psi3);
           sin(psi1), sin(psi1), sin(psi1), sin(psi2), sin(psi2), sin(psi3), sin(psi3), sin(psi3);
           h11, h12, h13, l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi2), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3), l_j1j2*sin(psi1-psi3);
           0, 0, 0, h21, h22, l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3), l_j2j3*sin(psi2-psi3);
           0, 0, 0, 0, 0, h31, h32, h33;]

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
Ru = diag(R);

Tdr = [Tdr1; Tdr2; Tdr3; Tdr4; Tdr5; Tdr6; Tdr7; Tdr8;];
Tdr_D = diag(Tdr);

ca = [ca1; ca2; ca3; ca4; ca5; ca6; ca7; ca8;];
Ca = diag(ca);

delta = [delta1; delta2; delta3; delta4; delta5; delta6; delta7; delta8;];
delta_s = [sin(delta1); sin(delta2); sin(delta3); sin(delta4); sin(delta5); sin(delta6); sin(delta7); sin(delta8);];
delta_c = [cos(delta1); cos(delta2); cos(delta3); cos(delta4); cos(delta5); cos(delta6); cos(delta7); cos(delta8);];
Delta_s = diag(delta_s);
Delta_c = diag(delta_c);

beta11 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a1*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h11*psi1d)); 
beta12 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a1*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h11*psi1d)); 
beta21 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a2*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h12*psi1d)); 
beta22 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a2*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h12*psi1d)); 
beta31 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a3*psi1d) / (xd*cos(psi1) + yd*sin(psi1) - h13*psi1d)); 
beta32 = atan((-xd*sin(psi1) + yd*cos(psi1) - l_j1a3*psi1d) / (xd*cos(psi1) + yd*sin(psi1) + h13*psi1d)); 

beta41 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a4*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d - h21*psi2d)); 
beta42 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a4*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d + h21*psi2d));
beta51 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d - h22*psi2d)); 
beta52 = atan((-xd*sin(psi2) + yd*cos(psi2) - l_j1j2*cos(psi1-psi2)*psi1d - l_j2a5*psi2d) / (xd*cos(psi2) + yd*sin(psi2) + l_j1j2*sin(psi1-psi2)*psi1d + h22*psi2d));

beta61 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a6*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h31*psi3d)); 
beta62 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a6*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h31*psi3d));
beta71 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a7*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h32*psi3d)); 
beta72 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a7*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h32*psi3d));
beta81 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a8*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d - h33*psi3d)); 
beta82 = atan((-xd*sin(psi3) + yd*cos(psi3) - l_j1j2*cos(psi1-psi3)*psi1d - l_j2j3*cos(psi2-psi3)*psi2d - l_j3a8*psi3d) / (xd*cos(psi3) + yd*sin(psi3) + l_j1j2*sin(psi1-psi3)*psi1d + l_j2j3*sin(psi2-psi3)*psi2d + h33*psi3d));

beta1 = [beta11; beta21; beta31; beta41; beta51; beta61; beta71; beta81]
beta2 = [beta12; beta22; beta32; beta42; beta52; beta62; beta72; beta82]

% M*qdd + C*qd + D*qd = L_long1*(Delta_c*inv(Ru)*Tdr - Delta_s*Ca*(delta-beta1)) + L_lat1*(Delta_c*inv(Ru)*Tdr_D*delta + Delta_s*Ca*(delta-beta1)) + L_long2*(Delta_c*inv(Ru)*Tdr - Delta_s*Ca*(delta-beta2)) + L_lat2*(Delta_c*inv(Ru)*Tdr_D*delta + Delta_s*Ca*(delta-beta2))

%%
Twv = [cos(psi1), -sin(psi1), 0, 0, 0;
     sin(psi1), cos(psi1), 0, 0, 0;
     0, 0, 1, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1;];
% need to check
Tdwv = diff(Twv) * psi1d;

Mv = Twv'*M*Twv;
Cv = Twv'*C*Twv + Twv'*M*Tdwv;
Dv = Twv'*D;
L_long1_v = Twv'*L_long1;
L_long2_v = Twv'*L_long2;
L_lat1_v = Twv'*L_lat1;
L_lat2_v = Twv'*L_lat2;
