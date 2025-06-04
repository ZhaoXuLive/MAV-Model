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
% basic equation
P = [1, 0, l_j1g1*sin(psi1), 0, 0;
     0, 1, -l_j1g1*cos(psi1), 0, 0;
     1, 0, l_j1j2*sin(psi1), l_j2g2*sin(psi2), 0;
     0, 1, -l_j1j2*cos(psi1), -l_j2g2*cos(psi2), 0;
     1, 0, l_j1j2*sin(psi1), l_j2j3*sin(psi2), l_j3g3*sin(psi3);
     0, 1, -l_j1j2*cos(psi1), -l_j2j3*cos(psi2), -l_j3g3*cos(psi3);];

Pd = [0, 0, l_j1g1*cos(psi1)*psi1d, 0, 0;
      0, 0, l_j1g1*sin(psi1)*psi1d, 0, 0;
      0, 0, l_j1j2*cos(psi1)*psi1d, l_j2g2*cos(psi2)*psi2d, 0;
      0, 0, l_j1j2*sin(psi1)*psi1d, l_j2g2*sin(psi2)*psi2d, 0;
      0, 0, l_j1j2*cos(psi1)*psi1d, l_j2j3*cos(psi2)*psi2d, l_j3g3*cos(psi3)*psi3d;
      0, 0, l_j1j2*sin(psi1)*psi1d, l_j2j3*sin(psi2)*psi2d, l_j3g3*sin(psi3)*psi3d;];

R = [0, 0, l_j1g1*cos(psi1), 0, 0;
     0, 0, l_j1g1*sin(psi1), 0, 0;
     0, 0, l_j1j2*cos(psi1), l_j2g2*cos(psi2), 0;
     0, 0, l_j1j2*sin(psi1), l_j2g2*sin(psi2), 0;
     0, 0, l_j1j2*cos(psi1), l_j2j3*cos(psi2), l_j3g3*cos(psi3);
     0, 0, l_j1j2*sin(psi1), l_j2j3*sin(psi2), l_j3g3*sin(psi3);];

% Rd = (R'.*qd)'

m = [m1; m1; m2; m2; m3; m3;];
D1 = diag(m);

i = [0; 0; iz1; iz2; iz3];
D2 = diag(i);

M = simplify(P'*D1*P + D2)
% problem 
Md = simplify(Pd'*D1*P + P'*D1*Pd);

Cq1 = simplify(Md*qd)
Cq2 = simplify((qd'*P'*D1*R)'.*qd)

Cq3 = simplify(qd'*P'*D1*Pd*qd)

C = simplify(Cq1 - Cq2)

% %%
% J1 = [1, 0, 1, 0, 1, 0;
%       0, 1, 0, 1, 0, 1;
%       l_j1a1*sin(psi1), -l_j1a1*cos(psi1), l_j1a2*sin(psi1), -l_j1a2*cos(psi1), l_j1a3*sin(psi1), -l_j1a3*cos(psi1);
%       0, 0, 0, 0, 0, 0;
%       0, 0, 0, 0, 0, 0;];
% J2 = [1, 0, 1, 0;
%       0, 1, 0, 1;
%       l_j1j2*sin(psi1), -l_j1j2*cos(psi1), l_j1j2*sin(psi1), -l_j1j2*cos(psi1);
%       l_j2a4*sin(psi2), -l_j2a4*cos(psi2), l_j2a5*sin(psi2), -l_j2a5*cos(psi2);
%       0, 0, 0, 0;];
% J3 = [1, 0, 1, 0, 1, 0;
%       0, 1, 0, 1, 0, 1;
%       l_j1j2*sin(psi1), -l_j1j2*cos(psi1), l_j1j2*sin(psi1), -l_j1j2*cos(psi1), l_j1j2*sin(psi1), -l_j1j2*cos(psi1);
%       l_j2j3*sin(psi2), -l_j2j3*cos(psi2), l_j2j3*sin(psi2), -l_j2j3*cos(psi2), l_j2j3*sin(psi2), -l_j2j3*cos(psi2);
%       l_j3a6*sin(psi3), -l_j3a6*cos(psi3), l_j3a7*sin(psi3), -l_j3a7*cos(psi3), l_j3a8*sin(psi3), -l_j3a8*cos(psi3);];
% J = [J1, J2, J3];
% 
% vx1 = xd + sin(psi1)*psi1d*l_j1a1;
% vy1 = yd - cos(psi1)*psi1d*l_j1a1;
% vxt1 = cos(psi1 + delta1) * vx1 + sin(psi1 + delta1) * vy1;
% vyt1 = -sin(psi1 + delta1) * vx1 + cos(psi1 + delta1) * vy1;
% alpha1 = atan(vyt1, vxt1);
% Fyt1 = ca1*alpha1; 
% Fxt1 = Tdr1 / r1;
% Fx1 = cos(psi1 + delta1) * Fxt1 - sin(psi1 + delta1) * Fyt1;
% Fy1 = sin(psi1 + delta1) * Fxt1 + cos(psi1 + delta1) * Fyt1;
% 
% vx2 = xd + sin(psi1)*psi1d*l_j1a2;
% vy2 = yd - cos(psi1)*psi1d*l_j1a2;
% vxt2 = cos(psi1 + delta2) * vx2 + sin(psi1 + delta2) * vy2;
% vyt2 = -sin(psi1 + delta2) * vx2 + cos(psi1 + delta2) * vy2;
% alpha2 = atan(vyt2, vxt2);
% Fyt2 = ca2*alpha2; 
% Fxt2 = Tdr2 / r2;
% Fx2 = cos(psi1 + delta2) * Fxt2 - sin(psi1 + delta2) * Fyt2;
% Fy2 = sin(psi1 + delta2) * Fxt2 + cos(psi1 + delta2) * Fyt2;
% 
% vx3 = xd + sin(psi1)*psi1d*l_j1a3;
% vy3 = yd - cos(psi1)*psi1d*l_j1a3;
% vxt3 = cos(psi1 + delta3) * vx3 + sin(psi1 + delta3) * vy3;
% vyt3 = -sin(psi1 + delta3) * vx3 + cos(psi1 + delta3) * vy3;
% alpha3 = atan(vyt3, vxt3);
% Fyt3 = ca3*alpha3; 
% Fxt3 = Tdr3 / r3;
% Fx3 = cos(psi1 + delta3) * Fxt3 - sin(psi1 + delta3) * Fyt3;
% Fy3 = sin(psi1 + delta3) * Fxt3 + cos(psi1 + delta3) * Fyt3;
% 
% vx4 = xd + sin(psi1)*psi1d*l_j1j2 + sin(psi2)*psi2d*l_j2a4;
% vy4 = yd - cos(psi1)*psi1d*l_j1j2 - cos(psi2)*psi2d*l_j2a4;
% vxt4 = cos(psi2 + delta4) * vx4 + sin(psi2 + delta4) * vy4;
% vyt4 = -sin(psi2 + delta4) * vx4 + cos(psi2 + delta4) * vy4;
% alpha4 = atan(vyt4, vxt4);
% Fyt4 = ca4*alpha4; 
% Fxt4 = Tdr4 / r4;
% Fx4 = cos(psi2 + delta4) * Fxt4 - sin(psi2 + delta4) * Fyt4;
% Fy4 = sin(psi2 + delta4) * Fxt4 + cos(psi2 + delta4) * Fyt4;
% 
% vx5 = xd + sin(psi1)*psi1d*l_j1j2 + sin(psi2)*psi2d*l_j2a5;
% vy5 = yd - cos(psi1)*psi1d*l_j1j2 - cos(psi2)*psi2d*l_j2a5;
% vxt5 = cos(psi2 + delta5) * vx5 + sin(psi2 + delta5) * vy5;
% vyt5 = -sin(psi2 + delta5) * vx5 + cos(psi2 + delta5) * vy5;
% alpha5 = atan(vyt5, vxt5);
% Fyt5 = ca5*alpha5; 
% Fxt5 = Tdr5 / r5;
% Fx5 = cos(psi2 + delta5) * Fxt5 - sin(psi2 + delta5) * Fyt5;
% Fy5 = sin(psi2 + delta5) * Fxt5 + cos(psi2 + delta5) * Fyt5;
% 
% vx6 = xd + sin(psi1)*psi1d*l_j1j2 + sin(psi2)*psi2d*l_j2j3 + sin(psi3)*psi3d*l_j3a6;
% vy6 = yd - cos(psi1)*psi1d*l_j1j2 - cos(psi2)*psi2d*l_j2j3 - cos(psi3)*psi3d*l_j3a6;
% vxt6 = cos(psi3 + delta6) * vx6 + sin(psi3 + delta6) * vy6;
% vyt6 = -sin(psi3 + delta6) * vx6 + cos(psi3 + delta6) * vy6;
% alpha6 = atan(vyt6, vxt6);
% Fyt6 = ca6*alpha6; 
% Fxt6 = Tdr6 / r6;
% Fx6 = cos(psi3 + delta6) * Fxt6 - sin(psi3 + delta6) * Fyt6;
% Fy6 = sin(psi3 + delta6) * Fxt6 + cos(psi3 + delta6) * Fyt6;
% 
% vx7 = xd + sin(psi1)*psi1d*l_j1j2 + sin(psi2)*psi2d*l_j2j3 + sin(psi3)*psi3d*l_j3a7;
% vy7 = yd - cos(psi1)*psi1d*l_j1j2 - cos(psi2)*psi2d*l_j2j3 - cos(psi3)*psi3d*l_j3a7;
% vxt7 = cos(psi3 + delta7) * vx7 + sin(psi3 + delta7) * vy7;
% vyt7 = -sin(psi3 + delta7) * vx7 + cos(psi3 + delta7) * vy7;
% alpha7 = atan(vyt7, vxt7);
% Fyt7 = ca7*alpha7; 
% Fxt7 = Tdr7 / r7;
% Fx7 = cos(psi3 + delta7) * Fxt7 - sin(psi3 + delta7) * Fyt7;
% Fy7 = sin(psi3 + delta7) * Fxt7 + cos(psi3 + delta7) * Fyt7;
% 
% vx8 = xd + sin(psi1)*psi1d*l_j1j2 + sin(psi2)*psi2d*l_j2j3 + sin(psi3)*psi3d*l_j3a8;
% vy8 = yd - cos(psi1)*psi1d*l_j1j2 - cos(psi2)*psi2d*l_j2j3 - cos(psi3)*psi3d*l_j3a8;
% vxt8 = cos(psi3 + delta8) * vx8 + sin(psi3 + delta8) * vy8;
% vyt8 = -sin(psi3 + delta8) * vx8 + cos(psi3 + delta8) * vy8;
% alpha8 = atan(vyt8, vxt8);
% Fyt8 = ca8*alpha8; 
% Fxt8 = Tdr8 / r8;
% Fx8 = cos(psi3 + delta8) * Fxt8 - sin(psi3 + delta8) * Fyt8;
% Fy8 = sin(psi3 + delta8) * Fxt8 + cos(psi3 + delta8) * Fyt8;
% 
% F = [Fx1; Fy1; Fx2; Fy2; Fx3; Fy3; Fx4; Fy4; Fx5; Fy5; Fx6; Fy6; Fx7; Fy7; Fx8; Fy8;];
% 
% Q = J*F

% M*qdd + C = Q