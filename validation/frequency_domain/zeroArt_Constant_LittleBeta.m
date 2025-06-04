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

vx = 20 / 3.6;

%% linear matrix

M = [m1 + m2 + m3, -(l_j1g1*m1 + l_j1j2*(m2 + m3)), -(l_j2g2*m2 + l_j2j3*m3), -l_j3g3*m3;
     -(l_j1g1*m1 + l_j1j2*(m2 + m3)), Iz1+m1*(l_j1g1^2)+(m2+m3)*(l_j1j2^2), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3), l_j3g3*l_j1j2*m3;
     -(l_j2g2*m2 + l_j2j3*m3), (l_j2g2*l_j1j2*m2+l_j1j2*l_j2j3*m3), Iz2+m2*(l_j2g2^2)+m3*(l_j2j3^2), l_j3g3*l_j2j3*m3;
     -l_j3g3*m3, l_j3g3*l_j1j2*m3, l_j3g3*l_j2j3*m3, Iz3+m3*(l_j3g3^2);];

C = vx * [0, m1 + m2 + m3, 0, 0;
          0, -(l_j1g1*m1 + l_j1j2*(m2 + m3)), 0, 0;
          0, -(l_j2g2*m2 + l_j2j3*m3), 0, 0;
          0, -l_j3g3*m3, 0, 0;];

D = [0, 0, 0, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, d1, -d1, 0;
     0, 0, -d1, d1+d2, -d2;
     0, 0, 0, -d2, d2;];

L_lat = [ 1, 1, 1, 1, 1, 1, 1, 1;
          -l_j1a1, -l_j1a2, -l_j1a3, -l_j1j2, -l_j1j2, -l_j1j2, -l_j1j2, -l_j1j2;
          0, 0, 0, -l_j2a4, -l_j2a5, -l_j2j3, -l_j2j3, -l_j2j3;
          0, 0, 0, 0, 0, -l_j3a6, -l_j3a7, -l_j3a8;];

beta = 1 / vx * [1, -l_j1a1, 0, 0;
                 1, -l_j1a2, 0, 0;
                 1, -l_j1a3, 0, 0;
                 1, -l_j1j2, -l_j2a4, 0;
                 1, -l_j1j2, -l_j2a5, 0;
                 1, -l_j1j2, -l_j2j3, -l_j3a6;
                 1, -l_j1j2, -l_j2j3, -l_j3a7;
                 1, -l_j1j2, -l_j2j3, -l_j3a8;];

ca = [ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8]';
Ca = diag(ca);

E = 2 * L_lat * Ca * beta;
F = 2 * L_lat * Ca;

%% transfer function

A1 = zeros(8, 8);
B1 = zeros(8, 8);
C1 = zeros(4, 8);

A1(1:4, 5:8) = eye(4, 4);
A1(5:8, 5:8) = -inv(M) * (C + E);
B1(5:8, :) = inv(M) * F;
C1(:, 5:8) = eye(4, 4);
D1 = zeros(4, 8);

% C1 for first unit latAcc and three yawrate
% C2 for second unit latAcc
% C3 for third unit latAcc

C2 = zeros(4, 8);

sys = ss(A1,B1,C1,D1);
tfun = tf(sys)
w = logspace(-2, 4);

%% plot

% bode(tfun);
% freqs(tfun);

% h = freqresp(sys, w);
% mag = abs(h);
% phase = angle(h);
% phasedeg = phase*180/pi;

% mag11 = mag(1, 1, :);
% mag11 = mag11(:);
% subplot(2,1,1)
% semilogx(w, mag11)
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('Magnitude')
% 
% phasedeg11 = phasedeg(1, 1, :);
% phasedeg11 = phasedeg11(:);
% subplot(2,1,2)
% semilogx(w, phasedeg11)
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('Phase (degrees)')

% mag11 = mag(1, 1, :);
% mag11 = mag11(:);
% mag41 = mag(4, 1, :);
% mag41 = mag41(:);
% RA = zeros(50, 1);
% for i = 1:50
%     RA(i) = mag41(i) / mag11(i);
% end
% semilogx(w, RA);
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('RA')