function qd_output = kinetic_model(q, vx, vy, delta)

psi1 = q(3);
psi2 = q(4);
psi3 = q(5);

% fixed parameter
l_g1j2 = 5.538;
l_j2g2 = 9.4 - 4.757;
l_g2j3 = 4.757;
l_a1a3 = 7.4;
l_a4a5 = 6.4;
l_a6a8 = 7.4;

% 假设所有车体都有侧偏，即速度方向与运动方向不同
% 输入为vx1, vy1, x, y, psi1, psi2, psi3, delta
vx1 = vx * cos(psi1) + vy * sin(psi1);
vy1 = -vx * sin(psi1) + vy * cos(psi1);

% 目标为xd, yd, psi1d, psi2d, psi3d
% 由于我们知道六个轴的输入转向角，根据公式可得到各车辆单元质心速度方向角
psi1d = vx1 * (tan(delta(1)) + tan(delta(3))) / l_a1a3;
vx2 = vx1 * cos(psi1 - psi2) + (-vy1 + l_g1j2 * psi1d) * sin(psi1 - psi2);
psi2d = vx2 * (tan(delta(4)) + tan(delta(5))) / l_a4a5;
vy2 = vx1 * sin(psi1 - psi2) - (-vy1 + l_g1j2 * psi1d) * cos(psi1 - psi2) - l_j2g2 * psi2d;
vx3 = vx2 * cos(psi2 - psi3) + (-vy2 + l_g2j3 * psi2d) * sin(psi2 - psi3);
psi3d = vx3 * (tan(delta(6)) + tan(delta(8))) / l_a6a8;

xd = vx - 5.862 * sin(psi1) * psi1d;
yd = vy + 5.862 * cos(psi1) * psi1d;
% xd = vx;
% yd = vy;

qd_output = [xd; yd; psi1d; psi2d; psi3d];