% the length of every vehicle unit
lu = [11.4, 9.4, 11.4];
% the length of every vehicle unit from the front joint to the mass center
lg = [5.862, 4.643, 5.538];

% 假设所有车体都有侧偏，即速度方向与运动方向不同
% 输入为vx1, vy1, x, y, psi1, psi2, psi3, delta
% 目标为xd, yd, psi1d, psi2d, psi3d
xd = vx1 * cos(psi1) - vy1 * sin(psi1);
yd = vx1 * sin(psi1) + vy1 * cos(psi1);
% 由于我们知道六个轴的输入转向角，根据公式可得到各车辆单元质心速度方向角
ps1d = vx1 * (tan(delta1) + tan(delta2)) / l_a1a2;
vx2 = vx1 * cos(psi1 - psi2) + (vy1 + l_g1j2 * psi1d) * sin(psi1 - psi2);
ps2d = vx2 * (tan(delta3) + tan(delta4)) / l_a3a4;
vy2 = vx1 * sin(psi1 - psi2) - (vy1 + l_g1j2 * psi1d) * cos(psi1 - psi2) - l_j2g2 * psi2d;
vx3 = vx2 * cos(psi2 - psi3) + (vy2 + l_g2j3 * psi2d) * sin(psi2 - psi3);
ps3d = vx3 * (tan(delta5) + tan(delta6)) / l_a5a6;


