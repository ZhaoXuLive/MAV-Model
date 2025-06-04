sign = 1;

% 车辆质心与铰接点之间纵向距离
a1 = 11.4 - 5.538;
b1 = 5.538;
a2 = 4.643;
b2 = 4.757;
a3 = 5.538;
b3 = 11.4 - 5.538;

% 各轴相对于前虚拟铰接点的距离
l1 = 2.5;
l2 = 4;
l3 = 9.9;
l4 = 1.5;
l5 = 7.9;
l6 = 1.5;
l7 = 7.4;
l8 = 8.9;

% load first virtual joint location
if sign == 1
    jx1 = out.q_out(:,1);
    jy1 = out.q_out(:,2);
else
    jx1 = out.q_out(1,:,:);
    jy1 = out.q_out(2,:,:);
end
jx1 = jx1(:);
jy1 = jy1(:);

% load first virtual joint velocity
if sign == 1
    djx1 = out.qd_out(:,1);
    djy1 = out.qd_out(:,2);
else
    djx1 = out.qd_out(1,:,:);
    djy1 = out.qd_out(2,:,:);
end
djx1 = djx1(:);
djy1 = djy1(:);

% load yaw angle 
if sign == 1
    psi1 = out.q_out(:,3);
    psi2 = out.q_out(:,4);
    psi3 = out.q_out(:,5);
else
    psi1 = out.q_out(3,:,:);
    psi2 = out.q_out(4,:,:);
    psi3 = out.q_out(5,:,:);
end
psi1 = psi1(:);
psi2 = psi2(:);
psi3 = psi3(:);

% load yaw rate 
if sign == 1
    dpsi1 = out.qd_out(:,3);
    dpsi2 = out.qd_out(:,4);
    dpsi3 = out.qd_out(:,5);
else
    dpsi1 = out.qd_out(3,:,:);
    dpsi2 = out.qd_out(4,:,:);
    dpsi3 = out.qd_out(5,:,:);
end
dpsi1 = dpsi1(:);
dpsi2 = dpsi2(:);
dpsi3 = dpsi3(:);

% load side slip angle
alpha1 = out.alpha_out(:, 1)*180/pi;
alpha2 = out.alpha_out(:, 2)*180/pi;
alpha3 = out.alpha_out(:, 3)*180/pi;
alpha4 = out.alpha_out(:, 4)*180/pi;
alpha5 = out.alpha_out(:, 5)*180/pi;
alpha6 = out.alpha_out(:, 6)*180/pi;
alpha7 = out.alpha_out(:, 7)*180/pi;
alpha8 = out.alpha_out(:, 8)*180/pi;

% load steer angle
steer1 = out.steer(1, :, :);
steer2 = out.steer(2, :, :);
steer3 = out.steer(3, :, :);
steer4 = out.steer(4, :, :);
steer5 = out.steer(5, :, :);
steer6 = out.steer(6, :, :);
steer7 = out.steer(7, :, :);
steer8 = out.steer(8, :, :);
steer1 = steer1(:)*180/pi;
steer2 = steer2(:)*180/pi;
steer3 = steer3(:)*180/pi;
steer4 = steer4(:)*180/pi;
steer5 = steer5(:)*180/pi;
steer6 = steer6(:)*180/pi;
steer7 = steer7(:)*180/pi;
steer8 = steer8(:)*180/pi;

% load AFG reference angle
refang1 = out.ANGref(:, 1)*180/pi;
refang2 = out.ANGref(:, 2)*180/pi;
refang3 = out.ANGref(:, 3)*180/pi;
refang4 = out.ANGref(:, 4)*180/pi;
refang5 = out.ANGref(:, 5)*180/pi;
refang6 = out.ANGref(:, 6)*180/pi;

% calculate three unit cg location
xg1 = jx1 - a1 * cos(psi1);
yg1 = jy1 - a1 * sin(psi1);
xg2 = jx1 - (a1 + b1) * cos(psi1) - a2 * cos(psi2);
yg2 = jy1 - (a1 + b1) * sin(psi1) - a2 * sin(psi2);
xg3 = jx1 - (a1 + b1) * cos(psi1) - (a2 + b2) * cos(psi2) - a3 * cos(psi3);
yg3 = jy1 - (a1 + b1) * sin(psi1) - (a2 + b2) * sin(psi2) - a3 * sin(psi3);

% calculate three unit cg velocity
dxg1 = djx1 + a1 * sin(psi1) .* dpsi1;
dyg1 = djy1 - a1 * cos(psi1) .* dpsi1;
dxg2 = djx1 + (a1 + b1) * sin(psi1) .* dpsi1 + a2 * sin(psi2) .* dpsi2;
dyg2 = djy1 - (a1 + b1) * cos(psi1) .* dpsi1 - a2 * cos(psi2) .* dpsi2;
dxg3 = djx1 + (a1 + b1) * sin(psi1) .* dpsi1 + (a2 + b2) * sin(psi2) .* ...
    dpsi2 + a3 * sin(psi3) .* dpsi3;
dyg3 = djy1 - (a1 + b1) * cos(psi1) .* dpsi1 - (a2 + b2) * cos(psi2) .* ...
    dpsi2 - a3 * cos(psi3) .* dpsi3;

% calculate four joint location
jointx1 = jx1;
jointy1 = jy1;
jointx2 = jointx1 - (a1 + b1) * cos(psi1);
jointy2 = jointy1 - (a1 + b1) * sin(psi1);
jointx3 = jointx1 - (a1 + b1) * cos(psi1) - (a2 + b2) * cos(psi2);
jointy3 = jointy1 - (a1 + b1) * sin(psi1) - (a2 + b2) * sin(psi2);
jointx4 = jointx1 - (a1 + b1) * cos(psi1) - (a2 + b2) * cos(psi2) - (a3 + b3) * cos(psi3);
jointy4 = jointy1 - (a1 + b1) * sin(psi1) - (a2 + b2) * sin(psi2) - (a3 + b3) * sin(psi3);

% calculate four joint velocity
djointx1 = djx1;
djointy1 = djy1;
djointx2 = djointx1 + (a1 + b1) * sin(psi1) .* dpsi1;
djointy2 = djointy1 - (a1 + b1) * cos(psi1) .* dpsi1;
djointx3 = djointx1 + (a1 + b1) * sin(psi1) .* dpsi1 + (a2 + b2) * ...
    sin(psi2) .* dpsi2;
djointy3 = djointy1 - (a1 + b1) * cos(psi1) .* dpsi1 - (a2 + b2) * ...
    cos(psi2) .* dpsi2;
djointx4 = djointx1 + (a1 + b1) * sin(psi1) .* dpsi1 + (a2 + b2) * ...
    sin(psi2) .* dpsi2 + (a3 + b3) * sin(psi3) .* dpsi3;
djointy4 = djointy1 - (a1 + b1) * cos(psi1) .* dpsi1 - (a2 + b2) * ...
    cos(psi2) .* dpsi2 - (a3 + b3) * cos(psi3) .* dpsi3;

% calculate eight axle location
axlex1 = jointx1 - l1 * cos(psi1);
axley1 = jointy1 - l1 * sin(psi1);
axlex2 = jointx1 - l2 * cos(psi1);
axley2 = jointy1 - l2 * sin(psi1);
axlex3 = jointx1 - l3 * cos(psi1);
axley3 = jointy1 - l3 * sin(psi1);
axlex4 = jointx2 - l4 * cos(psi2);
axley4 = jointy2 - l4 * sin(psi2);
axlex5 = jointx2 - l5 * cos(psi2);
axley5 = jointy2 - l5 * sin(psi2);
axlex6 = jointx3 - l6 * cos(psi3);
axley6 = jointy3 - l6 * sin(psi3);
axlex7 = jointx3 - l7 * cos(psi3);
axley7 = jointy3 - l7 * sin(psi3);
axlex8 = jointx3 - l8 * cos(psi3);
axley8 = jointy3 - l8 * sin(psi3);

% calculate Articulation Angle
arta1 = (pi - (psi1 - psi2))*180/pi;
arta2 = (pi - (psi2 - psi3))*180/pi;

% calculate unit lateral/longitudinal velocity
longitudinalu1 = dxg1 .* cos(psi1) + dyg1 .* sin(psi1);
lateralu1 = -dxg1 .* sin(psi1) + dyg1 .* cos(psi1);
longitudinalu2 = dxg2 .* cos(psi2) + dyg2 .* sin(psi2);
lateralu2 = -dxg2 .* sin(psi2) + dyg2 .* cos(psi2);
longitudinalu3 = dxg3 .* cos(psi3) + dyg3 .* sin(psi3);
lateralu3 = -dxg3 .* sin(psi3) + dyg3 .* cos(psi3);

psi1 = psi1(:)*180/pi;
psi2 = psi2(:)*180/pi;
psi3 = psi3(:)*180/pi;

% % calculate Fx
% Fx1 = out.Fx_out(:, 1);
% Fx2 = out.Fx_out(:, 2);
% Fx3 = out.Fx_out(:, 3);
% Fx4 = out.Fx_out(:, 4);
% Fx5 = out.Fx_out(:, 5);
% Fx6 = out.Fx_out(:, 6);
% Fx7 = out.Fx_out(:, 7);
% Fx8 = out.Fx_out(:, 8);