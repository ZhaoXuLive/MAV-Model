time = out.tout;

% trajectory

figure(1)
plot(xg1, yg1)
hold on 
plot(xg2, yg2)
hold on 
plot(xg3, yg3)
title('unit cg trajectory')
xlabel('X(m)')
ylabel('Y(m)')
legend('unit1', 'unit2', 'unit3')
axis equal
% axis([0 350 0 350]);

figure(2)
plot(jointx1, jointy1)
hold on 
plot(jointx2, jointy2)
hold on 
plot(jointx3, jointy3)
hold on 
plot(jointx4, jointy4)
title('joint trajectory')
xlabel('X(m)')
ylabel('Y(m)')
legend('joint1', 'joint2', 'joint3', 'joint4')
axis equal
% axis([0 350 0 350]);

figure(3)
plot(axlex1, axley1)
hold on
plot(axlex2, axley2)
hold on
plot(axlex3, axley3)
hold on
plot(axlex4, axley4)
hold on
plot(axlex5, axley5)
hold on
plot(axlex6, axley6)
hold on
plot(axlex7, axley7)
hold on
plot(axlex8, axley8)
title('axle trajectory')
xlabel('X(m)')
ylabel('Y(m)')
legend('axle1', 'axle2', 'axle3', 'axle4', 'axle5', 'axle6', 'axle7', 'axle8')
axis equal
% axis([0 350 0 350]);

% angle

figure(4)
plot(time, psi1)
hold on
plot(time, psi2)
hold on
plot(time, psi3)
title('yaw angle')
xlabel('t(s)')
ylabel('deg(°)')
legend('unit1', 'unit2', 'unit3')

figure(5)
plot(time, refang1)
hold on
plot(time, refang2)
hold on
plot(time, refang3)
hold on
plot(time, refang4)
hold on
plot(time, refang5)
hold on
plot(time, refang6)
title('AFG reference angle')
xlabel('t(s)')
ylabel('deg(°)')
legend('angref1', 'angref2', 'angref3', 'angref4', 'angref5', 'angref6')

figure(6)
plot(time, steer1)
hold on
plot(time, steer2)
hold on
plot(time, steer3)
hold on
plot(time, steer4)
hold on
plot(time, steer5)
hold on
plot(time, steer6)
hold on
plot(time, steer7)
hold on
plot(time, steer8)
title('axle steer angle')
xlabel('t(s)')
ylabel('deg(°)')
legend('axle1', 'axle2', 'axle3', 'axle4', 'axle5', 'axle6', 'axle7', 'axle8')

figure(7)
plot(time, arta1)
hold on
plot(time, arta2)
title('Articulation Angle')
xlabel('t(s)')
ylabel('deg(°)')
legend('first', 'second')

% velocity 

figure(8)
plot(time, dpsi1)
hold on
plot(time, dpsi2)
hold on
plot(time, dpsi3)
title('yaw rate')
xlabel('t(s)')
ylabel('w(°/s)')
legend('unit1', 'unit2', 'unit3')

figure(9)
plot(time, longitudinalu1)
hold on
plot(time, longitudinalu2)
hold on
plot(time, longitudinalu3)
title('unit longitudinal velocity')
xlabel('t(s)')
ylabel('v(m/s)')
legend('unit1', 'unit2', 'unit3')

figure(10)
plot(time, lateralu1)
hold on
plot(time, lateralu2)
hold on
plot(time, lateralu3)
title('unit lateral velocity')
xlabel('t(s)')
ylabel('v(m/s)')
legend('unit1', 'unit2', 'unit3')

figure(11)
plot(time, djointx1)
hold on
plot(time, djointx2)
hold on
plot(time, djointx3)
hold on
plot(time, djointx4)
title('joint global x velocity')
xlabel('t(s)')
ylabel('v(m/s)')
legend('joint1', 'joint2', 'joint3', 'joint4')

figure(12)
plot(time, djointy1)
hold on
plot(time, djointy2)
hold on
plot(time, djointy3)
hold on
plot(time, djointy4)
title('joint global y velocity')
xlabel('t(s)')
ylabel('v(m/s)')
legend('joint1', 'joint2', 'joint3', 'joint4')

% figure(13)
% plot(time, Fx1)
% hold on
% plot(time, Fx2)
% hold on
% plot(time, Fx3)
% hold on
% plot(time, Fx4)
% hold on
% plot(time, Fx5)
% hold on
% plot(time, Fx6)
% hold on
% plot(time, Fx7)
% hold on
% plot(time, Fx8)
% title('Fx')
% xlabel('t(s)')
% ylabel('F(N)')
% legend('Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Fx6', 'Fx7', 'Fx8')