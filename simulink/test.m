q = [1;2;3;4;5];
qd = [6;7;8;9;10;];

x = q(1);
y = q(2);
psi1 = q(3);

v_Rot_w = [cos(psi1), sin(psi1);
           -sin(psi1), cos(psi1)];

w_T_v = [cos(psi1), -sin(psi1), 0, 0, 0;
       sin(psi1), cos(psi1), 0, 0, 0;
       0, 0, 1, 0, 0;
       0, 0, 0, 1, 0;
       0, 0, 0, 0, 1;];

q_v = zeros(5, 1);
qd_v = zeros(5, 1);

q_v(1:2, :) = v_Rot_w * [x; y];
q_v(3:5, :) = q(3:5, :);

qd_v = inv(w_T_v) * qd;

x_v = q_v(1);
y_v = q_v(2);
psi1_v = q_v(3);

w_Rot_v = [cos(psi1_v), -sin(psi1_v);
           sin(psi1_v), cos(psi1_v);];

w_T_v = [cos(psi1_v), -sin(psi1_v), 0, 0, 0;
         sin(psi1_v), cos(psi1_v), 0, 0, 0;
         0, 0, 1, 0, 0;
         0, 0, 0, 1, 0;
         0, 0, 0, 0, 1;];

q_n = zeros(5, 1);
qd_n = zeros(5, 1);

q_n(1:2, :) = w_Rot_v * [x_v; y_v];
q_n(3:5, :) = q_v(3:5, :);

disp(q_n)
qd_n = w_T_v * qd_v

