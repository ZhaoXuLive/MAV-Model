% 魔术公式参数
% A = [1.65, -34, 1250, 3036, 12.8, 0.00501, -0.02103,  0.77394, ...
%      0.0022890, 0.013442, 0.003709, 19.1656, 1.21356, 6.26206]';
A = [1.65, -34, 1250, 4535, 12.8, 0.00501, -0.02103,  0.77394, ...
     0.0022890, 0.013442, 0.003709, 19.1656, 1.21356, 6.26206]';

% 轮胎垂向载荷 N
Fz_list = [10000, 15000, 20000, 25000, 30000] / 1000;

% 轮胎侧偏角 °
alpha_list = -8:0.25:8;

Fy = zeros(65, 5);
for i = 1:65
    alpha = alpha_list(i);
    for j = 1:5
        Fz = Fz_list(j);
        Fy(i, j) = latforce(Fz, 0, alpha, A);
    end
end

plot(alpha_list, Fy(:, 5));

% plot(alpha_list, Fy(:, 1));
% hold on
% plot(alpha_list, Fy(:, 2));
% hold on
% plot(alpha_list, Fy(:, 3));
% hold on
% plot(alpha_list, Fy(:, 4));
% hold on
% plot(alpha_list, Fy(:, 5));
% hold on
% title('magic formula')
% xlabel('alpha(°)')
% ylabel('latforce(KN)')
% legend('Fz=10KN', 'Fz=15KN', 'Fz=20KN', 'Fz=25KN', 'Fz=30KN')


%% lat force
function Fy = latforce(Fz, gamma, alpha, A)

% 水平方向漂移
Sh = A(10) * Fz + A(11) + A(9) * gamma;

% 垂直方向漂移
Sv = A(12) * Fz * gamma + A(13) * Fz + A(14);

% 形状因子
C = A(1);

% 峰值因子
D = A(2) * Fz ^ 2 + A(3) * Fz;

% 侧向力零点处的侧向刚度
BCD = A(4) * sin(2 * atan(Fz / A(5))) * (1 - A(6) * abs(gamma));

% 刚度因子
B = BCD / (C * D);

% 曲率因子
E = A(7) * Fz + A(8);

% 输入
x = alpha; % + Sh;

% 轮胎侧向力
Fy = D * sin(C * atan(B * x - E * (B * x - atan(B * x)))); % + Sv;

end