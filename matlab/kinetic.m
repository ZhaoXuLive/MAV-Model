% the length of every vehicle unit
lu = [11.4, 9.4, 11.4];
% the length of every vehicle unit from the front joint to the mass center
lg = [5.862, 4.643, 5.538];

% % formula explain
% dx = v1 * cos(psi1);
% dy = v1 * sin(psi1);
% dpsi1 = (v1 * sin(r1) + lg(2) * dr1) / ((lu(1) - lg(1)) * cos(r1) + lg(2));
% dpsi2 = (v2 * sin(r2) + lg(3) * dr2) / ((lu(2) - lg(2)) * cos(r2) + lg(3));
% v2 = v1 * cos(r1) - (lu(1) - lg(1)) * dpsi1 * sin(r1);
% dr2 = f(r1, r2, delta);