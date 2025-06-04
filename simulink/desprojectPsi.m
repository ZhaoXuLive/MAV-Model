function [Scur_des, psi_des, t] = desprojectPsi(Scur0_des, t, L, Srange, TRACK)

minDistance = 500;
if Scur0_des(2) == -1
    Scur0x = TRACK(2, Scur0_des(1));
    Scur0y = TRACK(3, Scur0_des(1));
else
    Scur0x = t * TRACK(2, Scur0_des(1)) + (1 - t) * TRACK(2, Scur0_des(2));
    Scur0y = t * TRACK(3, Scur0_des(1)) + (1 - t) * TRACK(3, Scur0_des(2));
end

for i = Scur0_des(1) - 2 * Srange : Scur0_des(1)
    if abs(hypot(TRACK(2, i) - Scur0x, TRACK(3, i) - Scur0y) - L) < minDistance 
        minDistance = abs(hypot(TRACK(2, i) - Scur0x, TRACK(3, i) - Scur0y) - L);
        Scur(1) = i;
    end
end

if hypot(TRACK(2, Scur(1)-1) - Scur0x, TRACK(3, Scur(1)-1) - Scur0y) - L > ...
        hypot(TRACK(2, Scur(1)+1) - Scur0x, TRACK(3, Scur(1)+1) - Scur0y) - L
    Scur(2) = Scur(1) + 1;
elseif hypot(TRACK(2, Scur(1)-1) - Scur0x, TRACK(3, Scur(1)-1) - Scur0y) - L < ...
        hypot(TRACK(2, Scur(1)+1) - Scur0x, TRACK(3, Scur(1)+1) - Scur0y) - L
    Scur(2) = Scur(1) - 1;
else
    Scur(2) = -1;
end

psi_des_1 = atan2(TRACK(5, Scur(1)), TRACK(4, Scur(1)));
if Scur(2) == -1
    psi_des = psi_des_1;
else
    psi_des_2 = atan2(TRACK(5, Scur(2)), TRACK(4, Scur(2)));
    t = hypot(TRACK(2, Scur(1)) - x, TRACK(3, Scur(1)) - y) - L \ ...
        (hypot(TRACK(2, Scur(1)) - x, TRACK(3, Scur(1)) - y) - L + ...
        hypot(TRACK(2, Scur(2)) - x, TRACK(3, Scur(2)) - y)) - L;
    psi_des = psi_des_1 * t + psi_des_2 * (1 - t);
end

end