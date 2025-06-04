function [Scur, psi_p, t] = projectPsi(x, y, S0, Srange, TRACK)

minDistance = 500;
t = 1;
for i = -Srange:+Srange
    if hypot(TRACK(2, S0+i) - x, TRACK(3, S0+i) - y) < minDistance
        minDistance = hypot(TRACK(2, S0+i) - x, TRACK(3, S0+i) - y);
        Scur(1) = S0+i;
    end
end

if hypot(TRACK(2, Scur(1)-1) - x, TRACK(3, Scur(1)-1) - y) > ...
        hypot(TRACK(2, Scur(1)+1) - x, TRACK(3, Scur(1)+1) - y)
    Scur(2) = Scur(1) + 1;
elseif hypot(TRACK(2, Scur(1)-1) - x, TRACK(3, Scur(1)-1) - y) < ...
        hypot(TRACK(2, Scur(1)+1) - x, TRACK(3, Scur(1)+1) - y)
    Scur(2) = Scur(1) - 1;
else
    Scur(2) = -1;
end

psi_pj1_1 = atan2(TRACK(5, Scur(1)), TRACK(4, Scur(1)));
if Scur(2) == -1
    psi_p = psi_pj1_1;
else
    psi_pj1_2 = atan2(TRACK(5, Scur(2)), TRACK(4, Scur(2)));
    t = hypot(TRACK(2, Scur(1)) - x, TRACK(3, Scur(1)) - y) \ ...
        (hypot(TRACK(2, Scur(1)) - x, TRACK(3, Scur(1)) - y) + ...
        hypot(TRACK(2, Scur(2)) - x, TRACK(3, Scur(2)) - y));
    psi_p = psi_pj1_1 * t + psi_pj1_2 * (1 - t);
end

end