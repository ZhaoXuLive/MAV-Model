function [psi1_des, psi2_des, psi3_des] = psi2ides(Scur_s1_des, Scur_s2_des, Scur_s3_des, Scur_s4_des, ts, TRACK)

ps1x = TRACK(2, Scur_s1_des(1)) * ts(1) + TRACK(2, Scur_s1_des(2)) * (1 - ts(1));
ps1y = TRACK(3, Scur_s1_des(1)) * ts(1) + TRACK(3, Scur_s1_des(2)) * (1 - ts(1));
ps2x = TRACK(2, Scur_s2_des(1)) * ts(2) + TRACK(2, Scur_s2_des(2)) * (1 - ts(2));
ps2y = TRACK(3, Scur_s2_des(1)) * ts(2) + TRACK(3, Scur_s2_des(2)) * (1 - ts(2));
ps3x = TRACK(2, Scur_s3_des(1)) * ts(3) + TRACK(2, Scur_s3_des(2)) * (1 - ts(3));
ps3y = TRACK(3, Scur_s3_des(1)) * ts(3) + TRACK(3, Scur_s3_des(2)) * (1 - ts(3));
ps4x = TRACK(2, Scur_s4_des(1)) * ts(4) + TRACK(2, Scur_s4_des(2)) * (1 - ts(4));
ps4y = TRACK(3, Scur_s4_des(1)) * ts(4) + TRACK(3, Scur_s4_des(2)) * (1 - ts(4));

psi1_des = atan2(ps2y - ps1y, ps2x - ps1x);
psi2_des = atan2(ps3y - ps2y, ps3x - ps2x);
psi3_des = atan2(ps4y - ps3y, ps4x - ps3x);

end