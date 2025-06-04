function [s,eflag]=x2s(x,track)
s=[0,0];
Ds=500;
small=1e-3;

test=norm(track(1,2:3)-track(end,2:3))+norm(track(1,4:5)-track(end,4:5));
if test<small 
    TRACK=track;
else    
    row0=track(1,:);x1=track(1,2:3);t1=track(1,4:5);
    row0(1)=row0(1)-Ds;row0(2:3)=x1-Ds*t1;row0(8)=0;
    rowz=track(end,:);xn=track(end,2:3);tn=track(end,4:5);
    rowz(1)=rowz(1)+Ds;rowz(2:3)=xn+Ds*tn;rowz(8)=0;
    track(end,8)=0;
    TRACK=[row0;track;rowz];
end

[s,eflag]=x2s0(x,TRACK);
end