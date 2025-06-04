function [row,k]=trkint(s,track) 


small=1e-3;
n=size(track,1);
s=s(1);

test=norm(track(1,2:3)-track(n,2:3))+norm(track(1,4:5)-track(n,4:5));
if test<small
    s=wrap_cts(s,track(n,1));
end
     
S=track(:,1);I=(1:n)';
j=interp1(S,I,s,'linear','extrap');

if j<1
    row=track(1,:);k=0;
    row(1)=s;ds=s-track(1,1);
    row(2:3)=row(2:3)+ds*row(4:5);
    row(8)=0;
    return
elseif j>n-small
    row=track(n,:);k=n;
    row(1)=s;ds=s-track(n,1);
    row(2:3)=row(2:3)+ds*row(4:5);
    row(8)=0;
    return
end
k=floor(j);
dj=j-k;
row0=track(k,:);
row1=track(k+1,:);
ds=dj*(row1(1)-row0(1));
curv=row0(8);

phi=ds*curv/2;
if abs(phi)>small
    fac=sin(phi)/phi;
else
    fac=1-phi^2/6;
end

d=fac*ds;

r0=row0(2:3);t0=row0(4:5);n0=row0(6:7);
r=r0+d*(cos(phi)*t0+sin(phi)*n0);

t=t0*cos(2*phi)+n0*sin(2*phi);
n=[-t(2),t(1)];
row=[s,r,t,n,curv];
end



