function [w,L]=w8s(x,track,U,s_in,a,U0)

y0=0.1;
Lmax=20;
U=max(U,U0);
b=0;

sxy=x2ss(x,track,s_in);
trk=trkint(sxy(1),track); 
t1=trk(4:5);n1=trk(6:7);
gamma=1-sxy(2)*trk(8);

y1=abs(sxy(2));
if y1<y0
    flag=1;y=y0;
else
    flag=0;y=y1;
end
den=sqrt(2*a*y+b);
L=U*y/den;
L=min(L,Lmax);

q=atan(y/L);
if flag==1
    q=q*(y1/y0);
end

ws=[0 0];
ws(1)=U*cos(q);
ws(2)=-sign(sxy(2))*U*sin(q);

w=gamma*ws(1)*t1+ws(2)*n1;

return