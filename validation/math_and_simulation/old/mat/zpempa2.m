function [y,dy]=zpempa2(x,n)


[n1,n2]=size(x);
flag=0;
if n1<n2, x=x';flag=1; end
[n1,n2]=size(x);

n3=n+3; 
one=ones(n3,1);

B=x(1:n3,:);
ts=(1:n3)'; 
A=[one,ts];
P=A\B; 
AA=[one, ts-n3];
sx=AA*P;

B=x(n1-n3+1:n1,:);
tf=(n1-n3+1:n1)';
A=[one,tf];
P=A\B;
AA=[one, tf+n3];
fx=AA*P;


X=[sx;x;fx];
b=ones(1,n)/n;
Xr=flipud(X);
y1= filter(b,1,X);
y2= filter(b,1,Xr);
y2=flipud(y2);
y=(y1+y2)/2;

dy=(y(3:end,:)-y(1:end-2,:))/2;
y=y(n3+1:n3+n1,:);
dy=dy(n3:n3+n1-1,:);

if flag==1
    y=y';
    dy=dy';
end

return