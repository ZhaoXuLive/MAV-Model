function track=build4(p0,t0,s,c,d0)

if ~exist('d0','var'), d0=[];end
if isempty(d0), d0=0; end

tiny=1e-6;
np=length(s);
ds=diff(s);
c=c(1:np-1); 
ha=cumsum(c.*ds); 
ha1=atan2(t0(2),t0(1));
ha=[ha1;ha+ha1];
ha=wrap_angle(ha);

t=[cos(ha) sin(ha)]; 
n=[-sin(ha) cos(ha)]; 

dq=diff(ha);
dq=wrap_angle(dq);

p=zeros(np,2);
p(1,:)=p0;
for i=1:np-1 
    if abs(c(i))<tiny 
        p(i+1,:)=p(i,:)+ds(i)*t(i,:);
    else
        R=1/c(i);
        p(i+1,:)=p(i,:)+R*n(i,:)-R*n(i+1,:);
    end 
end


if d0>0 
    s=[0;s+d0;s(end)+2*d0];
    t=[t(1,:);t;t(end,:)];
    n=[n(1,:);n;n(end,:)];
    c=[0;c;0];
    p=[p(1,:)-d0*t(1,:);p;p(end,:)+d0*t(end,:)];
else
    c=[c;0];
end
track=[s,p,t,n,c];
end

