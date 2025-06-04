function [s,e]=x2s0(x,track)
s=[0 0];
e=0;
N=size(track,1); 
one=ones(N,1);
dx=one*x-track(1:N,2:3);
dxt=sum(dx.*track(1:N,4:5),2);
test=(dxt>=0); 
temp=diff(test);
I=find(temp==-1);
if isempty(I)
    
    if test(1)==0
        i0=1;e=-1;
    else
        i0=N;e=1;
    end
    disp('warning - search error in x2s0 in UTILS, [x,e]= '),disp([x,e])
else
    X0=track(I,2:3);X1=track(I+1,2:3);dX=X1-X0;
    D=hypot(dX(:,1),dX(:,2));
    Tv=dX./(D*[1 1]);
    Nv=[-Tv(:,2),Tv(:,1)];
    offset=sum(Nv.*dx(I,:),2);
    [~,j]=min(abs(offset));
    i0=I(j);
end

ds=x2ds(x,track(i0,:));
s=[track(i0,1)+ds(1),ds(2)];

end





