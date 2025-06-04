function [s,eflag]=x2ss(x,track,s0,Tol)

s=[0,0];eflag=0;
s0=s0(1);
if nargin <4
    Tol=20;
end

small=1e-3;
N=size(track,1);

test=norm(track(1,2:3)-track(N,2:3))+norm(track(1,4:5)-track(N,4:5));
if test<small
    disp('closed track detected in x2ss')
    s=x2ssc(x,track,s0);
    return
end

indx=interp1(track(:,1),1:N,[s0-Tol,s0+Tol],'nearest','extrap');
a=indx(1)-1;a=max(a,1);
b=indx(2)+1;b=min(b,N);
subtrack=track(a:b,:);

[s,eflag]=x2s0(x,subtrack);

end