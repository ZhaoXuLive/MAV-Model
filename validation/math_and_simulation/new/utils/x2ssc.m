function s=x2ssc(x,track,s0)

small=1e-3;
s=[0,0];
N=size(track,1);
smax=track(N,1);
s0=s0(1);

one=ones(N,1);
dx=one*x-track(1:N,2:3);
dxt=sum(dx.*track(1:N,4:5),2);
test=(dxt>=0); 
temp=diff(test);
I=find(temp==-1);

dS=wrap_cts(s0-track(I,1),smax);
[~,j0]=min(dS);
i0=I(j0(1));

ds=x2ds(x,track(i0,:));
s=[track(i0,1)+ds(1),ds(2)];

end