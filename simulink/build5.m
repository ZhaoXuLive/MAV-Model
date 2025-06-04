function track=build5(L,Q,p0,t0)


tiny=1e-6;
if ~exist('p0','var'), p0=[];end
if isempty(p0), p0=[0 0]; end
if ~exist('t0','var'), t0=[];end
if isempty(t0), t0=[1 0]; end


ic=find(abs(Q)>tiny); 
C=zeros(size(L));
C(ic)=sign(Q(ic))./L(ic); 


dS=L;
dS(ic)=dS(ic).*abs(Q(ic));

C=[C,0];
S=[0,cumsum(dS)];
track=build4(p0,t0,S',C');