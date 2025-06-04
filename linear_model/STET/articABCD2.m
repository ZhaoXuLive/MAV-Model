function [A,B1,B2,C,D]=articABCD2(PARSu,PARSa,U)
% articulated vehicle linear model
% Linear Model of n units with N axles
% n and N determined from vehicle data
% states: body slip angle for unit 1, yaw rates, articulation angles (2n)
%   = [beta_1,r1,r2,...,rn,theta_1,theta_2,...theta_n-1]
% inputs: steer angles at each axle (N inputs) .. input matrix B1
%     AND yaw moments at each unit (n inputs)  .. input matrix B2
% 
% outputs: all states
%
% M= n-vector of masses (row or col, n elements, mass in kg)
% J= n-vector of CG yaw moments of inertia (ditto, kg m^2)
% L= n-element cell array of axle locations {L1,L2,...}
%   Li=[x1,x2,...] x coords relative to the mass centre (m) on each unit
%  (so the units may have different numbers of axles)
% Ca=cell array of axle cornering stiffnesses - same structure as L
%   units: N/rad
% Xf= n-vector of front joint x coords rel to CG (front of unit 1 is not used)
% Xr= n-vector of rear joints (similar)
% U=vehicle speed (m/s) ... assumed constant

% PARSu: [M J L xg Bj]
% PARSa: [xa u Ca DRIVE BRAKE STEER AUTOSTEER]
% xg: 各单元质心到后面的铰接点的距离 均为正 xr = -xg 均为负
% xa: 各轴中心到后面的铰接点的距离  
% xf: 各单元质心到前面铰接点的距离 均为正 xf = L - xg

n=size(PARSu,1); %number of units
N=size(PARSa,1); %number of axles

% convert from PARSu and PARSa to the above variables
M=PARSu(:,1)';
J=PARSu(:,2)';
Xr=-PARSu(:,4)';
Xf=PARSu(:,3)'+Xr;
% Bj=PARSu(:,5)'; 
L=cell(1,n);Ca=cell(1,n);
for k=1:n
    % ind: 第k单元对应的轴，10*1 logical
    ind=PARSa(:,2)==k;
    % L1:  各轴中心到单元质心的距离，前为正
    L1=PARSa(ind,1)'+Xr(k);
    % 对应轴侧偏刚度
    Ca1=PARSa(ind,3);
    L{k}=L1;Ca{k}=Ca1;
end

nx=3*n-1; %number of states
nv=4*n-2; %number of variables (including joint forces)

% dimensions of primary matrices
S=zeros(nv,nv);
P=zeros(nv,nx);
BY=zeros(nv,N);
K1=zeros(N,nx);
K2=zeros(N,N);
%new
BM=zeros(nv,n); %input matrix for external yaw moments

%populate S,P,B equation-by-equation

%first define end-index of sub-matrix blocks
%3n-1个状态，前n个是质心侧偏角，中间n个是横摆角速度，最后n-1个是铰接角
b1=n;    %last row of beta_dot
b2=2*n;  %last row of r_dot
b3=3*n-1; %last row of theta_dot

%also need number of 'preceding' axles at each unit
%每个单元前面所有的单元共有多少个轴
pax=zeros(1,n); %preceding axle numbers
for i=2:n
    pax(i)=pax(i-1)+length(L{i-1});
end

%equation (1)
for i=1:n 
    S(i,i)=M(i)*U;
    P(i,b1+i)=M(i)*U;
    if i>1,S(i,b3+i-1)=-1;end %R_(i-1) term in equ(1)
    if i<n,S(i,b3+i)=1;end %R_i term 
    for j=1:length(L{i})
       BY(i,pax(i)+j)=1; %Yij terms
    end
end

%equation (2)
for i=1:n 
    S(b1+i,b1+i)=J(i);
    if i>1,S(b1+i,b3+i-1)=-Xf(i);end %R_(i-1) term in equ(1)
    if i<n,S(b1+i,b3+i)=Xr(i);end %R_i term 
    for j=1:length(L{i})
       BY(b1+i,pax(i)+j)=L{i}(j); %Mij terms
    end
    BM(b1+i)=1; %direct yaw moment term
end

%equation (3)
for i=1:n-1 %joint index
    S(b2+i,b2+i)=1; %theta-dot term
    P(b2+i,b1+i)=-1; %r_i in notes = unit in front of the joint
    P(b2+i,b1+i+1)=1; %r_(i+1)
end

%equation (4) 
for i=1:n-1 %joint index, terms in order as per the notes
    S(b3+i,i+1)=U;
    S(b3+i,i)=-U; 
    S(b3+i,b1+i+1)=Xf(i+1);
    S(b3+i,b1+i)=-Xr(i);
    S(b3+i,b2+i)=-U; 
end

% K1 matrix,Y=K1*X for each axle force
for i=1:n
    for j=1:length(L{i})
        ia=pax(i)+j; %'flat' axle index
        K1(ia,i)=U;
        K1(ia,b1+i)=L{i}(j);
        K1(ia,:)=-K1(ia,:)*Ca{i}(j)/U;
    end
end

% K2 matrix, go axle by axle
for i=1:n
    for j=1:length(L{i})
        ia=pax(i)+j; %'flat' axle index
        K2(ia,ia)=Ca{i}(j);
    end
end

% disp(L);
% disp(pax);
% disp(BM);

AA=S\(BY*K1-P);
A=AA(1:nx,:);

BB=S\BY*K2;
B1=BB(1:nx,:);
BB=S\BM;
B2=BB(1:nx,:);

% BM ？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？

% eliminate the constraints to use a smaller set of independent states xa = beta_1,r,theta [2*n states]
% constraint is of the form Qx=0
Q=S(b3+1:end,1:nx); %same terms as in equation (4)
temp=eye(nx);
T1=temp([1,b1+1:nx],:); %selects the independent states
T2=temp(2:b1,:); %selects the states to be eliminated
T=[T1;T2];
Tin=inv(T);Tin1=Tin(:,1:2*n);Tin2=Tin(:,2*n+1:end);

E=-inv(Q*Tin2)*(Q*Tin1); 
Umat=Tin1+Tin2*E;
Am=T1*A*Umat;
B1m=T1*B1;
B2m=T1*B2;

%outputs = all states
Cm=eye(2*n);
Dm=zeros(2*n,N+n); %dimensioned for B=[B1 B2], steer and direct yaw moment

%finally overwrite A etc. for function output
A=Am;B1=B1m;B2=B2m;C=Cm;D=Dm;
end









