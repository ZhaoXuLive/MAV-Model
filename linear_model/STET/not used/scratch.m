% convert PARSu and PARSa to/from parameters for LinearN

%% LinearN to PARSu / PARSa 
Nu=length(M);
PARSu=zeros(Nu,5);
PARSu(:,1)=M;
PARSu(:,2)=J;
PARSu(:,3)=Xf-Xr; %unit lengths
PARSu(:,4)=-Xr; %CG location relative to rear joint
PARSu(:,5)=Bj;

Na=sum(cellfun(@numel,L));
PARSa=zeros(Na,7);
PARSa(:,1)=[L1-Xr(1),L2-Xr(2),L3-Xr(3)];
PARSa(:,2)=[1 1 2 2 3 3];
PARSa(:,3)=[Ca{:}];
PARSa(:,4)=[1 0 0 0 0 1];
PARSa(:,5)=ones(Na,1);
PARSa(1,6)=1;


%% PARSu / PARSa to LinearN

M=PARSu(:,1);
J=PARSu(:,2);
Xr=-PARSu(:,4);
Xf=PARSu(:,3)+Xr;
Bj=PARSu(:,5);


PARSa(:,1)=[L1-Xr(1),L2-Xr(2),L3-Xr(3)];
PARSa(:,2)=[1 1 1 2 2 3 3 3];
PARSa(:,3)=[Ca1,Ca2,Ca3];
PARSa(:,4)=[0 1 1 0 0 0 0 0];
PARSa(:,5)=ones(Na,1);
PARSa(1,6)=1;