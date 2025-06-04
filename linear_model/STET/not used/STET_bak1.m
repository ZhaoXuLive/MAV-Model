% STability Evaluation Tool (STET.m)
% Use articABCD to build state space (linear) model of ART
clear
close all
g=9.81;deg=pi/180;kph=1/3.6;
% data (ART, initial estimates)
M=[1 1 1]*10000; 
J=[1 1 1]*6.84e4;

%axle locations rel to mass centres
L1=[2.82,-3.68];
L2=[3.0,-3.0];
L3=[3.68,-2.82];
L={L1,L2,L3};

%joint locations (front, rear)
Xf=[4.84 5.02 5.7];
Xr=[-5.7 -5.02 -4.84];

%cornering stiffnesses (axle total, N/rad, all assumed equal)
% Ca=3.5e5;
% axle loads (used in the Ca estimate)
AL1=M(1)*g*[-L1(2) L1(1)]/(L1(1)-L1(2));
AL2=M(2)*g*[-L2(2) L2(1)]/(L2(1)-L2(2));
AL3=M(3)*g*[-L3(2) L3(1)]/(L3(1)-L3(2));

Ca=3.5e5;
Ca=[Ca Ca];Ca={Ca,Ca,Ca}; %same for all units


%speed
U=50*kph;

[A,B,C,D]=articABCD(M,J,L,Ca,Xf,Xr,U);
[eigvec,eigval]=eig(A);
eigval=diag(eigval);


% animation of least stable mode - remove real part of eigenvector for animation if motion is sinusoidal
% the time constant is the measure of stability of instability (time taken
% to reduce or increase amplitude by a factor of e = 2.72 (approx)

[s,k]=max(real(eigval));
if s>0 
    disp('unstable dynamics')
end
tc=1/abs(s); tc=round(100*tc)/100;
disp(['time constant = ',num2str(tc),' s'])
disp([])
disp(['mode number = ',num2str(k)])

close all
tiny=1e-5;
Nu=length(M); %number of units
Q=zeros(1,Nu+2);

%construct PARSu
PARSu=zeros(Nu,5);
PARSu(:,1)=M;
PARSu(:,2)=J;
PARSu(:,3)=Xf-Xr; %unit lengths
PARSu(:,4)=-Xr; %CG location relative to rear joint


t=(0:0.1:30)';%column vector of times
w=imag(eigval(k));

X=eigvec(:,k)'/5;%row vector
% time history has the form: real(v*T)

if abs(w)>tiny
    T=exp(1j*w*t);
else
    tau=5; %arbitrary time constant for visualization
    T=exp(-t/tau); 
end


%convert eigenvector X to standard q vector

a=Xf(1); %distance from CG1 to front joint
if abs(w)>tiny
    Q(2)=(U*X(1)+a*X(2))/(1j*w); %integration of lateral velocity at joint
    Q(3)=X(2)/(1j*w); %integration of front yaw rate.
else
    Q(2)=(U*X(1)+a*X(2))*(-tau); %integration of lateral velocity at joint
    Q(3)=X(2)*(-tau); %integration of front yaw rate.
end
for i=2:Nu %use articulation angles to find other yaw angles
    Q(2+i)=Q(1+i)+X(Nu+i); %e.g. i=2,nu=3: Q(4)=Q(3)+X(5)
end
q=real(T*Q);
q(:,1)=U*t; %longitudinal distance


% animate_artic_01(q,nstep,PARSu,TRACK)
animate_artic_01(q,1,PARSu)








