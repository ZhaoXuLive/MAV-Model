% Q_ART_par.m
% data from CRRC Qingdao
% FORMAT
%    PARSu=[M J L xg Bj] (Nu-by-5) = mass, moment of inertia, length, cg
%    location measured from rear joint, rear joint damping (at rear joint)
%
%    PARSa=[xa u Ca DRIVE BRAKE STEER AUTOSTEER]  (Na-by-7) =
%       col 1 = axle location from rear joint
%       col 2   unit number
%       col 3 = cornering stiffness
%       col 4 = drive factor (normally 1 or 0, 1 indicated drive axle)
%       col 5 = brake factor (normally all ones)
%       col 6 = manual steer factor (here 1 for front two axles, zero on others)
 


% axle-related parameters (see above) - geometry from drawing
% use stated body lengths to position virtual front and rear joints at end
% of body structures => 2.5 m from nearest axle 
Na=8; %number of axles
str0=zeros(Na,1); %sets dimension in simulink diagram
Ca=400e3; %rough estimate, and assume all axles have equal cornering stiffness
PARSa=[
    8.9  1  Ca  1  1  1
    7.4  1  Ca  1  1  1
    1.5  1  Ca  0  1  0
    7.9  2  Ca  0  1  0
    1.5  2  Ca  0  1  0
    9.9  3  Ca  0  1  0
    4.0  3  Ca  1  1  0
    2.5  3  Ca  1  1  0];

% unit-related parameters
if 1 %unloaded
    M=[11.142,7.651,11.142]'*1000 ; %pure estimate for unit 2
    J=[143.57,37.8,143.57]'*1000; %ditto
else %fully loaded
    M=[16.910,13.470,16.910]'*1000 ; %pure estimate for unit 2
    J=[225.57,63.071,225.57]'*1000; %ditto
end
TOTALMASS=sum (M,'all'); 

Rgyr=(J./M).^0.5; %radius of gyration, test reasonableness

Nu=3;
q0=zeros(Nu+2,1); %number of generalized coordinates (for simulink)

%unit lengths (from drawing) including (unit 1) 2.5 m from axle 1 to front
%of body = virtual front 'joint'
L=[11.4,9.4,11.4]';

%CG locations relative to rear joint
Xg=zeros(Nu,1);
Xg(1)=5.538;
Xg(2)=L(2)/2; %assume geometric centre
Xg(3)=L(3)-Xg(1);

%assemble into parameter matrix
PARSu=zeros(Nu,5);
PARSu(:,1)=M;
PARSu(:,2)=J;
PARSu(:,3)=L; 
PARSu(:,4)=Xg;
% column 5, joint damping left at zero










