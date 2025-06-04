function []=artic_shape_01(q,PARSu,color)
%  []=ARTshape(q,scale,color)
%     plot the shape of the ART using configuration variables:
%    [q1,q2]=x,y coords of the front (virtual) hitch point of unit_1
%     q3=yaw angle of unit 1, q4 for unit 2, q5 for unit 3 etc.
%     WORKS FOR N UNITS, N>=1 and with variable lengths, showing axle
%     locations
L=PARSu(:,3);Nu=size(PARSu,1);
W=2.5; %nominal width
d=0.5; %distance between joint and body shape


if ~exist('color','var'),color=[];end
if isempty(color)
color=[1 0 0]; %body color (RGB)
end


%standard body shape is a unit square, clockwise from top left, cenred at
%the origin
car0=[
    -1,1,1,-1,-1
     1,1,-1,-1,1 ]/2;


%first unit
scale=[L(1)-d,W];
theta=q(3);
rot=[
   cos(theta) -sin(theta)
   sin(theta)  cos(theta)];
mat=rot*diag(scale);
ux=rot*[1 0]';%uy=rot*[0 1]';
R1=[q(1) q(2)]'; %coords of reference point on first unit (front virtual hitch)
%find new reference (coupling) point at rear of first unit, R2
R2=R1-L(1)*ux;
Rc=(R1+R2)/2; %centre of box
car=Rc*ones(1,size(car0,2),1)+mat*car0;
patch(car(1,:),car(2,:),color,'EdgeColor',[0 0 0]);
hold on
plot(R1(1),R1(2),'ko')
plot(R2(1),R2(2),'ko')

%repeat for the following units
%slight difference ... shorten the body by 2*d for interior units
%the joint

for i=2:Nu
    R1=R2;
    scale=[L(i)-d,W];
    theta=q(i+2);
    rot=[
        cos(theta) -sin(theta)
        sin(theta)  cos(theta)];
    mat=rot*diag(scale);
    ux=rot*[1 0]';
    R2=R1-L(i)*ux;
    Rc=(R1+R2)/2; %centre of box
    car=Rc*ones(1,size(car0,2),1)+mat*car0;
    patch(car(1,:),car(2,:),color,'EdgeColor',[0 0 0]);
    plot(R2(1),R2(2),'ko')

end





