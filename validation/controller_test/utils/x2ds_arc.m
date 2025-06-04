function ds=x2ds_arc(x,x0,t,c)
ds=[0 0];
tiny=1e-6;
n=[-t(2),t(1)];
R=1/c;
cen=x0+R*n;

dx=x-cen;
d=norm(dx);
q=atan2(dx(2),dx(1));

if d<tiny
    ds=[0,R];
    return
end

dx0=x0-cen;
q0=atan2(dx0(2),dx0(1));


dq=wrap_angle(q-q0)*sign(c);
dsx=abs(R)*dq;
dsy=(abs(R)-d)*sign(c);
ds=[dsx,dsy];

end

