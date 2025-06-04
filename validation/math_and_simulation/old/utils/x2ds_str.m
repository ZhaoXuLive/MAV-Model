function ds=x2ds_str(x,x0,t)
ds=[0 0];
dx=(x-x0);
n=[-t(2),t(1)];
ds1=dot(dx,t);
ds2=dot(dx,n);
ds=[ds1,ds2];
end
