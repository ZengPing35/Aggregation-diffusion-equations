[x,y,z] =cylinder(pi:-pi/5:0,10)
colormap(lines)
subplot(1,3,1)
surf(x,y,z);
shading flat
subplot(1,3,2)
surf(x,y,z);
shading interp
subplot(1,3,3)
surf(x,y,z)