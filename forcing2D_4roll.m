function Fhat=forcing2D_4roll(grid)

% This creates the fourier transform of the force field which is periodic
% in x and y directions


Lx = grid.Lx;
Ly = grid.Ly;
xmin = grid.xmin;
ymin = grid.ymin;

dx = grid.dx;


[x,y]=ndgrid(xmin:dx:xmin+Lx-dx,ymin:dx:ymin+Ly-dx);

F(:,:,1) =  2*sin(2*pi/Lx*x).*cos(2*pi/Ly*y);
F(:,:,2) = -2*cos(2*pi/Lx*x).*sin(2*pi/Ly*y);

Fhat = fft2(F);



