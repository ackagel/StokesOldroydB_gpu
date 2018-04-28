
addpath('./src');

Lx = 1;
Ly = 1;
xmin=-.5;
ymin=-.5;
lam = .3;
tin = 10;
dc = 0;
ns = [16 32 64 128];

for n=ns
    fname = sprintf('./SOB_4roll/lam%1.1f/4roll__n%03d_lam%1.2f_dc%d_t%1.2f.mat',lam,n,lam,dc,tin);
    load(fname);
    S=real(ifft2(Shat));
    
    xx=[-.5:1/n:.5-1/n];
    plot(xx,S(n/2+1,:,1),'-o')
    hold all
end


figure

gUh = matrix_derivative_fourier(Uhat,Lx,Ly);
gU = real(ifft2(gUh));
vorticity = gU(:,:,3)-gU(:,:,2);
dx=1/n;
[x,y]=ndgrid(xmin:dx:xmin+Lx-dx,ymin:dx:ymin+Ly-dx);
contourf(x,y,vorticity,25);

