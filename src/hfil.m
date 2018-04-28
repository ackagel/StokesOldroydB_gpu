function Ffil=hfil(F,Lx,Ly)

% Compute the Hou-filter to cut-off high fourier modes

z=size(F);

if ndims(F)<3, k=1; else
   k=z(3);end

nx=z(1);  % try again
ny=z(2);

kx = (1/Lx)*(-nx/2+1:1:nx/2);
ky = (1/Ly)*(-ny/2+1:1:ny/2);

[k1, k2]=meshgrid(kx,ky);


cutoff=exp(-36*(sqrt((k1/(nx/2)).^2+(k2/(ny/2)).^2)).^36);

for i=1:k
    fftcutoff(:,:,i)=fftshift(cutoff);
end

Ffil=fftcutoff.*F;