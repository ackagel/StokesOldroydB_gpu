
addpath('./src');


Nys =[16 32 64 128];           % number of grid points in y-direction

lam = 0.3;        % relaxation time
diffconst = 0;   % diffusion =  (diffconst*dx)^2
xi=.5;           % polymer viscosity/solvent viscosity
t0 = 0;          % initial time
tend = 10.0;     % end time
savetime=1;      % when to save




for Ny=Nys
    tic
    [Uhat,Shat,newRHS]=update_4roll(Ny,lam,xi,diffconst,t0,tend,savetime);
    comp_time = toc;
    fprintf('total computation time  = %g \n',comp_time);
end


