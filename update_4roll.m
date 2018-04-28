function [Uhat,Shat,newRHS]=update_4roll(Ny,lam,xi,diffconst,t0,tend,savetime)
%
% [U,X,Shat,newRHS]=update_4roll_nt(Ny,lam,xi,nu,t0,tend,savetime,pertid,conf)
%

[grid,params,Shat,newRHS]=get_4roll_inputs(Ny,lam,xi,diffconst);

% unpack parameters

lam = params.lam;
xi = params.xi;
nu = params.nu;
dc = params.diffconst;
dt = params.dt;

Lx = grid.Lx;
Ly = grid.Ly;
xmin = grid.xmin;
ymin = grid.ymin;
Nx = grid.Nx;
Ny = grid.Ny;
dx = grid.dx;


datadir=sprintf('./SOB_4roll/lam%1.1f',lam);

runname    = '4roll_';

fileprefix = sprintf('%s_n%03d_lam%1.2f_dc%d',runname,Ny,lam,dc);

[status,message,messageid] = mkdir(datadir);  % make directory if it is not already present

paramfile  = sprintf('%s/PARAMS_%s.txt',datadir,fileprefix)  % parameter file


% write parameters to file
%
fileID = fopen(paramfile,'w');
fprintf(fileID,'lam = %f\n',lam);
fprintf(fileID,'xi = %f\n',xi);
fprintf(fileID,'nu = %f\n',nu);
fprintf(fileID,'Lx = %f\n',Lx);
fprintf(fileID,'Ly = %f\n',Ly);
fprintf(fileID, ' xmin = %f\n',xmin);
fprintf(fileID, ' ymin = %f\n',ymin);
fprintf(fileID,'dx = %f\n',dx);
fprintf(fileID, 'Ny = %d\n',Ny);
fprintf(fileID,'dt = %0.8f\n',dt);
fprintf(fileID, ' End time = %f\n',tend);

fclose(fileID);


if t0~=0   % load in old data for restart
    fin = sprintf('%s/%s_t%1.2f.mat',datadir_in,fileprefix_in,t0);
    load(fin);
end

frollFhat = forcing2D_4roll(grid);  % compute the constant (in time) forcing

% Start the time stepping

for t=t0:dt:tend
    
    if(isnan(Shat(1,1,1)))  % if the stress blows up stop running
        break
    end
    
    
    if(lam==0)  % Stokes solve
        
        fbhat = frollFhat;   
        
    else  % Compute VE force
        
        veforce_hat = get_veforcehat(Shat,xi,grid);
        
        fbhat = veforce_hat+frollFhat;
        
    end
    
    % Find Velocity
    % Solve Stokes equation with VE stress
    
    Uhat = stokes_solve_fourier(fbhat,grid.Lx,grid.Ly);
    
    % Update the stress tensor if not solving Stokes
    if(lam~=0)
        
        [Shat, newRHS] = update_Shat(Uhat,grid,Shat,nu,dt,lam,newRHS);
        
    end
    
    if( mod(t,savetime)==0  )  % save data when specified
        
        foutw = sprintf('%s/%s_t%1.2f.mat',datadir,fileprefix,t);
        
        save(foutw,'Uhat','Shat','newRHS');
        
        fprintf('time=%1.2f \n',t);
    end
    
end



