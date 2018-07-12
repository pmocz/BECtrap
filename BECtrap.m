close all
clear all
clc

%% Simulate a BEC with SI in a trap
% Philip Mocz (2018)

% units hbar = m = 1
% [i d_t + nabla^2/2 + V_trap - lambda*|psi|^2]  = 0

%% parameters
Lx = 10;
Ly = 80;
Lz = 10;

Nx = 32;
Ny = Nx * Ly/Lx;   assert(Ny==round(Ny));
Nz = Nx * Lz/Lx;   assert(Nz==round(Nz));

lambda = 0.3;   % attractive SI strength

t_final = 10;

%% some more setup
dx = Lx/Nx;
xlin = linspace(0.5*dx,Lx-0.5*dx,Nx);
ylin = linspace(0.5*dx,Ly-0.5*dx,Ny);
zlin = linspace(0.5*dx,Lz-0.5*dx,Nz);
[x, y, z] = meshgrid(xlin, ylin, zlin);

% potential
V = 10*0.5*(((x-Lx/2)/Lx).^2+((y-Ly/2)/Ly).^2+((z-Lz/2)/Lz).^2);
clear x;
clear y;
clear z;

t = 0;
dt = min( dx^2/6, 1/max(abs(V(:))));

%% fourier space variables

fftw('planner','measure');
kxlin = ((-Nx/2:Nx/2-1)') * (2*pi/Lx);
kylin = ((-Ny/2:Ny/2-1)') * (2*pi/Ly);
kzlin = ((-Nz/2:Nz/2-1)') * (2*pi/Lz);

[kx, ky, kz] = meshgrid(kxlin, kylin, kzlin);
kSq = fftshift(kx.^2 + ky.^2 + kz.^2);
clear kx;
clear ky;
clear kz;

%% initial contition

psi = exp(-V/(2*0.5^2));
%psi = exp(-1.i.*8*y*2*pi/Ly) .* exp(-1.i.*x*2*pi/Lx);


%% time evolution
fh = figure('position',[10 10 200 200*Ly/Lx]);
while t < t_final
    t
    
    % kick
    psi = exp(-1.i * dt/2 * (V-lambda*abs(psi).^2)).*psi;
    
    % drift
    psi = ifftn(exp(dt * (-1.i*kSq/2)).*fftn(psi));
    
    % kick
    psi = exp(-1.i * dt/2 * (V-lambda*abs(psi).^2)).*psi;
    
    t = t + dt;
    
    % make plots / diagnostics
    %imagesc([0 Lx], [0 Ly], log10(abs(psi(:,:,end/2+1)).^2)); % slice -log
    %imagesc([0 Lx], [0 Ly], real(psi(:,:,end/2))); % real part 
    imagesc([0 Lx], [0 Ly], log10(mean(abs(psi).^2,3))); % projection -log
    caxis([-4,1])
    pbaspect([1 Ly/Lx 1])
    axis off
    pause(0.0001)
    
end

