function [GridObj] = GridMakerPBCxk(Nx,Ny,Nm,Lx,Ly)

% Make pos grid spacings
dx   = Lx/Nx;              
dy   = Ly/Ny;              
dphi = 2*pi/Nm;
% Make vectors and grids
x                = ( 0: dx : Lx - dx );
% x                = ( -Lx/2 : dx: Lx/2 - dx); 
y                = ( 0: dy : Ly - dy );
% y                = ( -Ly/2 : dy: Ly/2 - dy);
phi              = ( 0: dphi: (2*pi - dphi) );
[y2D,x2D]        = meshgrid(y,x);                          
[y3D,x3D,phi3D]  = meshgrid(y,x,phi);            

% Make k-space spacings
dkx          = 2*pi/Lx;
dky          = 2*pi/Ly;
dkm          = 1;
% Make k vectors and grids
kx_max           = pi/dx;            %Maximum spatial k-vector allowed by grid
ky_max           = pi/dy;            %Maximum spatial k-vector allowed by grid
km_max           = pi / dphi;        %Maximum angular k-vector
kx               = ( -kx_max: dkx: (kx_max - dkx) );
ky               = ( -ky_max: dky: (ky_max - dky) );
km               = ( -km_max: dkm: (km_max - dkm) );
[ky2D,kx2D]      = meshgrid(ky,kx);                  %2D k-vector grid
[ky3D,kx3D,km3D] = meshgrid(ky,kx,km);          %3D k-vector grid

%Put it all in an object
GridObj = struct('x',x,'y',y,'phi',phi,'dx',dx,'dy',dy,'dphi',dphi,'x2D',x2D,'y2D',y2D,'x3D',x3D,'y3D',y3D, 'phi3D',phi3D,...
                     'kx',kx,'ky',ky,'km',km,'kx2D',kx2D,'ky2D',ky2D,'kx3D',kx3D,'ky3D',ky3D, 'km3D',km3D );

end
