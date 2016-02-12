% Function: dRhoInterCalcVcId.m

% Description: Wrapper function to calcuate the interaction contribution
% to the PDE in k-space. This calls functions that calculate the excess 
% chemical potential calculated by the virial expansion. It can call
% multiple expansions, but it should not be trusted above 2nd order. From
% excess chemical potential, it calculates the flux, then take the
% divergence of it, assuming isotropic diffusion.
% 
% ISOTROPIC Diffusion
% 
% Called by:  HR2DrotDenEvolverFTBodyIDCube.m - main body of isotropic
% diffusion

% Calls: 
% MuExCalcVc2Id(rho_FT, Fm_FT, ParamObj)- Calculates excess Chem pot. Using
% 2nd virial approx.

function [NegDivFluxExcess_FT] = ...
       dRhoIntCalc2D(rho,rho_FT,V_FT,ParamObj,GridObj,Mob)
%%%%%%%%%%%%%%%%%%%Hard rod %%%%%%%%%%%%%%%%


%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile
% keyboard
%Now includes the correct scale

[MuEx_FT] = ParamObj.Lx * ParamObj.Ly / ...
(ParamObj.Nx * ParamObj.Ny) * rho_FT * V_FT;

% MuEx_FT    = MuEx2_FT+MuEx3_FT;
%Takes its derivative in k-space
dMuEx_dx_FT   =  sqrt(-1) .* GridObj.kx2D .*  MuEx_FT;
dMuEx_dy_FT   =  sqrt(-1) .* GridObj.ky2D .*  MuEx_FT;

%Excess chemical potential derivative in real space
%Mayer function derivative in real-space
dMuEx_dx   =  real(ifftn(ifftshift(dMuEx_dx_FT)));
dMuEx_dy   =  real(ifftn(ifftshift(dMuEx_dy_FT)));

%Do the hard disk interaction portion of the PDE in real space then FT it
% Isolate the seperate parts and call them some arbitrary function. We
% will Fourier transform these functions to solve this in Fourier space
%
% Take the divergence of the product of functions. Call these products
% random variables

jx = - Mob .* rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - Mob .* rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));

% Calculate the - divergence of the interaction flux
NegDivFluxExcess_FT = - sqrt(-1) .* ( GridObj.kx2D .* Jx_FT + ...
    GridObj.ky2D .* Jy_FT  );

