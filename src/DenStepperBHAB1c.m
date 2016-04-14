% Better Hybrid AB.
% Integrate using trapezoid method, 
% say Gamma_(n+1) = Gamma_n

function [rho_FT_next] = ...
    DenStepperBHAB1c( Prop, rho_FT, GammaEx_FT,dt)


rho_FT_next = Prop .* (rho_FT) + dt / 2 * ( 1 + Prop) .* GammaEx_FT;

end