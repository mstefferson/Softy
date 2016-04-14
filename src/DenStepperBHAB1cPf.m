% Better Hybrid AB.
% Integrate using trapezoid method, 
% say Gamma_(n+1) = Gamma_n

function [rho_FT_next] = ...
    DenStepperBHAB1cPf( Prop, rho_FT, GammaEx_FT, NlPf)


rho_FT_next = Prop .* (rho_FT) + NlPf .* GammaEx_FT;

end