% Better Hybrid AB2 with prefactor
% Integrate using trapezoid method, say 
% Gamma_(n+1) = 2 * Gamma_n - Gamma_(n-1)

function [rho_FT_next] = DenStepperBHAB2cPf( ...
    Prop, rho_FT, GammaEx_FT, GammaEx_FTprev, NlPf, NlPfprev)


rho_FT_next = Prop .* (rho_FT) + NlPf .* GammaEx_FT - NlPfprev .* GammaEx_FTprev  ;

end