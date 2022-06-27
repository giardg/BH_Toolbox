function [phi, dphi, ddphi] = resoudPhi(xvec, H, dH, x_Pjoule, P_Joule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resoud_phi
% 
% Resolution de l'equation en phi: 
% phi' = -(1/H)*sqrt((2/rho)*P_Joule - (H')^2)
%
% Inputs: - xvec: vecteur de coordonnees pour lesquelles les fonctions
%           H et dH sont evaluees
%         - H et dH: resultats du probleme elements finis resolu dans
%           fem_1D
%         - x_Pjoule: vecteur de coordonnees associe a P_Joule
%         - P_joule: vecteur de puissances associees aux pertes par
%           courants de Foucault
%
% Outputs: - phi, dphi et ddphi: vecteur de la phase et de ses derivee
%            (premiere et seconde)
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_num

P_Joule_interp = spline(x_Pjoule, P_Joule,xvec);
h = (xvec(2)-xvec(1)); %espacement entre les points
nmax = length(xvec)-1; %nombre de points pour la resolution
x0 = xvec(1); %point initial

% Resolution par Runge-Kutta d'ordre 4
[x, phi] = solve_rk4(x0, param_num.phi0, h, nmax, H, dH, P_Joule_interp);

[dphi, ddphi] = diffO2(x,phi);