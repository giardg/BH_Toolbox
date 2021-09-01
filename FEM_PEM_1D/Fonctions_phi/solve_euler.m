function  [x, phi] = solve_euler(x0, phi0, h, nmax, H, dH, P_Joule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve_euler
% 
% Applique la methode d'Euler explicite d'ordre 4 pour resoudre l'equation
% en phi: phi' = -(1/H)*sqrt((2/rho)*P_Joule - (H')^2)
%
% Inputs: - x0: coordonnee initiale
%         - phi0: phase initiale
%         - h: espacement entre les coordonnees
%         - nmax: nombre d'espacements
%         - H et dH: resultats du probleme elements finis resolu dans
%           fem_1D
%         - P_joule: vecteur de puissances associees aux pertes par
%           courants de Foucault (defini sur les memes points que H et dH)
%
% Outputs: - x: vecteur de points sur lesquels phi est evalue
%          - phi: phase (solution de l'equation en phi)
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_phy

nbeq = size(phi0,2);
phi0 = phi0';
k1 = zeros(nbeq,1);
phi = zeros(1,nbeq);
n = 1;
x(1) = x0;
phi(1,:) = phi0';

while(n <= nmax)
    k1 = (-1/H(n))*sqrt(max((2/param_phy.rho)*P_Joule(n)-dH(n)^2,0));
    phi0 = phi0 + h * k1;
    phi(n+1,:) = phi0';
    x(n+1) = x(n) + h;
    n = n+1;
end


