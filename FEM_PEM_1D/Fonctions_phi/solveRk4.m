function  [x, phi] = solveRk4(x0, phi0, h, nmax, H, dH, P_Joule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve_rk4
% 
% Applique la methode de Runge-Kutta d'ordre 4 pour resoudre l'equation en
% phi: phi' = -(1/H)*sqrt((2/rho)*P_Joule - (H')^2)
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
k2 = zeros(nbeq,1);
k3 = zeros(nbeq,1);
k4 = zeros(nbeq,1);
phi = zeros(1,nbeq);
n = 1;
x(1) = x0;
phi(1,:) = phi0';

while(n <= nmax)
    
    % Comme on connait les fonctions P_Joule, H et dH seulement pour des
    % points precis, je ne suis pas certain qu'il soit valide d'evaluer
    % les fonctions en x+h/2 seulement en prenant la moyenne de f(n) et
    % f(n+1)... (voir k2 et k3)
    k1 = h * (-1/H(n))*(sqrt(max((2/param_phy.rho)*P_Joule(n)-dH(n)^2,0)));
    k2 = h * (-2/(H(n)+H(n+1)))*(sqrt(max((2/param_phy.rho)*((P_Joule(n)+P_Joule(n+1))/2)-((dH(n)+dH(n+1))/2)^2,0)));
    k3 = k2;
    k4 = h * (-1/H(n+1))*(sqrt(max((2/param_phy.rho)*P_Joule(n+1)-dH(n+1)^2,0)));
    
    phi0 = phi0 + (1/6) * (k1+2*(k2+k3) + k4);
    phi(n+1,:) = phi0';
    x(n+1) = x(n) + h;
    n = n+1;
end
    
    