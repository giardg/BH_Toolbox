function H = foucaultSolve( PQ, Mu )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FoucaultSolve.m
% 
% Resoud numeriquement l'equation: rho*H'' = j*omega*mu*H a l'aide de la
% permeabilite equivalente obtenue par le Power Equivalent Model (inspire
% du code de Maxime Tousignant)
%
% Inputs: - PQ: structure avec les parametres physiques (rho, omega, dx,
%               nx)
%         - Mu: donnees de permeabilite equivalente (complexe) pour chaque
%               points du domaine
%
% Outputs: - H: amplitude du champ magnetique (solution de l'equation
%               differentielle)
% 
% Aout 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = PQ.nx;

%% Construction de la matrice
A1 = zeros(N);

% terme de derivee seconde
A1(1,1:2) = [2,0];
for i = 2:N-1
    A1(i,i-1:i+1) = [1,-2,1];
end
A1(N,(N-2):N) = [1,-4,3];
A1 = PQ.rho/PQ.dx^2*A1;

% terme de la permeabilite
A2 = diag(-1i*PQ.omega*Mu);

% combinaison des deux matrices
A = A1+A2;

%% terme source
b = zeros(N,1);
b(1) = PQ.rho*(2/PQ.dx^2)*PQ.Hmax;

%% resolution
H = A\b;

    
