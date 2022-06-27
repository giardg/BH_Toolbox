function [M,F,S] = assemblage(coord, connec, bord, numer, bord_ess, xvec, Hk, dHk, P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage
% 
% Assemble les matrices M, F, S pour une iteration du probleme elements
% finis par la methode de Newton
%
% Inputs: - tableaux coord, connec et numer,
%         - vecteurs de points de la fonction Hk, de sa derivee (dHk) et de
%           la puissance (P) aux points de xvec
%         - bord et bord_ess pour mieux gerer F et S eventuellement
%
% Outputs: - M[ndll x ndll] matrice globale du systeme
%          - F[nddl x 1] composant du vecteur de sollicitations
%          - S[nddl x 1] composant du vecteur de sollicitation (pas encore
%            bien implemente)
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_num
global param_phy
deg = param_num.deg;

[nel, ndk] = size(connec);

nddl = size(numer,1)*size(numer,2);

M = zeros( nddl, nddl );
F = zeros( nddl, 1 );
S = zeros( nddl, 1 );

% Construction de la matrice d'adressage (associer les ddl aux noeuds)
adress = numer(connec);

for K=1:nel %boucle sur les elements
    
    % Coordonnees des noeuds geometriques
    intervalle = coord(connec(K,:));
    x1 = intervalle(1);
    x2 = intervalle(2);
    
    % Jacobien pour l'element d'interet
    jac = detTk(intervalle);
    
    Ladress   = adress(K,:); %ligne associee a l'element K dans adress (ndk x 1)
    
    % Initialisation des matrices elementaires
    FK = zeros(ndk,1);
    SK = zeros(ndk,1);
    MK = zeros(ndk,ndk);
    
    [pts,poids] = ptsGauss1d(param_num.nquad); %points d'evaluation de l'integrale et poids
    for k = 1:param_num.nquad 
        
        %Points et poids de gauss
        xi = pts(k);
        w = poids(k);
        xq=(x2-x1)/2*xi+(x1+x2)/2;   %global coordinate of xi
        [phi, dphi] = shape(xi,param_num.deg); %evaluation des fonctions de base aux points de gauss
        
        %Evaluation des fonctions de x par interpolation aux points d'integration
        HkQ = spline(xvec,Hk,xq);
        dHkQ = spline(xvec,dHk,xq);
        PQ = spline(xvec,P,xq);
        
        for i = 1:ndk
            phi_i = phi(i);
            dphi_i = dphi(i)/jac;
            
            % Evaluation des elements de vecteurs elementaires
            FK(i) = FK(i)-w*PQ*phi_i*jac-(param_phy.rho/2)*w*HkQ*dHkQ*dphi_i*jac;
            SK(i) = 0; %Temporaire vecteur de sollicitations 
            
            for j = 1:ndk  
                phi_j = phi(j);
                dphi_j = dphi(j)/jac;
                
                % Evaluation des elements de matrices elementaires
                MK(i,j) = MK(i,j)+(param_phy.rho/2)*w*HkQ*dphi_j*dphi_i*jac;
                MK(i,j) = MK(i,j)+(param_phy.rho/2)*w*dHkQ*phi_j*dphi_i*jac;
            end
        end
    end
    
    % Ajout a la matrice globale selon l'adressage
    M(Ladress,Ladress) = M(Ladress,Ladress) + MK; 
    F(Ladress) = F(Ladress) + FK;
    S(Ladress) = S(Ladress) + SK;

end