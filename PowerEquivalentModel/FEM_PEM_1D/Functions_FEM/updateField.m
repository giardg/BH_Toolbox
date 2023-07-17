function [Hk_new, dHk_new, ddHk_new, xvec_new, delta_tot, flag] = updateField(coord, connec, numer, xvec, Hk, U)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateField
% 
% Met a jour le champ et sa derivee (Hk et dHk) selon la solution du
% probleme elements finis
%
% Inputs: - tableaux coord, connec et numer,
%         - xvec: vecteur initial des coordonnees en x
%         - vecteurs de points de la fonction Hk et de sa derivee (dHk) aux
%           points de xvec
%         - U: coefficients qui permettent de determiner les corrections
%           delta (trouve dans resoud_systeme)
%
% Outputs: - vecteurs mis a jour (Hk_new, dHk_new, xvec_new)
%          - delta_tot: corrections aux points de xvec
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag = false;

global param_num
[nel, ndk] = size(connec);

adress = numer(connec);
nx=param_num.deg+1; %nombre de points par element ou on evalue H
xvec_new = [];
Hk_new = [];
delta_tot = [];

for K=1:nel
    for i=1:ndk
        c=adress(K,i);
        uL(K,i)=U(c);
    end
end

for K=1:nel %boucle sur les elements
    
    % Coordonnees des noeuds geometriques
    intervalle = coord(connec(K,:));
    x1 = intervalle(1);
    x2 = intervalle(2);
    
    % Jacobien pour l'element d'interet
    jac = detTk(intervalle);
    
    delta=zeros(nx,1);
    ddelta=zeros(nx,1);
    dddelta=zeros(nx,1);
    
    x=(x1:(x2-x1)/(nx-1):x2);
    xi=(-1:2/(nx-1):1);
    
    for j=1:nx
        [phi, dphi] = shape(xi(j),param_num.deg);
        for i=1:ndk
            delta(j)=delta(j)+phi(i)*uL(K,i);
        end
    end
    
    % Ajout dans les donnees globales
    if K == 1
        xvec_new = [xvec_new ; x'];
        delta_tot = [delta_tot;delta];
    else
        xvec_new = [xvec_new ; x(2:end)'];
        delta_tot = [delta_tot;delta(2:end)];
    end
    
end

Hk_interp = spline(xvec,Hk,xvec_new);

% Amortissement
if param_num.amortissement < 1
    while any(delta_tot+Hk_interp<0)
        delta_tot = param_num.amortissement*delta_tot;
        if all(abs(delta_tot)./Hk(1) < 1e-3)
            flag = true; %arret manuel de l'amortissement
            break
        end
    end
end

if ~flag
    Hk_new = Hk_interp+delta_tot;
else
    Hk_new = Hk_interp;
end

% Derivee 1re et 2nd de H
[dHk_new, ddHk_new] = diffO2(xvec_new,Hk_new);
