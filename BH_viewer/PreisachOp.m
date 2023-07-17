function M = PreisachOp(E,H)
% -------------------------------------------------------------------------
% Implementation du modele d'hysteresis scalaire de Preisach
% Maxime Tousignant 2019
%
% Entrees :
%   H       Champ magnetique au temps t
%   E(a,b)  Fonction d'Everett
%
% Sortie  :
%   M       Magnetisation au temps t
% -------------------------------------------------------------------------

% Declaration des variables internes
persistent Hp Mpartiel n alpha beta

% Initialisation des variables internes lors du premier appel
if isempty(Hp)
    Hp       = 0.0;     % Champ magnetique au pas de temps precedent
    Mpartiel = 0.0;     % Somme partielle de la magnetisation
    n        = 1;       % Nombre de points dans l'historique au temps t
    alpha    = 0.0;     % Vecteur d'historique des maxima
    beta     = 0.0;     % Vecteur d'historique des minima
end

% Nombre de points dans l'historique au pas precedent
np = n;

% Mise-a-jour de la frontiere de Preisach
% Champ croissant
if H > Hp           
    for k = np:-1:0
        if k == 0           % Nouveau champ maximal
            n = 1;
            alpha(n) =  abs(H);
            beta(n)  = -abs(H);
            
        elseif H < alpha(k) % Wipe-out ou nouveau point
            n = k+1;
            alpha(n) = H;
            beta(n)  = beta(k);
            break
        end
    end
% Champ decroissant
elseif H < Hp       
    for k = np:-1:0
        if k == 0           % Nouveau champ maximal
            n = 1;
            alpha(n) =  abs(H);
            beta(n)  = -abs(H);
        elseif H > beta(k)  % wipe-out ou nouveau point
            n = k+1;
            alpha(n) = alpha(k);
            beta(n)  = H;
            break
        end
    end
end

% Mise-a-jour de la valeur precedente de champ
Hp = H;

% Calcul de la somme partielle si necessaire
if n ~= np
    Mpartiel = 0.0;
    for k = 2:(n-1)
        if beta(k-1) ~= beta(k)
            Mpartiel = Mpartiel + ...
                2.0*(E(alpha(k),beta(k-1))-E(alpha(k),beta(k)));
        end
    end
end

% Calcul de la magnetisation
M = -E(alpha(1),beta(1)) + 2.0*E(H,beta(n));
if n > 1
    M = M + Mpartiel;
    if beta(n-1) ~= beta(n)
        M = M + 2.0*(E(alpha(n),beta(n-1))-E(alpha(n),beta(n)));
    end
end

end % Fin de l'implementation