function [numer, bord_ess] = numbering( coord, bord, CF )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerotation
% 
% Cree la matrice numer qui numerote les degres de libertes du systeme (en
% terminant par les bords avec conditions essentielles)
%
% Inputs: - coord, bord: vecteur de coordonnes et des bords cree dans
%           maillage
%         - CF [2 x 1]: type de condition frontiere pour chacun des bords
%           (1: condition essentielle, 2: condition naturelle)
%
% Outputs: - numer [nnoeuds x 1]: vecteur de numerotation des degres de
%            libertes (on considere seulement un degre de liberte par noeud
%            a cause du probleme a resoudre)
%          - bord_ess: liste des bords avec conditions essentielles
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nnoeuds      = size(coord,1);
nnoeuds_bord = size(bord,1);
numer        = zeros(nnoeuds,1); %1 degre de liberte par noeud

bord_Neumann = CF==2;
bord_ess = bord(~bord_Neumann);

iddl = 1;
for k=1:nnoeuds
    if (~ismember(k,bord_ess))  % si le noeud k n'est pas un bord essentiel
        numer(k) = iddl;     
        iddl=iddl+1;
    end
end

for k=1:nnoeuds_bord
   if (~bord_Neumann(k))
       numer(bord(k)) = iddl;
       iddl = iddl + 1;
   end
end



end