function [coord, connec, bord] = mesh(xMin,xMax,nEls)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh
% 
% Cree le maillage selon les bornes, le nombre d'elements et le degre des
% fonctions de formes (implemente jusqu'au degre 3)
%
% Inputs: - xMin, xMax
%         - nEls: nombre d'elements
%
% Outputs: - coord [nng x 1]: vecteur des coordonnes des noeuds (les premiers noeuds
%            sont les noeuds geometriques)
%          - connec [nEls x nnoeuds/element]: matrice de connectivite
%          - bord [1 x 2]: noeuds associes aux bords
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_num
deg = param_num.deg;

coord = (xMin: (xMax-xMin)/nEls: xMax)';
connec = zeros(nEls,2);  %init matrice connectivite

nng = size(coord,1);                     % nombre de noeuds geometriques

bord = [1;nng];

i=1;
for j=1:nEls
    connec(j,1)=i;
    i=i+1;
    connec(j,2)=i;
end

if deg == 2
   coord = [coord;zeros(nEls,1)];
   connec = [connec,zeros(nEls,1)];
   for k=1:nEls
       coord(nng+k) = (coord(k+1)-coord(k))/2 + coord(k); 
       connec(k,3) = nng+k;
   end
    
end

if deg == 3
   coord = [coord;zeros(2*nEls,1)];
   connec = [connec,zeros(nEls,2)];
   for k=1:2*nEls
       ele = ceil(k/2);
       if rem(k,2) == 0
           coord(nng+k) = 2*(coord(k/2+1)-coord(k/2))/3 + coord(k/2);
           connec(k/2,4) = nng+k;
       else
           coord(nng+k) = (coord(ceil(k/2)+1)-coord(ceil(k/2)))/3 + coord(ceil(k/2));
           connec(ceil(k/2),3) = nng+k;
       end
   end
    
end



end