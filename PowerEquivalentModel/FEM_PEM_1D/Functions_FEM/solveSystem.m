function U = solveSystem( M, F, S, coord, bord_ess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solveSystem
% 
% Resoud le systeme matriciel afin d'obtenir les coefficients de la
% solution delta de la methode iterative
%
% 
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_num;

nnoeuds = size(coord,1);
nM11 = nnoeuds - size(bord_ess,1);

M11 = M(1:nM11,1:nM11);
M12 = M(1:nM11,((nM11+1):nnoeuds));

if param_num.sparsity   M11 = sparse(M11);   end

UC =  zeros(size(bord_ess,1),1); %delta sera toujours nul sur un bord avec condition essentielle
UI =  M11 \( F(1:nM11) + S(1:nM11) - M12*UC  );

U = [UI;UC];


