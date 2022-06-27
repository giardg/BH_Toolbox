function [Hk, dHk, ddHk, xvec] = fem1D( CF, xvec_init, Hk, dHk, P_init )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fem_1D
%
% Applique la methode des elements finis pour resoudre l'equation de
% puissance: P = (rho/2)*[H*H'' + (H')^2].
%
% *Semble bien fonctionner, mais je dois encore lire pour mieux
% comprendre comment gerer les conditions frontieres avec l'algorithme
% iteratif*
%
% Inputs: - nombre d'elements (nel), nombre d'iteration de la methode
%           iterative (nbIter), critere d'arret (critere)
%         - CF [2 x 1]: type de condition frontiere pour chacun des bords
%           (1: condition essentielle, 2: condition naturelle)
%         - xvec_init: vecteur de coordonnees pour lesquelles les fonctions
%           Hk, dHk et P_init sont evaluees
%         - Hk et dHk: guess initiaux du champ et de sa derivee
%         - P_init: donnees de puissance aux points de xvec_init
%
% Outputs: - vecteurs mis a jour (Hk, dHk, ddHk) aux points de xvec
%            (determine par le nombre d'elements)
%
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param_phy
global param_num

nel = param_num.nel;
nbIter = param_num.nbIter;
critere = param_num.critere;

%% Maillage
xMin = 0;
xMax = param_phy.longueur;
[coord, connec, bord] = maillage(xMin,xMax,nel);

%% Numerotation
[numer, bord_ess] = numerotation( coord, bord, CF );

%% Boucle d'iterations pour la methode de Newton
xvec = xvec_init;
delta  =zeros(size(xvec));
iter = 1;
critere_vec = 1000*ones(size(Hk));

h=figure(1);
x0=400;
y0=150;
width=800;
height=625;
set(gcf,'position',[x0,y0,width,height])
pH1 =plot(xvec,Hk,'Linewidth',2);
hold on
pH2 =plot(xvec,delta,'Linewidth',2);
lgd=legend('$H^i$','$\delta_H^i$');
lgd.FontSize=20;
lgd.Interpreter='latex';
ax = gca;
ax.FontSize = 22; 
xlabel('$x$ (m)','Interpreter','latex','Fontsize',28)
ylabel('$H$ (A/m)','Interpreter','latex','Fontsize',28)



while iter < nbIter && max(critere_vec) > critere
    
    P = spline(xvec_init,P_init,xvec);
    
    %% Assemblage des matrices globales (M,F,S)
    [M,F,S] = assemblage(coord, connec, bord, numer, bord_ess, xvec, Hk, dHk, P);
    
    %% Resolution du systeme matriciel
    [U] = resoudSysteme( M, F, S, coord, bord_ess );
    
    %% Mise a jour de H et de sa derivee avec delta (H_{k+1} = H_k+delta)
    Hk_old = Hk;
    xvec_old = xvec;
    [Hk, dHk, ddHk, xvec, delta, flag] = updateChamp(coord, connec, numer, xvec, Hk, U);
    
    %% Update de la figure (H et delta vs x)
    set(pH1,'Ydata',Hk_old);
    set(pH1,'Xdata',xvec_old);
    set(pH2,'Ydata',delta);
    set(pH2,'Xdata',xvec);
    drawnow;

    %% Verification de l'atteinte de convergence
    if flag
        break %arret manuel de l'amortissement, on considere la convergence
    end
    
    critere_vec = abs(delta)./abs(Hk); %verification avec le critere d'arret
    iter = iter+1;
end

if iter == nbIter
    fprintf('Convergence pas atteinte, arret a %d iterations \n',nbIter);
else
    fprintf('Convergence apres %d iterations \n',iter);
end