function [] = PowerEquivalentModel(freq, L, rho, H0_list, mattype, currentFolder, output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PowerEquivalentModel.m
%
% Script principal pour resoudre le PEM en 1D et determiner mu a partir de
% donnees de puissance
%
% Etapes: 1) Resolution de P = (rho/2)*[H*H'' + (H')^2] dans fem_1D.m
%         2) Resolution de phi' = -(1/H)*sqrt((2/rho)*P_Joule - (H')^2)
%            dans resoud_phi.m
%         3) Calcul de mu' et mu''
%
% Juillet 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parametres numeriques
global param_num

param_num.nel           = 200; %nombre d'elements
param_num.nbIter        = 50; %nombre d'iterations maximal de la methode (Newton)
param_num.critere       = 1e-4; %critere d'arret de la methode (Newton) |delta|/|Hk|
param_num.deg           = 2; %degre des fonctions de base (Lagrange degre 1, 2 et 3 implementes)
param_num.nquad         = 4; %nombre de quadrature de Gauss
param_num.phi0          = 0; %phase en x = 0
param_num.sparsity      = false;
param_num.amortissement = .1; %facteur d'amortissement (convergence plus lente, mais plus stable: default = 1)
param_num.correction    = true; %Appliquer la correction sur mu (lissage et correction a faible champ)
CF                      = [1,2]; %Conditions frontieres pour les deux bords (1 = Dirichlet, 2 = Neumann)


%% Parametres physiques
mu0 = 4*pi*1E-7;
global param_phy
param_phy.rho           = rho;
param_phy.mattype       = mattype; 
param_phy.longueur      = L; %longueur du domaine
param_phy.omega         = 2*pi*freq; %frequence

%Boucle sur les H0
for H0 = H0_list
    param_phy.H0 = H0;
    
    %Boucle sur les temperatures (mesures magnetiques)
    for i = 1:length(param_phy.Temp_list)
        
        close all
        param_phy.Temp = param_phy.Temp_list(i);
        
        % 1 : lineaire permeabilite relative constante
        %     B = mu0*mur*H
        if param_phy.mattype == 1
            typename = 'Linear';
            param_phy.mur     = param_phy.mur_list(i);
            
            % 2 : non lineaire modele arctangente
            %     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
        elseif param_phy.mattype == 2
            typename = 'Foucault';
            param_phy.murmax  = param_phy.murmax_list(i);
            param_phy.Bsat = param_phy.Bsat_list(i);
            
            % 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
            %     B+ = Sum ai*atan((H+ci)/bi)
            %     B- = Sum ai*atan((H-ci)/bi)
        elseif param_phy.mattype == 3
            typename = 'Hysteresis';
            param_phy.natan   = param_phy.natan(i);
            param_phy.ai = param_phy.ai_list(i);
            param_phy.bi = param_phy.bi_list(i);
            param_phy.ci = param_phy.ci_list(i);
            
            % 4 : Cas limite a fort champ (on veut Hsat le plus petit possible, qui permet la convergence)
            %     B = mu0*H+Bsat si H > Hsat
            %     B = mu0*H-Bsat si H < Hsat
            %     polynome de degre 3 qui assure la continuité de B et de dBdH
        elseif param_phy.mattype == 4
            typename = 'Limites';
            param_phy.Bsat = param_phy.Bsat_list(i);
            param_phy.Hsat = param_phy.Hsat_list(i);
            
            % 5 : Modele EFG avec cycle majeur elliptique
        elseif param_phy.mattype == 5
            typename = 'Ellipse';
            param_phy.mur_real     = param_phy.mur_real_list(i);
            param_phy.mur_im       = param_phy.mur_im_list(i);
            
        elseif param_phy.mattype == 6
            typename = 'Ellipse2';
            param_phy.mur_real     = param_phy.mur_real_list(i);
            param_phy.mur_im       = param_phy.mur_im_list(i);
            
        elseif param_phy.mattype == 7
            typename = 'Hysteresis2';
            param_phy.Hc = param_phy.Hc_list(i);
            
        else
            error('Le type de materiau choisi n''est pas implemente.')
        end
        
        %% Chercher P_Joule et P_Hyst
        
        %addpath('Resultats_transitoires_fortran','Resultats_transitoires_COMSOL','Fonctions_FEM','Fonctions_phi')
        addpath([currentFolder output '\Resultats_transitoire_fortran'],[currentFolder '\FEM_PEM_1D\Fonctions_FEM'],[currentFolder '\FEM_PEM_1D\Fonctions_phi'])
        
        
        % Aller chercher dans les fichiers (FEM_transitoire_fortran)
        [xvec_init, P_Joule, P_Hyst] = read_puissance(currentFolder, output);
        P_exp = P_Joule+P_Hyst;
        
        % Approximation initiale de la solution (doit respecter les CF)
        syms x
        Hk=@(x) param_phy.H0+0*x;
        dHk=@(x) 0*x;
        ddHk=@(x) 0*x;
        
        Hk = Hk(xvec_init);
        dHk = dHk(xvec_init);
        ddHk = ddHk(xvec_init);
        
        
        %% Probleme elements finis pour obtenir H (Etape 1)
        
        % Resolution du probleme elements finis
        [H, dH, ddH, xvec] = fem_1D( CF, xvec_init, Hk, dHk, P_exp );
        P = (param_phy.rho/2)*((H).*ddH+dH.^2);
        
        %% Resolution de l'equation en phi (Etape 2)
        [phi, dphi, ddphi] = resoud_phi(xvec, H, dH, xvec_init, P_Joule);
        
        %% Trouver mu' et mu'' (Etape 3)
        mu_real = (1/mu0)*(param_phy.rho/param_phy.omega)*((2./H).*dH.*dphi+ddphi);
        mu_im = (1/mu0)*(-2*spline(xvec_init,P_Hyst,xvec)/param_phy.omega)./(H.^2);
        if (param_num.correction) [mu_real,mu_im] = corriger_mu(mu_real, mu_im, H, P); end %Correction de mu à faible champ
        mu = mu0*(mu_real+1i*mu_im);
        
        %% Verification (FoucaultSolve.m)
        PQ = struct('rho',param_phy.rho,'omega',param_phy.omega,'x',xvec,'dx',xvec(2)-xvec(1),'nx',length(xvec),'Hmax', param_phy.H0);
        H2 = FoucaultSolve( PQ, transpose(mu) );
        [Ptot_approx,Pjoule_approx,Physt_approx] = PowerLosses(PQ,H2,mu);
        
        
        %% Affichage des resultats
        afficher_resultats(xvec_init, P_Joule, P_Hyst, Hk, xvec, H, ...
            dH, ddH, P_exp, P, phi, mu_real, mu_im, ...
            Pjoule_approx, Physt_approx, Ptot_approx)
        
        %% Enregistrer les resultats
        filename = strcat('H0',{' '},string(param_phy.H0),'_mu',typename,'_',string(param_phy.longueur*1e3),'mm_',num2str(param_phy.Temp),'deg');
        if param_num.correction filename = strcat(filename,'_corrected'); end
        filename = strcat(filename,'.txt');
        dlmwrite(strcat(currentFolder, output, 'Resultats_mu\', filename), [H mu_real mu_im], '\t');
    end
end
