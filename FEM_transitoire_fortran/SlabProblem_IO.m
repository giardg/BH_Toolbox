% -------------------------------------------------------------------------
% Programme matlab qui lance le calcul elements finis en fortran pour la
% resolution du slab problem.
% Ce programme gere les Inputs/Outputs via des fichiers textes
% dans le repertoire ou se trouve l'executable du slab problem.
%
% Input files
% Parameters.txt    : Contient les parametres physiques
% Dirichelet.txt    : Contient la condition de Dirichelet H(x=0,t)
%
% Output files
% Results.txt       : Contient H(x,t), B(x,t), J(x,t), Pj(x,t) et Ph(x,t)
%
% -------------------------------------------------------------------------
clear
clc
close all

H0_vec = [50]*1e3; %A/m Liste de H0

%% Parametres a modifier par l'utilisateur

% Path jusqu'au programme fortran et repertoire I/O
%path = 'C:/Dropbox/Polytechnique/Cycle Sup/Projet de maitrise/Simulations/Calcul FEM transitoire fortran_new/SlabProblem/';
currentFolder = pwd;
path =strcat(currentFolder,'\SlabProblem\');

% Discretisation du domaine de calcul :
% Taille du domaine de calcul (m)
L           = 1.5e-03;
% Nombre d'elements
ne          = 200;
% Duree d'un cycle (s)
T           = 1/10.0e+03;
% Nombre de pas de temps par periode
ntpc        = 2000;

% Proprietes du materiau :
% Resistivite (Ohms/m)
rho         = 2.5e-7;
% Type de courbe BH
mattype     = 3;
% 1 : lineaire permeabilite relative constante
%     B = mu0*mur*H
if mattype == 1
    mur     = 0.10e+03;
    % 2 : non lineaire modele arctangente
    %     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
elseif mattype == 2
    murmax  = 6.0e+02;
    bsat    = 1.8e+00;
    % 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
    %     B+ = Sum ai*atan((H+ci)/bi)
    %     B- = Sum ai*atan((H-ci)/bi)
elseif mattype == 3
    natan   = 1;
    ai = 3.6/pi;
    bi = 1.7e3/(tan(1.4/ai));
    ci = 1.7e3;
    % 4 : Cas limite a fort champ (on veut Hsat le plus petit possible, qui permet la convergence)
    %     B = mu0*H+Bsat si H > Hsat
    %     B = mu0*H-Bsat si H < Hsat
    %     polynome de degre 3 qui assure la continuité de B et de dBdH sinon
elseif mattype == 4
    Bsat = 1;
    Hsat = 20;
else
    error('Le type de materiau choisi n''est pas implemente.')
end

% Champ magnetique applique (condition de Dirichelet) :
% Amplitude 1 (A/m)
%hamp1       = 100.0e+03;

for hamp1 = H0_vec
    
    % Frequence 1 (Hz)
    freq1       = 1/T;
    % Vecteur temps (s)
    t = (0:ntpc-1)*T/ntpc;
    % Champ magnetique applique
    
    % Sin
    H_0_t = hamp1*sin(2*pi*freq1*t);
    
    % Square
    % H_0_t = -hamp1*ones(size(t));
    % H_0_t(t>=0.5*T) = hamp1;
    
    %% Ecriture des inputs :
    % Ouverture du fichier Parameters.txt
    fio = fopen(strcat(path,'Parameters.txt'),'wt');
    
    % Ecriture
    fprintf(fio,'Duree d''un cycle (s)\n');
    fprintf(fio,'%1.15e\n',T);
    fprintf(fio,'Nombre de pas de temps par cycle\n');
    fprintf(fio,'%u\n',ntpc);
    fprintf(fio,'Taille du domaine de calcul (m)\n');
    fprintf(fio,'%1.15e\n',L);
    fprintf(fio,'Nombre d''elements\n');
    fprintf(fio,'%u\n',ne);
    fprintf(fio,'Resistivite (Ohms/m)\n');
    fprintf(fio,'%1.15e\n',rho);
    fprintf(fio,'Type de courbe BH\n');
    fprintf(fio,'%u\n',mattype);
    fprintf(fio,'Parametres de la courbe BH\n');
    if mattype == 1
        fprintf(fio,'%1.15e\n',mur);
    elseif mattype == 2
        fprintf(fio,'%1.15e\n',murmax);
        fprintf(fio,'%1.15e\n',bsat);
    elseif mattype == 3
        fprintf(fio,'%1.15e\n',natan);
        for i = 1:natan
            fprintf(fio,'%1.15e\n',ai(i));
        end
        for i = 1:natan
            fprintf(fio,'%1.15e\n',bi(i));
        end
        for i = 1:natan
            fprintf(fio,'%1.15e\n',ci(i));
        end
    elseif mattype == 4
        fprintf(fio,'%1.15e\n',Bsat);
        fprintf(fio,'%1.15e\n',Hsat);
    else
        error('Le type de materiau choisi n''est pas implemente.')
    end
    
    % Fermeture du fichier
    fclose(fio);
    
    % Ouverture du fichier Dirichelet.txt
    fio = fopen(strcat(path,'Dirichelet.txt'),'wt');
    
    % Ecriture
    for n = 1:ntpc
        fprintf(fio,'%1.15e\n',H_0_t(n));
    end
    
    % Fermeture du fichier
    fclose(fio);
    
    %% Lancement de la simulation
    % system(char(strcat(path,{'SlabProblem '},path)));
    pathSys = [path '/SlabProblem/'];
    cd(pathSys)
    system('SlabProblem.exe '+string(path))
    % system('SlabProblem_tour.exe')
    
    %% Lecture des resultats
    % ouverture du fichier
    fio = fopen(strcat(path,'Results.txt'),'rt');
    
    % Lecture de l'en-tete
    fgetl(fio); fgetl(fio); fgetl(fio);
    
    % Recuperation des dimensions
    nt = fscanf(fio,'Nombre de pas de temps : %i\n',1);
    nx = fscanf(fio,'Nombre de points en x  : %i\n',1);
    
    % Initialisation des vecteurs de donnes
    t  = zeros(1,nt);
    x  = zeros(nx,1);
    H  = zeros(nx,nt);
    B  = zeros(nx,nt);
    J  = zeros(nx,nt);
    Pj = zeros(nx,nt);
    Ph = zeros(nx,nt);
    
    % Lecture des donnees
    for n = 1:nt
        
        % Lecture de l'en-tete du pas de temps
        fgetl(fio);
        t(n) = fscanf(fio,'Time = %e\n',1);
        fgetl(fio);
        
        % Lecture de la solution au pas de temps n
        data = fscanf(fio,'%e',[6,nx]);
        fgetl(fio); fgetl(fio);
        
        % Copie des donnes dans les structures
        x(1:nx)    = data(1,1:nx)';
        H(1:nx,n)  = data(2,1:nx)';
        B(1:nx,n)  = data(3,1:nx)';
        J(1:nx,n)  = data(4,1:nx)';
        Pj(1:nx,n) = data(5,1:nx)';
        Ph(1:nx,n) = data(6,1:nx)';
        
    end
    
    % fermeture du fichier
    fclose(fio);
    
    Pj_x = mean(Pj.').';
    Ph_x = mean(Ph.').';
    
    if mattype == 1
        nomFich = ['PertesLinear_' num2str(L*1E3) 'mm.txt'];
    elseif  mattype == 2
        nomFich = ['PertesFoucault_' num2str(L*1E3) 'mm.txt'];
    elseif mattype == 3
        nomFich = ['PertesHysteresis_' num2str(L*1E3) 'mm.txt'];
    elseif mattype == 4
        nomFich = ['PertesLimites_' num2str(L*1E3) 'mm.txt'];
    end
    nomFich = ['H0 ' num2str(hamp1) '_' nomFich];
    dlmwrite([currentFolder '\..\Resultats_transitoire_fortran_aisi4340\' nomFich],[x Pj_x Ph_x], '\t')
    
    delete([path 'Results.txt'])
    
end







