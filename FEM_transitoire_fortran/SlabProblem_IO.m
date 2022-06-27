function [] = slabProblem_IO(freq, L, rho, H0_list, mattype, currentFolder, output)

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

global param_phy

%Boucle sur les H0
for H0 = H0_list
    
    
    %% Parametres a modifier par l'utilisateur
    
    % Path jusqu'au programme fortran et repertoire I/O
    path =strcat(currentFolder,'\FEM_transitoire_fortran\SlabProblem\');
    
    % Nombre d'elements
    ne          = 200;
    % Duree d'un cycle (s)
    T           = 1/freq;
    % Nombre de pas de temps par periode
    ntpc        = 2000;
    
    %Boucle sur les temp
    for k = 1:length(param_phy.Temp_list)
        
        
        % 1 : lineaire permeabilite relative constante
        %     B = mu0*mur*H
        if mattype == 1
            mur     = param_phy.mur_list(k);
        % 2 : non lineaire modele arctangente
        %     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
        elseif mattype == 2
            murmax  = param_phy.murmax_list(k);
            bsat    = param_phy.bsat_list(k);
        % 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
        %     B+ = Sum ai*atan((H+ci)/bi)
        %     B- = Sum ai*atan((H-ci)/bi)
        elseif mattype == 3
            natan   = param_phy.natan_list(k);
            ai = param_phy.ai_list(k);
            ci = param_phy.ci_list(k);
            bi = param_phy.bi_list(k);
            
         % 4 : Cas limite a fort champ (on veut Hsat le plus petit possible, qui permet la convergence)
         %     B = mu0*H+Bsat si H > Hsat
         %     B = mu0*H-Bsat si H < Hsat
         %     polynome de degre 3 qui assure la continuité de B et de dBdH sinon
        elseif mattype == 4
            Bsat = param_phy.Bsat_list(k);
            Hsat = param_phy.Hsat_list(k);
            
        % 5: Modele Preisach en definissant les fonctions F et G a partir d'un
        % cycle majeur en forme d'ellipse defini par la permeabilite constante
        % complexe
        elseif mattype == 5
            mu_real = param_phy.mu_real_list(k);
            mu_imag = param_phy.mu_imag_list(k);
            
        % 6: Modele Preisach avec l'hysteresis elliptique (mineur et majeur, pas )
        elseif mattype == 6
            mu_real = param_phy.mu_real_list(k);
            mu_imag = param_phy.mu_imag_list(k);
        
        % 7: Modele Preisach a 4 parametres
        elseif mattype == 7
            Br = param_phy.Br_list(k);
            Bsat = param_phy.Bsat_list(k);
            Hc = param_phy.Hc_list(k);
            s = param_phy.s_list(k);
            a = param_phy.a_list(k);
            b = param_phy.b_list(k);
            mur_max = param_phy.mur_max_list(k);
            
        else
            error('Le type de materiau choisi n''est pas implemente.')
        end
        
        % Champ magnetique applique (condition de Dirichelet) :
        % Amplitude 1 (A/m)
        %hamp1       = sqrt(1/(A-(B^2)/(4*C)));
        hamp1       = H0;
        
        % Frequence 1 (Hz)
        freq1       = 1/T;
        % Amplitude 2 (A/m)
        hamp2       = 0.0e+03;
        % Frequence 2 (Hz)
        freq2       = 190.0e+03;
        % Vecteur temps (s)
        t = (0:ntpc-1)*T/ntpc;
        % Champ magnetique applique
        
        % Sin
        H_0_t = hamp1*sin(2*pi*freq1*t) + hamp2*sin(2*pi*freq2*t);
        
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
        elseif mattype == 5
            fprintf(fio,'%1.15e\n',mu_real);
            fprintf(fio,'%1.15e\n',mu_imag);
            fprintf(fio,'%1.15e\n',hamp1);
        elseif mattype == 6
            fprintf(fio,'%1.15e\n',mu_real);
            fprintf(fio,'%1.15e\n',mu_imag);
            fprintf(fio,'%1.15e\n',hamp1);
        elseif mattype == 7
            fprintf(fio,'%1.15e\n',Br);
            fprintf(fio,'%1.15e\n',Bsat);
            fprintf(fio,'%1.15e\n',Hc);
            fprintf(fio,'%1.15e\n',s);
            fprintf(fio,'%1.15e\n',a);
            fprintf(fio,'%1.15e\n',b);
            fprintf(fio,'%1.15e\n',mur_max);
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
            nomFich = ['PertesLinear_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
        elseif  mattype == 2
            nomFich = ['PertesFoucault_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
        elseif mattype == 3
            nomFich = ['PertesHysteresis_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
        elseif mattype == 4
            nomFich = ['PertesLimites_' num2str(L*1E3) 'mm.txt'];
        elseif mattype == 5
            nomFich = ['PertesEllipse_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
        elseif mattype == 6
            nomFich = ['PertesEllipse2_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
        elseif mattype == 7
            nomFich = ['PertesHysteresis2_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
        end
        nomFich = ['H0 ' num2str(hamp1) '_' nomFich];
        dlmwrite([currentFolder output 'Resultats_transitoire_fortran\' nomFich], [x Pj_x Ph_x], '\t')
        
        delete([path 'Results.txt'])
        
        
    end
end






