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

%% Parametres a modifier par l'utilisateur

% Path jusqu'au programme fortran et repertoire I/O
path = 'C:/Users/grego/Documents/Université/Maitrise/Projet/mu_1D/FEM_transitoire_fortran/SlabProblem/';

% Discretisation du domaine de calcul :
% Taille du domaine de calcul (m)
L           = 0.25e-03;
% Nombre d'elements
ne          = 100;
% Duree d'un cycle (s)
T           = 1/10.0e+03;
% Nombre de pas de temps par periode
ntpc        = 2000;

% Proprietes du materiau :
% Resistivite (Ohms/m)
rho         = 4.60e-7;
% Type de courbe BH
mattype     = 3;
% 1 : lineaire permeabilite relative constante
%     B = mu0*mur*H
if mattype == 1
    mur     = 1.0e+03;
% 2 : non lineaire modele arctangente
%     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
elseif mattype == 2
    murmax  = 6.0e+02;
    bsat    = 1.8e+00;
% 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
%     B+ = Sum ai*atan((H+ci)/bi)
%     B- = Sum ai*atan((H-ci)/bi)
elseif mattype == 3
    natan   = 2;
    ai = [0.970970, 0.117365];
    bi = [623.481589, 2523.624533 ];
    ci = [2114.306176, -6270.991790];
else
    error('Le type de materiau choisi n''est pas implemente.')
end

% Champ magnetique applique (condition de Dirichelet) :
% Amplitude 1 (A/m)
hamp1       = 50.0e+03;
% Frequence 1 (Hz)
freq1       = 1/T;
% Amplitude 2 (A/m)
hamp2       = 50.0e+03;
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
system(char(strcat(path,{'SlabProblem '},path)));

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

%% Affichage des resultats

% Temps d'affichage d'un pas de temps
frame = 10/ntpc;

% Champ magnetique
figure(1)
pH = plot(x,H(1:nx,1));
title('Champ magnetique')
xlabel('x [m]')
ylabel('H [A/m]')
axis([0,max(x),min(H(1,:)),max(H(1,:))])
for n = 2:nt
    set(pH,'Ydata',H(1:nx,n));
    drawnow;
    pause(frame)
end

% Induction
figure(2)
pB = plot(x,B(1:nx,1));
title('Induction')
xlabel('x [m]')
ylabel('B [T]')
axis([0,max(x),min(B(1,:)),max(B(1,:))])
for n = 2:nt
    set(pB,'Ydata',B(1:nx,n));
    drawnow;
    pause(frame)
end

% Densite de courant
figure(3)
pJ = plot(x,J(1:nx,1));
title('Densite de courant')
xlabel('x [m]')
ylabel('J [A/m^2]')
axis([0,max(x),min(J(1,:)),max(J(1,:))])
for n = 2:nt
    set(pJ,'Ydata',J(1:nx,n));
    drawnow;
    pause(frame)
end

% Pertes
Pj_x = mean(Pj.').';
Ph_x = mean(Ph.').';
figure(4)
plot(x,Pj_x,x,Ph_x)
title('Pertes moyennes')
legend('Pertes Joule','Pertes par hysteresis')
xlabel('x [m]')
ylabel('P [W/m^3]')
axis([0,max(x),0,1.1*max([Pj_x;Ph_x])])

% Homogeneisation
d = 2*L;
Hs = H_0_t;
Ha = trapz(H)/(2*ne);
Ba = trapz(B)/(2*ne);
dH2a = trapz((H-ones(size(x))*Ha).^2)/(2*ne);
dH3a = trapz((H-ones(size(x))*Ha).^3)/(2*ne);
figure(5)
plot(Hs,Ba)
title('Courbe Ba-Hs')
xlabel('Hs [A/m]')
ylabel('Ba [T]')








