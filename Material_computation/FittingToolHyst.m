%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outil de curve fitting pour déterminer les 3n coefficients d'une courbe %
% d'hysétéris expérimentale.                                              %
% Auteur : Jaël Giguère                                                   %
% Dernière modification : 2021-11-04                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure du fichier à fournir :                                        %
% .txt, PAS d'entête, 2 colonnes                                          %
% [ H appliqué (A/m)    M mesuré (T) ]                                    %
% H doit être Hmax -> 0 -> -Hmax -> 0 -> Hmax                             %
%         ou -Hmax -> 0 -> Hmax -> 0 -> -Hmax                             %
%                  COURBE HYST COMPLÈTE                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forme du fit : sum [a_i * atan((H+c_i)/b_i)]                            %   
% Voir thèse de Maxime Tousignant p.52                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
%

clear all; clc; close all;

cd DonneesExpHyst %change folder pour choisir fichier
file = uigetfile('*.txt','Choose the file to fit'); %choisit le fichier
mkdir(file(1:end-4)); %crée un folder pour mettre les résultats du fit
id = fopen(file,'r');
fgetl(id);
dataFile = textscan(id,'%f %f');
fclose(id);
cd .. %revient au directory contenant les fcts Matlab
ChampExterne = dataFile{1}; %en A/m
Aimantation = dataFile{2}; %en Tesla
%% Traitement initial données
%

%Coupe les données pour sélectionner la courbe descendante
choppedAimant = Aimantation(1:floor(length(Aimantation)/2));
choppedHext = ChampExterne(1:floor(length(ChampExterne)/2));

choppedAimant2 = Aimantation(floor(length(Aimantation)/2)+1:end);
choppedHext2 = ChampExterne(floor(length(ChampExterne)/2)+1:end);

%Demande au user si c'est la courbe descendante
figure(1)
hold on;
plot(ChampExterne,Aimantation,'LineWidth',8);
plot(choppedHext,choppedAimant,'LineWidth',1.5);
title('Est-ce la courbe descendante?');
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hystérésis expérimentale','Courbe descendante?');

reponse = input('Est-ce la courbe descendante? \n 1 = oui, 2 = non\n Réponse = ');

close(figure(1));

%% Interpolation et scaling
%Fais une interpolation des données pour avoir plus de points pour faire le
%fit sur la courbe descendate
nbrePtsInterp = 1000;

if reponse == 1
    minHext = min(choppedHext);
    maxHext = max(choppedHext);
    maxAimant = max(choppedAimant);
    minAimant = min(choppedAimant);
    interpHext = linspace(minHext,maxHext,nbrePtsInterp);
    interpAimant = interp1(choppedHext,choppedAimant,interpHext);
    
elseif reponse == 2
    minHext = min(choppedHext2);
    maxHext = max(choppedHext2);
    maxAimant = max(choppedAimant2);
    minAimant = min(choppedAimant2);
    interpHext = linspace(minHext,maxHext,nbrePtsInterp);
    interpAimant = interp1(choppedHext2,choppedAimant2,interpHext);
end

figure(1001)
disp('Voici la courbe descendante interpolée pour le fit');
hold on;
plot(ChampExterne,Aimantation,'LineWidth',2);
plot(interpHext,interpAimant,'o');
title('Courbe interpolée pour le fit');
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hytérésis expérimentale','Interpolation');

%Scale données x et y pour que valeur maximale soit +-1
factorHext = max(maxHext,abs(minHext));
factorAimant = max(maxAimant,abs(minAimant));
interpScaledHext = interpHext./factorHext;
interpScaledAimant = interpAimant./factorAimant;

%% Fit données avec n = 1
%
[fitresultn1, gofn1, coefsFitn1, fitAimantationn1] = FitExpHystn1(interpScaledHext,interpScaledAimant);

finalHext = interpScaledHext*factorHext;
finalAimant = fitAimantationn1*factorAimant;

%Figure globale courbe hystérésis fittée
figure(1)
hold on;
plot(ChampExterne,Aimantation,'LineWidth',2);
plot([finalHext -finalHext],[finalAimant -finalAimant],'LineWidth',2);
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hystérésis expérimentale','Courbe hystérésis fittée','Location','northwest');
title(['Fit n = 1 avec Rsquare = ',num2str(gofn1.rsquare),' et RMSE = ',num2str(gofn1.rmse)]);

%Sauvegarde les données si demandé
displayn1 = input('Infos sur fit n = 1?\n 1 = oui, 2 = non\nRéponse = ');
if displayn1 == 1
    cd DonneesExpHyst
    cd(file(1:end-4))
    fid= fopen(strcat('Fitn=1',file(1:end-4),'.txt'),'w');
    fprintf(fid,'Fit n = 1 avec Rsquare = %f et RMSE = %f\n\nn = 1\ta_1=%e\tb_1=%f\tc_1=%f',gofn1.rsquare,gofn1.rmse,coefsFitn1(1)*factorAimant,coefsFitn1(3)*factorHext,coefsFitn1(2)*factorHext);
    fclose(fid);
    cd ..
    cd ..
    fprintf('Équation du fit n = 1: %e*atan((x+%f)/%f))\n',coefsFitn1(1)*factorAimant,coefsFitn1(2)*factorHext,coefsFitn1(3)*factorHext);
end

%% Fit données avec n = 2
%
[fitresultn2, gofn2, coefsFitn2, fitAimantationn2] = FitExpHystn2(interpScaledHext,interpScaledAimant);

finalHext2 = interpScaledHext*factorHext;
finalAimant2 = fitAimantationn2*factorAimant;

figure(2)
hold on;
plot(ChampExterne,Aimantation,'LineWidth',2);
plot([finalHext2 -finalHext2],[finalAimant2 -finalAimant2],'LineWidth',2);
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hystérésis expérimentale','Courbe hystérésis fittée','Location','northwest');
title(['Fit n = 2 avec Rsquare = ',num2str(gofn2.rsquare),' et RMSE = ',num2str(gofn2.rmse)]);

%Sauvegarde les données si demandé
displayn2 = input('Infos sur fit n = 2?\n 1 = oui, 2 = non\nRéponse = ');
if displayn2 == 1
    cd DonneesExpHyst
    cd(file(1:end-4))
    fid= fopen(strcat('Fitn=2',file(1:end-4),'.txt'),'w');
    fprintf(fid,'Fit n = 2 avec Rsquare = %f et RMSE = %f\n\n n = 1\ta_1=%e\tb_1=%f\tc_1=%f\n n = 2\ta_2=%e\tb_2=%f\tc_2=%f',gofn2.rsquare,gofn2.rmse,coefsFitn2(1)*factorAimant,coefsFitn2(3)*factorHext,coefsFitn2(2)*factorHext,coefsFitn2(4)*factorAimant,coefsFitn2(6)*factorHext,coefsFitn2(5)*factorHext);
    fclose(fid);
    cd ..
    cd ..
    fprintf('Équation du fit n = 2 : %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f))\n',coefsFitn2(1)*factorAimant,coefsFitn2(2)*factorHext,coefsFitn2(3)*factorHext,coefsFitn2(4)*factorAimant,coefsFitn2(5)*factorHext,coefsFitn2(6)*factorHext);
end

%% Fit données avec n = 3
%
[fitresultn3, gofn3, coefsFitn3, fitAimantationn3] = FitExpHystn3(interpScaledHext,interpScaledAimant);

finalHext3 = interpScaledHext*factorHext;
finalAimant3 = fitAimantationn3*factorAimant;

figure(3)
hold on;
plot(ChampExterne,Aimantation,'LineWidth',2);
plot([finalHext3 -finalHext3],[finalAimant3 -finalAimant3],'LineWidth',2);
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hystérésis expérimentale','Courbe hystérésis fittée','Location','northwest');
title(['Fit n = 3 avec Rsquare = ',num2str(gofn3.rsquare),' et RMSE = ',num2str(gofn3.rmse)]);

%Sauvegarde les données si demandé
displayn3 = input('Infos sur fit n = 3?\n 1 = oui, 2 = non\nRéponse = ');
if displayn3 == 1
    cd DonneesExpHyst
    cd(file(1:end-4))
    fid= fopen(strcat('Fitn=3',file(1:end-4),'.txt'),'w');
    fprintf(fid,'Fit n = 3 avec Rsquare = %f et RMSE = %f\n\n n = 1\ta_1=%e\tb_1=%f\tc_1=%f\n n = 2\ta_2=%e\tb_2=%f\tc_2=%f\n n = 3\ta_3=%e\tb_3=%f\tc_3=%f',gofn3.rsquare,gofn3.rmse,coefsFitn3(1)*factorAimant,coefsFitn3(3)*factorHext,coefsFitn3(2)*factorHext,coefsFitn3(4)*factorAimant,coefsFitn3(6)*factorHext,coefsFitn3(5)*factorHext,coefsFitn3(7)*factorAimant,coefsFitn3(9)*factorHext,coefsFitn3(8)*factorHext);
    fclose(fid);
    cd ..
    cd ..
    fprintf('Équation du fit n = 3 : %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f))\n',coefsFitn3(1)*factorAimant,coefsFitn3(2)*factorHext,coefsFitn3(3)*factorHext,coefsFitn3(4)*factorAimant,coefsFitn3(5)*factorHext,coefsFitn3(6)*factorHext,coefsFitn3(7)*factorAimant,coefsFitn3(8)*factorHext,coefsFitn3(9)*factorHext);
end

%% Fit données avec n = 4
%
[fitresultn4, gofn4, coefsFitn4, fitAimantationn4] = FitExpHystn4(interpScaledHext,interpScaledAimant);

finalHext4 = interpScaledHext*factorHext;
finalAimant4 = fitAimantationn4*factorAimant;

figure(4)
hold on;
plot(ChampExterne,Aimantation,'LineWidth',2);
plot([finalHext4 -finalHext4],[finalAimant4 -finalAimant4],'LineWidth',2);
xlabel('Champ externe (A/m)');
ylabel('Aimantation (T)');
legend('Courbe hystérésis expérimentale','Courbe hystérésis fittée','Location','northwest');
title(['Fit n = 4 avec Rsquare = ',num2str(gofn4.rsquare),' et RMSE = ',num2str(gofn4.rmse)]);

%Sauvegarde les données si demandé
displayn4 = input('Infos sur fit n = 4?\n 1 = oui, 2 = non\nRéponse = ');
if displayn4 == 1
    cd DonneesExpHyst
    cd(file(1:end-4))
    fid= fopen(strcat('Fitn=4',file(1:end-4),'.txt'),'w');
    fprintf(fid,'Fit n = 4 avec Rsquare = %f et RMSE = %f\n\n n = 1\ta_1=%e\tb_1=%f\tc_1=%f\n n = 2\ta_2=%e\tb_2=%f\tc_2=%f\n n = 3\ta_3=%e\tb_3=%f\tc_3=%f\n n = 4\ta_4=%e\tb_4=%f\tc_4=%f',gofn4.rsquare,gofn4.rmse,coefsFitn4(1)*factorAimant,coefsFitn4(3)*factorHext,coefsFitn4(2)*factorHext,coefsFitn4(4)*factorAimant,coefsFitn4(6)*factorHext,coefsFitn4(5)*factorHext,coefsFitn4(7)*factorAimant,coefsFitn4(9)*factorHext,coefsFitn4(8)*factorHext,coefsFitn4(10)*factorAimant,coefsFitn4(12)*factorHext,coefsFitn4(11)*factorHext);
    fclose(fid);
    cd ..
    cd ..
    fprintf('Équation du fit n = 4 : %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f)) + %e*atan((x+%f)/%f))\n',coefsFitn4(1)*factorAimant,coefsFitn4(2)*factorHext,coefsFitn4(3)*factorHext,coefsFitn4(4)*factorAimant,coefsFitn4(5)*factorHext,coefsFitn4(6)*factorHext,coefsFitn4(7)*factorAimant,coefsFitn4(8)*factorHext,coefsFitn4(9)*factorHext,coefsFitn4(10)*factorAimant,coefsFitn4(11)*factorHext,coefsFitn4(12)*factorHext);
end
