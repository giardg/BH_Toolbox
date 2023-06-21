%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve fitting tool to determine the 3n coefficients of the scalar       %
% Preisach model from experimental hysteresis data                        %
% Author : Jaël Giguère                                                   %
% Last modified : 2023-06-20 (Gregory Giard)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File structure :                                                        %
% .txt, NO BUFFER, 2 columns                                              %
% [ H applied (A/m)    M measured (T) ]                                   %
% H must follow Hmax -> 0 -> -Hmax -> 0 -> Hmax                           %
%                  COMPLETE HYST CURVE                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit equation : sum [a_i * atan((H+c_i)/b_i)]                            %
% See Maxime Tousignant's thesis p.52                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
%

clear all; clc; close all;

currentFolder = pwd;
file = 'NewFactorHystPastilleCSalu.txt'; %choose file containing experimental data
n = 4; %Number of atan

%% Retrieve data (descending hysteresis curve first in the file)
id = fopen(file,'r');
fgetl(id);
dataFile = textscan(id,'%f %f');
fclose(id);
ChampExterne = dataFile{1}; %A/m
Aimantation = dataFile{2}; %T
%% Initial data processing
%

%Select descending part
choppedAimant = Aimantation(1:floor(length(Aimantation)/2));
choppedHext = ChampExterne(1:floor(length(ChampExterne)/2));

choppedAimant2 = Aimantation(floor(length(Aimantation)/2)+1:end);
choppedHext2 = ChampExterne(floor(length(ChampExterne)/2)+1:end);

%reponse = input('Est-ce la courbe descendante? \n 1 = oui, 2 = non\n Réponse = ');
reponse = 1; %Let's just assume the data is the descending curve


%% Interpolation and scaling
%Interpolations of the data
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

%Scale x and y data
factorHext = max(maxHext,abs(minHext));
factorAimant = max(maxAimant,abs(minAimant));
interpScaledHext = interpHext./factorHext;
interpScaledAimant = interpAimant./factorAimant;

%% Fit with n = 1
%
switch n
    case 1
        [fitresultn1, gofn1, coefsFitn1, fitAimantationn1] = FitExpHystn1(interpScaledHext,interpScaledAimant);
        
        finalHext = interpScaledHext*factorHext;
        finalAimant = fitAimantationn1*factorAimant;
        
        %Figure globale courbe hystérésis fittée
        figure
        hold on;
        plot(ChampExterne,Aimantation,'LineWidth',2);
        plot([finalHext -finalHext],[finalAimant -finalAimant],'LineWidth',2);
        xlabel('Champ externe (A/m)');
        ylabel('Aimantation (T)');
        legend('Experimental hysteresis','Fit','Location','northwest');
        title(['Fit n = 1 with Rsquare = ',num2str(gofn1.rsquare),' and RMSE = ',num2str(gofn1.rmse)]);
   
        a = [coefsFitn1(1)*factorAimant ];
        b = [coefsFitn1(3)*factorHext ];
        c = [coefsFitn1(2)*factorHext ];
        
        fprintf('a = [%f]\n', a);
        fprintf('b = [%f]\n', b);
        fprintf('c = [%f]\n', c);
        
    case 2
        %% Fit données avec n = 2
        %
        [fitresultn2, gofn2, coefsFitn2, fitAimantationn2] = FitExpHystn2(interpScaledHext,interpScaledAimant);
        
        finalHext2 = interpScaledHext*factorHext;
        finalAimant2 = fitAimantationn2*factorAimant;
        
        figure
        hold on;
        plot(ChampExterne,Aimantation,'LineWidth',2);
        plot([finalHext2 -finalHext2],[finalAimant2 -finalAimant2],'LineWidth',2);
        xlabel('Champ externe (A/m)');
        ylabel('Aimantation (T)');
        legend('Experimental hysteresis','Fit','Location','northwest');
        title(['Fit n = 2 with Rsquare = ',num2str(gofn2.rsquare),' and RMSE = ',num2str(gofn2.rmse)]);
        
        a = [coefsFitn2(1)*factorAimant, coefsFitn2(4)*factorAimant ];
        b = [coefsFitn2(3)*factorHext, coefsFitn2(6)*factorHext];
        c = [coefsFitn2(2)*factorHext, coefsFitn2(5)*factorHext ];
        
        fprintf('a = [%f %f]\n', a);
        fprintf('b = [%f %f]\n', b);
        fprintf('c = [%f %f]\n', c);
        
    case 3
        
        %% Fit données avec n = 3
        %
        [fitresultn3, gofn3, coefsFitn3, fitAimantationn3] = FitExpHystn3(interpScaledHext,interpScaledAimant);
        
        finalHext3 = interpScaledHext*factorHext;
        finalAimant3 = fitAimantationn3*factorAimant;
        
        figure
        hold on;
        plot(ChampExterne,Aimantation,'LineWidth',2);
        plot([finalHext3 -finalHext3],[finalAimant3 -finalAimant3],'LineWidth',2);
        xlabel('Champ externe (A/m)');
        ylabel('Aimantation (T)');
        legend('Experimental hysteresis','Fit','Location','northwest');
        title(['Fit n = 3 with Rsquare = ',num2str(gofn3.rsquare),' and RMSE = ',num2str(gofn3.rmse)]);
        
        a = [coefsFitn3(1)*factorAimant, coefsFitn3(4)*factorAimant, coefsFitn3(7)*factorAimant ];
        b = [coefsFitn3(3)*factorHext, coefsFitn3(6)*factorHext, coefsFitn3(9)*factorHext];
        c = [coefsFitn3(2)*factorHext, coefsFitn3(5)*factorHext, coefsFitn3(8)*factorHext];
        
        fprintf('a = [%f %f %f]\n', a);
        fprintf('b = [%f %f %f]\n', b);
        fprintf('c = [%f %f %f]\n', c);
        
    case 4
        %% Fit données avec n = 4
        %
        [fitresultn4, gofn4, coefsFitn4, fitAimantationn4] = FitExpHystn4(interpScaledHext,interpScaledAimant);
        
        finalHext4 = interpScaledHext*factorHext;
        finalAimant4 = fitAimantationn4*factorAimant;
        
        figure
        hold on;
        plot(ChampExterne,Aimantation,'LineWidth',2);
        plot([finalHext4 -finalHext4],[finalAimant4 -finalAimant4],'LineWidth',2);
        xlabel('H (A/m)');
        ylabel('B (T)');
        legend('Experimental hysteresis','Fit','Location','northwest');
        title(['Fit n = 4 with Rsquare = ',num2str(gofn4.rsquare),' and RMSE = ',num2str(gofn4.rmse)]);
        
        a = [coefsFitn4(1)*factorAimant, coefsFitn4(4)*factorAimant, coefsFitn4(7)*factorAimant, coefsFitn4(10)*factorAimant];
        b = [coefsFitn4(3)*factorHext, coefsFitn4(6)*factorHext, coefsFitn4(9)*factorHext, coefsFitn4(12)*factorHext];
        c = [coefsFitn4(2)*factorHext, coefsFitn4(5)*factorHext, coefsFitn4(8)*factorHext, coefsFitn4(11)*factorHext];
        
        fprintf('a = [%f %f %f %f]\n', a);
        fprintf('b = [%f %f %f %f]\n', b);
        fprintf('c = [%f %f %f %f]\n', c);
        
end
