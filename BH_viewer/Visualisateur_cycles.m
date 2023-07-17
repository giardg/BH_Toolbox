clc
clear all
close all

%Choix du type de materiau
disptype = 2; %1-2

if disptype == 1
    
    %Preisach 3n coeff
    a = [0.6569, 0.3038, 0.0417];
    b = [466.3, 3712.2, 4243.8 ];
    c = [1627.1, 1651.9, -9026.3];
    
    %Affichage des courbes B-H (mineurs/majeur)
    DisplayMinorLoops_3N(a,b,c)
  
    
elseif disptype == 2
    
    %Preisach 4 param
    Br = 1.4472;
    Bsat = 1.8156;
    s = 0.9130;
    Hc = 1250;
    
    %Affichage des courbes B-H (mineurs/majeur)
    DisplayMinorLoops_4P(Hc,Br,Bsat,s);
    
end
