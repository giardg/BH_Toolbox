clc
clear all
close all

%% Creation de la courbe H(t) (peut aussi etre importee)
freq = 10*1e3;
period = 1/freq;
t = linspace(0,6*period,1000);
ramp = ones(size(t));
n = 5; %Nombre de periodes pour le ramp-up
ramp(t<n*period) = (freq/n)*t(t<n*period);
H0 = 5e3;

Ht = H0*sin(2*pi*freq*t).*ramp;


%% B(t) avec le modele de Preisach voulu
%Choix du type de matériau
disptype = 2; %1-2

if disptype == 1
    
    %Preisach 3n coeff
    a = [0.6569, 0.3038, 0.0417];
    b = [466.3, 3712.2, 4243.8 ];
    c = [1627.1, 1651.9, -9026.3];
    
    % Calcul de B(t)
    Bt = BHcurve_3n(Ht,a,b,c);
    
    
elseif disptype == 2
    
    %Preisach 4 param
    Br = 1.4472;
    Bsat = 1.8156;
    s = 0.9130;
    Hc = 1250;
    
    % Calcul de B(t)
    Bt = BHcurve_4P(Ht,Hc,Br,Bsat,s);

    

end

%% Figures

[haxes,hline1,hline2] = plotyy(t,Ht,t,Bt,'plot');
xlabel('t (s)')
set(hline1,'Color','b')
set(hline2,'Color','r')
axes(haxes(1))
ylabel('H (A/m)','Color','b')
axes(haxes(2))
ylabel('B (T)','Color','r')
set(haxes(1), 'YColor', 'b');
set(haxes(2), 'YColor', 'r');

%B(H)
figure
plot(Ht, Bt,'b')
hold off
grid on
xlabel('H (A/m)')
ylabel('B (T)')
