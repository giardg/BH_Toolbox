clc;
clearvars
close all

%% Parametres
mu0 = 4*pi*1E-7;
H0_list = [10,25,50,100]*1e3;
longueur = 1.5e-3;
correction = true;
mattype = 3;

if mattype == 1
    typename = 'Linear';
elseif mattype == 2
    typename = 'Foucault';
elseif mattype == 3
    typename = 'Hysteresis';
elseif mattype == 4
    typename = 'Limites';
else
    error('Le type de materiau choisi n''est pas implemente.')
end

%% Affichage des resultats
figure(1)
figure(2)
figure(3)

for H0 = H0_list
    filename = strcat('H0',{' '},string(H0),'_mu',typename,'_',string(longueur*1e3),'mm');
    if correction filename = strcat(filename,'_corrected'); end
    filename = strcat(filename,'.txt');
    
    if isfile(filename)
        
        fid = fopen(filename,'r');
        data = textscan(fid, '%f%f%f', 'CollectOutput',1);
        data = cell2mat(data);
        H = real(data(:,1));
        mu_real = (data(:,2));
        mu_im = (data(:,3));
    else
     error('Le fichier %s n''existe pas',filename)
    end
    
    figure(1)
    semilogx(H*1e-3,mu_real,'Linewidth',2)
    %plot(H*1e-3,mu_real,'Linewidth',2)
    hold on
    
    figure(2)
    semilogx(H*1e-3,mu_im,'Linewidth',2)
    %plot(H*1e-3,mu_im,'Linewidth',2)
    hold on
    
    figure(3)
    plot(H*1e-3,(mu_real-1).*H*4*pi*1e-7,'Linewidth',2)
    hold on
    
end

figure(1)
lgd = legend('H_0 = '+string(H0_list*1e-3)+{' '}+'kA/m');
lgd.FontSize = 12;
lgd.Location = 'southwest';
xlabel('H (kA/m)','Fontsize',16)
ylabel('Re(\mu_r)','Fontsize',16)
ylim([0 inf])
xlim([1e-2 inf])


figure(2)
lgd = legend('H_0 = '+string(H0_list*1e-3)+{' '}+'kA/m');
lgd.FontSize = 8;
lgd.Location = 'southwest';
xlabel('H (kA/m)','Fontsize',16)
ylabel('Im(\mu_r)','Fontsize',16)
ylim([-inf 0])
xlim([1e-2 inf])

figure(3) 
lgd = legend('H_0 = '+string(H0_list*1e-3)+{' '}+'kA/m');
lgd.FontSize = 8;
lgd.Location = 'southwest';
xlabel('H (kA/m)', 'Fontsize',16)
ylabel('M(T)', 'Fontsize',16)
