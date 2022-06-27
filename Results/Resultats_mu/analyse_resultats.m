clc;
clearvars
close all
fclose all;

x0=400;
y0=400;
width=700;
height=550;

%% Parametres
mu0 = 4*pi*1E-7;
H0 = 50*1e3; %liste des fichiers disponible
freq_list = [10,20,50,100,150,200]; %kHz
Temp = 25;
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
elseif param_phy.mattype == 5
    typename = 'Ellipse';
elseif param_phy.mattype == 6
    typename = 'Ellipse2';
elseif param_phy.mattype == 7
    typename = 'Hysteresis2';
    
else
    error('Le type de materiau choisi n''est pas implemente.')
end


figure(1)
set(gcf,'position',[x0,y0,width,height])
figure(2)
set(gcf,'position',[x0,y0,width,height])
figure(3)

cpt = 1;


%% Ouverture de tous les fichiers
for freq = freq_list
    filename = strcat('H0',{' '},string(H0),'_mu',typename,'_',string(Temp),'deg_',string(freq),'kHz');
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
    
    %     if rem(H0*1e-3,5) ~= 0
    %         cpt = cpt+1;
    %         continue
    %     end
    
    figure(1)
    semilogx(H*1e-3,mu_real,'Linewidth',2)
    %plot(H*1e-3,mu_real,'Linewidth',2)
    hold on
    
    figure(2)
    semilogx(H*1e-3,-mu_im,'Linewidth',2)
    %plot(H*1e-3,mu_im,'Linewidth',2)
    hold on
    
    figure(3)
    plot(H*1e-3,(mu_real-1).*H*4*pi*1e-7,'Linewidth',2)
    hold on
    cpt = cpt+1;
end

figure(1)

lgd = legend('10 kHz', '20 kHz', '50 kHz', '100 kHz', '150 kHz', '200 kHz', 'Location','Best');
lgd.FontSize = 14;
ax = gca;
ax.FontSize = 14; 
xlabel('$H$ (kA/m)','Interpreter','latex','Fontsize',20)
ylabel('$\Re(\mu_r)$','Interpreter','latex','Fontsize',20)

figure(2)

lgd = legend('10 kHz', '20 kHz', '50 kHz', '100 kHz', '150 kHz', '200 kHz', 'Location','Best');
lgd.FontSize = 14;
ax = gca;
ax.FontSize = 14; 
xlabel('$H$ (kA/m)','Interpreter','latex','Fontsize',20)
ylabel('$\Im(\mu_r)$','Interpreter','latex','Fontsize',20)
