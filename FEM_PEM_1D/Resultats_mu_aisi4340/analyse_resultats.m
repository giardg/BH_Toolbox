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
Tempvec = [24.7297,113.3906,246.0131,409.2073,461.1946,506.9757,522.4292,547.1615,560.8959,584.9624,607.1725,620.1256,632.5545,645.9273,663.0356]; %liste des fichiers disponible
%Tempvec = [547.1615,560.8959,584.9624,607.1725,620.1256,632.5545,645.9273,663.0356];
colors = {'b', 'r', 'g', 'c', 'm','y' , [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], ...
    [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840], [0.75, 0, 0.75] ,'k'};
longueur = 2.5e-3;
correction = true;
mattype = 7;
H0 = 50e3;

Temp_grid = [];
H_grid = [];
mu_real_grid = [];
mu_im_grid = [];

Temp_export = [];
H_export = [];
mu_real_export = [];
mu_im_export = [];

if mattype == 1
    typename = 'Linear';
elseif mattype == 2
    typename = 'Foucault';
elseif mattype == 3
    typename = 'Hysteresis';
elseif mattype == 4
    typename = 'Limites';
elseif mattype == 7
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
for Temp = Tempvec
    filename = strcat('H0',{' '},string(H0),'_mu',typename,'_',string(longueur*1e3),'mm_',num2str(Temp),'deg');
    if correction filename = strcat(filename,'_corrected'); end
    filename = strcat(filename,'.txt');
    
    if isfile(filename)
        
        fid = fopen(filename,'r');
        data = textscan(fid, '%f%f%f', 'CollectOutput',1);
        data = cell2mat(data);
        H = real(data(:,1));
        mu_real = (data(:,2));
        mu_im = (data(:,3));
        
        Temp_grid = [Temp_grid,ones(size(H)).*Temp];
        H_grid = [H_grid,H];
        mu_real_grid = [mu_real_grid,mu_real];
        mu_im_grid = [mu_im_grid,mu_im];
        
        Temp_export = [Temp_export;ones(size(H)).*Temp];
        
        H_export = [H_export;H];
        mu_real_export = [mu_real_export;mu_real];
        mu_im_export = [mu_im_export;mu_im];
        
        
    else
     error('Le fichier %s n''existe pas',filename)
    end
    
    figure(1)
    semilogx(H*1e-3,mu_real,'Color',colors{cpt},'Linewidth',2)
    %plot(H*1e-3,mu_real,'Linewidth',2)
    hold on
    
    figure(2)
    semilogx(H*1e-3,-mu_im,'Color',colors{cpt},'Linewidth',2)
    %plot(H*1e-3,mu_im,'Linewidth',2)
    hold on
    
    figure(3)
    plot(H*1e-3,(mu_real-1).*H*4*pi*1e-7,'Linewidth',2)
    hold on
    cpt = cpt+1;
end


%% Figure

x0=400;
y0=150;
width=1000;
height=850;

Tempvec = [24.7297,113.3906,246.0131,409.2073,461.1946,506.9757,522.4292,547.1615,560.8959,584.9624,607.1725,620.1256,632.5545,645.9273,663.0356];
Temp_str = string(Tempvec);
formatSpec = '%.0f';
for i = 1:length(Tempvec)
   Temp_str(i) = num2str(Tempvec(i),formatSpec); 
end


figure(1)
set(gcf,'position',[x0,y0,width,height])
xlim([1e-2 inf])
lgd = legend('$T$ = '+Temp_str+{' '}+'$^{o} C$','Interpreter','latex');
lgd.FontSize = 16;
lgd.Location = 'Best';
ax = gca;
ax.FontSize = 14; 
%lgd.Location = 'southwest';
xlabel('$H$ (kA/m)','Interpreter','latex','Fontsize',20)
ylabel('$\Re(\overline{\mu}_r)$','Interpreter','latex','Fontsize',20)
%ylim([0 inf])
xlim([1e-1 inf])
%saveas(gcf,'MU_REAL.png')

figure(2)
set(gcf,'position',[x0,y0,width,height])
lgd = legend('$T$ = '+Temp_str+{' '}+'$^{o} C$','Interpreter','latex');
lgd.FontSize = 16;
lgd.Location = 'Best';
xlim([1e-1 inf])
ax = gca;
ax.FontSize = 14; 
xlabel('$H$ (kA/m)','Interpreter','latex','Fontsize',20)
ylabel('$-\Im(\overline{\mu}_r)$','Interpreter','latex','Fontsize',20)
%ylim([-inf 0])
%saveas(gcf,'MU_IMAG.png')

figure(3) 
%lgd = legend('H_0 = '+string(H0_list*1e-3)+{' '}+'kA/m');
%lgd.FontSize = 8;
%lgd.Location = 'southwest';
xlabel('H (kA/m)', 'Fontsize',16)
ylabel('M(T)', 'Fontsize',16)


figure
h = surf(Temp_grid*1e-3,H_grid*1e-3,mu_real_grid);
set(h,'edgecolor','none')
xlabel('H_0 (kA/m)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)
zlabel('Re(\mu)','FontSize', 14)

figure
h = surf(Temp_grid*1e-3,H_grid*1e-3,mu_im_grid);
set(h,'edgecolor','none')
xlabel('T (deg)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)
zlabel('Im(\mu)','FontSize', 14)

figure
h = pcolor(Temp_grid*1e-3,H_grid*1e-3,mu_real_grid);
set(h,'edgecolor','none')
colorbar
xlabel('T (deg)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)

figure
h = pcolor(Temp_grid*1e-3,H_grid*1e-3,mu_im_grid);
set(h,'edgecolor','none')
colorbar
xlabel('T (deg)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)


%% Enregistrement pour lookup table (changement de grille)
[Xq, Yq] = meshgrid(linspace(700,0,500),linspace(100E3,0,500));
Vq_real = griddata(Temp_export, H_export, mu_real_export, Xq, Yq);
Vq_imag = griddata(Temp_export, H_export, mu_im_export, Xq, Yq);

figure
h = surf(Xq*1e-3,Yq*1e-3,Vq_real);
set(h,'edgecolor','none')
xlabel('T (deg)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)
zlabel('Re(\mu)','FontSize', 14)

figure
h = surf(Xq*1e-3,Yq*1e-3,Vq_imag);
set(h,'edgecolor','none')
xlabel('T (deg)','FontSize', 14)
ylabel('H (kA/m)','FontSize', 14)
zlabel('Im(\mu)','FontSize', 14)

Temp_export = reshape(Xq',[],1);
H_export = reshape(Yq',[],1);
mu_real_export = reshape(Vq_real',[],1);
mu_im_export = reshape(Vq_imag',[],1);

mat_real = [Temp_export H_export mu_real_export];
mat_real(any(isnan(mat_real), 2), :) = [];

mat_im = [Temp_export H_export mu_im_export];
mat_im(any(isnan(mat_im), 2), :) = [];

%dlmwrite(strcat('table_mu_real.txt'), mat_real, '\t');
%dlmwrite(strcat('table_mu_imag.txt'), mat_im, '\t');
