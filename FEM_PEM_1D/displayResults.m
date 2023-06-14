function displayResults(xvec_init, P_Joule, P_Hyst, Hk, xvec, H, ...
                            dH, ddH, P_exp, P, mu_real, mu_im, ...
                            Pjoule_approx, Physt_approx, Ptot_approx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% afficher_resultats
% 
% Affiche les resultats du probleme PEM 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
figure
plot(xvec_init*1e3,P_Joule*(100^(-3)),'Linewidth',2)
hold on
plot(xvec_init*1e3,P_Hyst*(100^(-3)),'Linewidth',2)
xlabel('x (mm)', 'Fontsize',18)
ylabel('Pertes (W/cm^3)', 'Fontsize',18)
lgd = legend('Pertes par courants de Foucault','Pertes par hysteresis');
lgd.FontSize = 10;

figure
plot(xvec_init*1e3,Hk,'Linewidth', 2)
hold on
plot(xvec*1e3,H,'Linewidth', 2)
xlabel('x (mm)', 'Fontsize',18)
ylabel('H (A/m)', 'Fontsize',18)
lgd = legend('Approximation initiale','Solution');
lgd.FontSize = 10;

figure
plot(xvec*1e3,dH,'Linewidth', 2)
xlabel('x (mm)', 'Fontsize',18)
ylabel('dH/dx (A/m^2)', 'Fontsize',18)

figure
plot(xvec*1e3,ddH,'Linewidth', 2)
xlabel('x (mm)', 'Fontsize',18)
ylabel('d^2H/dx^2 (A/m^3)', 'Fontsize',18)

figure
plot(xvec_init*1e3,P_exp*(100^(-3)), 'Linewidth', 2)
hold on
plot(xvec*1e3,P*(100^(-3)), '--', 'Linewidth', 2);
lgd=legend('Pertes calculees dans le domaine temporel','Pertes calculees par le PEM');
lgd.FontSize = 10;
xlabel('x (mm)', 'Fontsize',18)
ylabel('Pertes (W/cm^3)', 'Fontsize',18)

figure
plot(H*1e-3, mu_real, 'Linewidth', 2)
xlabel('H (kA/m)', 'Fontsize',18)
ylabel('Re(\mu_r)', 'Fontsize',18)
xlim([min(H)*1e-3, max(H)*1e-3])

figure
plot(H*1e-3,mu_im, 'Linewidth', 2)
xlabel('H (kA/m)', 'Fontsize',18)
ylabel('Im(\mu_r)', 'Fontsize',18)
xlim([min(H)*1e-3, max(H)*1e-3])

figure
semilogx(H*1e-3,mu_real, 'Linewidth', 2)
hold on
semilogx(H*1e-3,-mu_im, 'Linewidth', 2)
lgd = legend('Re(\mu_r)', 'Im(-\mu_r)');
lgd.FontSize = 10;
xlabel('H (kA/m)', 'Fontsize',18)
ylabel('\mu_r', 'Fontsize',18)
xlim([min(H)*1e-3, max(H)*1e-3])
ylim([0 inf])

figure
plot(H*1e-3,(mu_real-1).*H*4*pi*1e-7,'Linewidth',2)
xlabel('H (kA/m)', 'Fontsize',18)
ylabel('M(T)', 'Fontsize',18)

x0=400;
y0=150;
width=800;
height=625;

figure
set(gcf,'position',[x0,y0,width,height])
plot(xvec_init*1e3,P_Joule*(100^(-3)),'Linewidth',3)
hold on
plot(xvec*1e3,Pjoule_approx*(100^(-3)),'--','Linewidth',3)
lgd = legend('Régime temporel','Régime harmonique');
lgd.FontSize = 22;
lgd.Location = 'northeast';
ax = gca;
ax.FontSize = 18;
xlabel('Profondeur (mm)', 'Fontsize',20)
ylabel('Pertes par courants de Foucault (W/cm^3)', 'Fontsize',20)

figure
set(gcf,'position',[x0,y0,width,height])
plot(xvec_init*1e3,P_Hyst*(100^(-3)),'Linewidth',3)
hold on
plot(xvec*1e3,Physt_approx*(100^(-3)),'--','Linewidth',3)
lgd = legend('Régime temporel','Régime harmonique');
lgd.FontSize = 22;
lgd.Location = 'northeast';
ax = gca;
ax.FontSize = 18; 
xlabel('Profondeur (mm)', 'Fontsize',20)
ylabel('Pertes par hysteresis (W/cm^3)', 'Fontsize',20)

figure
set(gcf,'position',[x0,y0,width,height])
plot(xvec_init*1e3,(P_Hyst+P_Joule)*(100^(-3)),'Linewidth',3)
hold on
plot(xvec*1e3,(Ptot_approx)*(100^(-3)),'--','Linewidth',3)
xlabel('x (mm)', 'Fontsize',28)
ylabel('Pertes totales (W/cm^3)', 'Fontsize',28)
lgd = legend('Pertes dans le domaine temporel','Pertes dans le domaine harmonique');
lgd.FontSize = 18;
lgd.Location = 'northeast';
ax = gca;
ax.FontSize = 22; 