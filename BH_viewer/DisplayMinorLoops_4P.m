function DisplayMinorLoops_4P(Hc,Br,Bs,s)
% Displays a series of concentric minors hysteresis loops associated with the
% parameters Hc, Br, Bs and s of the Preisach model with 4 parameters.

% Calcul des parametres internes a et b
[a,b] = nra(Hc,Br,Bs,s);

% H saturation
Hsat = Hc + 10*b;

% Fonctions nonlineaires
s1 = s+1;
F = @(h) (Bs-Br).*(h/b)./(1+abs(h/b).^s1).^(1/s1);
G = @(h) sign(h).*(Br-Br./(1+abs(h/a).^(1+s1)));

% Boucle pour l'affichage
nbCycles = 10;

figure(1)
plot(nan,nan);
hold on

for i = 1:nbCycles
  
  % B(H=0)
  B0 = Br*(i-0.5)/nbCycles;
  
  % Point maximal sur la courbe de premiere aimantation
  G1 = @(h) sign(h).*G(h).^2.0/Br-B0;
  H1 = fzero(G1,[0,Hsat]);
  
  % Affichage H < 0
  Hi  = linspace(-H1,H1,201);
  Bdi = [2.0*G(Hi(1:100))*G(H1)/Br,zeros(1,101)] + B0 + F(Hi);
  Bai = [zeros(1,100),2.0*G(Hi(101:201))*G(H1)/Br] - B0 + F(Hi);
  plot(Hi,Bdi,'b',Hi,Bai,'b')

end

hold off
grid on

% axes
xlim(1.1*[-H1,H1])
xlabel('H (A/m)')
ylabel('B (T)')