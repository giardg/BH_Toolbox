function DisplayMinorLoops_3N(a,b,c)
% Displays a series of concentric minors hysteresis loops associated with the
% coefficients ai,bi,ci of the Preisach model with 3*n coefficients.

% Nombre de triplets
n = length(a);

if     n==1
    Bd = @(h) a*atan((h+c)/b);
elseif n==2
    Bd = @(h) a(1)*atan((h+c(1))/b(1))  + ...
              a(2)*atan((h+c(2))/b(2))  ;
elseif n==3
    Bd = @(h) a(1)*atan((h+c(1))/b(1))  + ...
              a(2)*atan((h+c(2))/b(2))  + ...
              a(3)*atan((h+c(3))/b(3))  ;
elseif n==4
    Bd = @(h) a(1)*atan((h+c(1))/b(1))  + ...
              a(2)*atan((h+c(2))/b(2))  + ...
              a(3)*atan((h+c(3))/b(3))  + ...
              a(4)*atan((h+c(4))/b(4))  ;
              
elseif n==5
    Bd = @(h) a(1)*atan((h+c(1))/b(1))  + ...
              a(2)*atan((h+c(2))/b(2))  + ...
              a(3)*atan((h+c(3))/b(3))  + ...
              a(4)*atan((h+c(4))/b(4))  + ...
              a(5)*atan((h+c(5))/b(5))  ;
else
    error('n must be an integer between 1 and 5.')
end

% B remanent et H saturation
Br = Bd(0);
Hsat = 10*max(b)+max(c);

% Fonctions nonlineaires
F = @(h) sign(h).*(Bd(abs(h))-Br);
G = @(h) sign(h).*(Br-0.5*(Bd(abs(h))+Bd(-abs(h))));

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

end
