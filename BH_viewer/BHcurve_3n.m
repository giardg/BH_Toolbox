function Bt = BHcurve_3n(Ht,a,b,c)
% Calculate the B(t) curve from the H(t) curve and the parameters Hc, Br,
% Bs and s of the Preisach model with 3n coefficients.

mu0_const = 4*pi*1e-7;
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

E_FG = @(alpha,beta) (F(alpha)-F(beta))/2 - (1/Br)*(alpha*beta<0)*G(alpha)*G(beta);

Bt = zeros(size(Ht));
for i = 1:length(Ht)
   H = Ht(i);
   M = PreisachOp(E_FG,H);
   Bt(i) = mu0_const*H+M;
end
