function Bt = BHcurve_4P(Ht,Hc,Br,Bs,s)
% Calculate the B(t) curve from the H(t) curve and the parameters Hc, Br,
% Bs and s of the Preisach model with 4 parameters.

mu0_const = 4*pi*1e-7;
% Calcul des parametres internes a et b
[a,b] = nra(Hc,Br,Bs,s);

% H saturation
Hsat = Hc + 10*b;

% Fonctions nonlineaires
s1 = s+1;
F = @(h) (Bs-Br).*(h/b)./(1+abs(h/b).^s1).^(1/s1);
G = @(h) sign(h).*(Br-Br./(1+abs(h/a).^(1+s1)));

E_FG = @(alpha,beta) (F(alpha)-F(beta))/2 - (1/Br)*(alpha*beta<0)*G(alpha)*G(beta);

Bt = zeros(size(Ht));
for i = 1:length(Ht)
   H = Ht(i);
   M = PreisachOp(E_FG,H);
   Bt(i) = mu0_const*H+M;
end

