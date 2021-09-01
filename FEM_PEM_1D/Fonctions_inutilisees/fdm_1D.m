function [H, dH, ddH, x] = fdm_1D( H0, dH0, xvec, P_exp )

% Test de resolution par deux méthodes rk4 si on a H0 et dH0 (fonction pas
% utilisee pour le moment)

global param_phy

u0 = H0*dH0;
x0 = xvec(1);
nmax = length(xvec)-1;
h = xvec(2)-xvec(1);
nbeq = size(u0,2);

u = zeros(1,nbeq);
n = 1;
x(1) = x0;
u(1,:) = u0';

while(n <= nmax)
    
    k1 = h * P_exp(n)*2/param_phy.rho;
    k2 = h * (P_exp(n)+P_exp(n+1))/param_phy.rho;
    k3 = k2;
    k4 = h * P_exp(n+1)*2/param_phy.rho;
    
    u0 = u0 + (1/6) * (k1+2*(k2+k3) + k4);
    u(n+1,:) = u0';
    x(n+1) = x(n) + h;
    n = n+1;
end

H = zeros(1,nbeq);
n = 1;
H(1,:) = H0';
while(n <= nmax)
    
    k1 = h * u(n)/H0;
    k2 = h * ((u(n)+u(n+1))/2)/(H0+k1/2);
    k3 = h * ((u(n)+u(n+1))/2)/(H0+k2/2);
    k4 = h * u(n+1)/(H0+k3);
    
    H0 = H0 + (1/6) * (k1+2*(k2+k3) + k4);
    H(n+1,:) = H0';
    n = n+1;
end

[dH, ddH] = diffO2(x,H);
x = x';