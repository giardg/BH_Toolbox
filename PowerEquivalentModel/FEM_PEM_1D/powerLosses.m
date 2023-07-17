function [Ptot,Pjoule,Physt] = powerLosses(PQ,H,mu)

[dabsHdz,d2absHdz2] = diffO2(PQ.x,abs(H));
Ptot = 0.5*PQ.rho*(abs(H).*d2absHdz2+dabsHdz.^2);

if nargin == 3
    dHdz = diffO2(PQ.x,H);
    Pjoule = 0.5*PQ.rho*abs(dHdz).^2;
    Physt  = -0.5*PQ.omega*imag(mu).*abs(H).^2;
    
else
    Pjoule = nan(size(H));
    Physt  = nan(size(H));
end


end