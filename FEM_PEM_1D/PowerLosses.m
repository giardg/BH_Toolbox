function [Ptot,Pjoule,Physt,delta] = powerLosses(PQ,H,mu)

[dabsHdz,d2absHdz2] = diffO2(PQ.x,abs(H));
Ptot = 0.5*PQ.rho*(abs(H).*d2absHdz2+dabsHdz.^2);

if nargin == 3
    dHdz = diffO2(PQ.x,H);
    Pjoule = 0.5*PQ.rho*abs(dHdz).^2;
    Physt  = -0.5*PQ.omega*imag(mu).*abs(H).^2;
    
    % Longueur de penetration effective
    surfJ = abs(dHdz(1)); %Courant maximal à la surface
    flag_skindepth = false;
    i = 1;
    while ~flag_skindepth
        if i == length(H)
            flag_skindepth = true;
            delta = Inf;
        elseif abs(dHdz(i)) <= surfJ*exp(-1)
            flag_skindepth = true;
            delta = PQ.x(i);
        end
        i = i+1;
    end
else
    Pjoule = nan(size(H));
    Physt  = nan(size(H));
end


end