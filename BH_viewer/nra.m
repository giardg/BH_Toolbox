function [a,b] = nra(Hc,Br,Bs,s)
  
  % constantes
  lambda = sqrt((Bs-Br)/Br);
  
  % parametres N-R
  tol = 1e-9;
  maxIter = 25;
  
  % initialisation
  a = 1.1*Hc;
  da = Hc;
  iter = 0;
  
  % boucle Newton-Raphson
  while abs(da)/Hc > tol && iter < maxIter
    
    % incrementation
    iter = iter + 1;
    
    % iteration
    R  = Br + (Bs-Br)*Hc/(a*(s+lambda))/(1+(Hc/(a*(s+lambda)))^(1+s))^(1/(1+s))...
            - 2*Br/(1+(Hc/a)^(2+s));
    dR = (Bs-Br)/(s+lambda)/(1+(Hc/(a*(s+lambda)))^(1+s))^(1+1/(1+s))...
         + 2*Br*(2+s)*(Hc/a)^(1+s)/(1+(Hc/a)^(2+s))^2;
    da = a^2/Hc*R/dR;
    a  = a + da;
  
  end
  
  % algorithme de bissection si N-R echoue
  if abs(da)/Hc > tol
    
    warning('Newton failed!')
    
    % initialisation des bornes
    a1 = (1+eps)*Hc;
    R1 = Br + (Bs-Br)*Hc/(a1*(s+lambda))/(1+(Hc/(a1*(s+lambda)))^(1+s))^(1/(1+s))...
            - 2*Br/(1+(Hc/a1)^(2+s));
    p2 = -4;  
    a2 = (1+pow2(p2))*Hc;
    R2 = Br + (Bs-Br)*Hc/(a2*(s+lambda))/(1+(Hc/(a2*(s+lambda)))^(1+s))^(1/(1+s))...
            - 2*Br/(1+(Hc/a2)^(2+s));
    while R1*R2 >= 0
      p2 = p2 + 1;
      a2 = (1+pow2(p2))*Hc;
      R2 = Br + (Bs-Br)*Hc/(a2*(s+lambda))/(1+(Hc/(a2*(s+lambda)))^(1+s))^(1/(1+s))...
            - 2*Br/(1+(Hc/a2)^(2+s));
    end
    da = a2-a1;
    
    % verification a1 ou a2 = 0
    if R1 == 0
      a = a1;
      return
    elseif R2 == 0
      a = a2;
      return
    end
    
    % boucle bissection
    while abs(da)/Hc > tol
      am = 0.5*(a1+a2);
      Rm = Br + (Bs-Br)*Hc/(am*(s+lambda))/(1+(Hc/(am*(s+lambda)))^(1+s))^(1/(1+s))...
            - 2*Br/(1+(Hc/am)^(2+s));
      if Rm == 0
        a = am;
        return
      elseif R1*Rm < 0
        a2 = am;
        R2 = Rm;
      elseif Rm*R2 < 0
        a1 = am;
        R1 = Rm;
      end
      da = a2-a1;
    end
    a = 0.5*(a1+a2);
    
  end
  
  % Calcul de b
  b = a*(s+lambda);
  
end