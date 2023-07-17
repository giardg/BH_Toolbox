function [N,dN, ddN]=shape(xi,pDeg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shape
% 
% Retourne la valeur des fonctions de forme et de leur derivee en 1D au
% point xi sur l'element de reference ([-1,1]). Seuls les polynomes de
% Lagrange de degre 1, 2 et 3 ont ete implementes pour le moment (issu du
% cours MTH8207 A2019) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pDeg==1  %linear (only one form)
    N(1)=.5*(1-xi);
    N(2)=.5*(1+xi);
    dN(1)=-.5;
    dN(2)=.5;
    ddN(1) = 0;
    ddN(2) = 0;
end
if pDeg==2  %quadratic
    %if pType==1   %Lagrangian
        N(1)=.5*xi*(xi-1.);
        N(3)=1-xi^2;
        N(2)=.5*xi*(xi+1.);
        dN(1)=xi-.5;
        dN(3)=-2*xi;
        dN(2)=xi+.5;
        ddN(1) = 1;
        ddN(3) = -2;
        ddN(2) = 1;
    %end
    %if pType(ne)==2;   %hierarchical
    %    N(1)=.5*(1-xi);
    %    N(2)=.5*(1+xi);
    %    N(3)=1-xi^2;
    %    dN(1)=-.5;
    %    dN(2)=.5;
    %    dN(3)=-2*xi;
    %end
end

if pDeg==3
    N(1) = -(9/16)*(xi+1/3)*(xi-1/3)*(xi-1);
    N(2) = (9/16)*(xi+1/3)*(xi-1/3)*(xi+1);
    N(3) = (27/16)*(xi+1)*(xi-1/3)*(xi-1);
    N(4) = -(27/16)*(xi+1)*(xi+1/3)*(xi-1);
    
    dN(1) = -(9/16)*((xi-1/3)*(xi-1)+(xi+1/3)*(xi-1)+(xi+1/3)*(xi-1/3));
    dN(2) = (9/16)*((xi-1/3)*(xi+1)+(xi+1/3)*(xi+1)+(xi+1/3)*(xi-1/3));
    dN(3) = (27/16)*((xi-1/3)*(xi-1)+(xi+1)*(xi-1)+(xi+1)*(xi-1/3));
    dN(4) = -(27/16)*((xi+1/3)*(xi-1)+(xi+1)*(xi-1)+(xi+1)*(xi+1/3));

end

end

    