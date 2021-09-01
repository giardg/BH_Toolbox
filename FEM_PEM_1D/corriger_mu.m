function [mu_real_out,mu_im_out] = corriger_mu(mu_real_in, mu_im_in, H, P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corriger_mu
% 
% Applique une correction sur la partie reelle de la permeabilite a faibles
% champs selon le type de materiau et ce qui est attendu theoriquement. On
% interpole les valeurs de permeabilite a faibles champs avec:
% 1) un polynome de degre 3 dans le cas lineaire (a partir de Pmax*2.5%)
% 2) un polynome de degre 3 dans le cas du modele avec arctangente (a 
%    partir de Pmax*2.5%)
% 3) un polynome de degre 3 dans le cas du hysteretique a partir de mu_max
%    pour H > Hc (si H0 < Hc, on coupe simplement a partir de mu_max )
%
% Un lissage de la courbe est aussi effectue.
% 
% Aout 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4*pi*1E-7;
mu_im_out = mu_im_in;
global param_phy

if param_phy.mattype == 1
    
    number_removed = sum(P < 0.025*max(P));
    mu_real_new = mu_real_in(1:end-number_removed);
    H_new = H(1:end-number_removed);
    H_removed = H(end-number_removed+1:end);
    mu_real_new = smooth(H_new,mu_real_new);
    
    mu_H0 = param_phy.mur;
    x2 = H_new(end);
    y2 = mu_real_new(end);
    dy2 = (mu_real_new(end)-mu_real_new(end-1))/(H_new(end)-H_new(end-1));
    
    mat = [3*x2^2,2*x2;x2^3,x2^2];
    b = [dy2;y2-mu_H0];
    coeff = mat\b;
    
    syms x
    fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2+mu_H0;
    mu_real_out = [mu_real_new;fct(H_removed)];
    
elseif param_phy.mattype == 2
    
    number_removed = sum(P < 0.025*max(P));
    mu_real_new = mu_real_in(1:end-number_removed);
    H_new = H(1:end-number_removed);
    H_removed = H(end-number_removed+1:end);
    mu_real_new = smooth(H_new,mu_real_new);
    
    mu_H0 = param_phy.murmax;
    
    x2 = H_new(end);
    y2 = mu_real_new(end);
    dy2 = (mu_real_new(end)-mu_real_new(end-1))/(H_new(end)-H_new(end-1));
    
    mat = [3*x2^2,2*x2;x2^3,x2^2];
    b = [dy2;y2-mu_H0];
    coeff = mat\b;
    
    syms x
    fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2+mu_H0;
    mu_real_out = [mu_real_new;fct(H_removed)];
    
elseif param_phy.mattype == 3
    
    if max(H) > param_phy.ci
        number_removed = sum(H < param_phy.ci);
    else
        number_removed = 0;
    end
    
    mu_real_new = mu_real_in(1:end-number_removed);
    
   
    
    [val,loc] = max(mu_real_new);
    number_removed = sum(H < H(loc));
   
    
    mu_real_new = mu_real_in(1:end-number_removed);
    H_new = H(1:end-number_removed);
    H_removed = [H(end-number_removed+1:end)];
    
    y1  = mu_real_new(end);
    x2 = H_new(end);
    y2 = mu_real_new(end);
    dy2 = 0;
    
    mat = [3*x2^2,2*x2;x2^3,x2^2];
    b = [dy2;y1];
    coeff = mat\b;
    
    syms x
    fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2;
    mu_real_out = [mu_real_new;fct(H_removed)];
    mu_real_out = smooth(H,mu_real_out);
   
%     mat = [3*x2^2,2*x2;x2^3,x2^2];
%     b = [dy2;y2-y1];
%     coeff = mat\b;
%     
%     syms x
%     fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2+y1;
%     mu_real_out = [mu_real_new;fct(H_removed)];
%     mu_real_out = smooth(H,mu_real_out);

    
elseif param_phy.mattype == 4
    number_removed = sum(P < 0.025*max(P));
    mu_real_new = mu_real_in(1:end-number_removed);
    H_new = H(1:end-number_removed);
    H_removed = H(end-number_removed+1:end);
    
    mu_H0 = (mu0+(3.0*param_phy.Bsat)/(2.0*param_phy.Hsat))/mu0;
    
    mu_real_new(end+1) = mu_H0;
    H_new(end+1) = 0;
    mu_real_out = spline(H_new,mu_real_new,H);
    
elseif param_phy.mattype == 6
    
    number_removed = sum(P < 0.025*max(P));
    mu_real_new = mu_real_in(1:end-number_removed);
    mu_im_new = mu_im_in(1:end-number_removed);
    H_new = H(1:end-number_removed);
    H_removed = H(end-number_removed+1:end);
    mu_real_new = smooth(H_new,mu_real_new);
    mu_im_new = smooth(H_new,mu_im_new);
    
    mu_H0 = param_phy.mur_real;
    x2 = H_new(end);
    y2 = mu_real_new(end);
    dy2 = (mu_real_new(end)-mu_real_new(end-1))/(H_new(end)-H_new(end-1));
    
    mat = [3*x2^2,2*x2;x2^3,x2^2];
    b = [dy2;y2-mu_H0];
    coeff = mat\b;
    
    syms x
    fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2+mu_H0;
    mu_real_out = [mu_real_new;fct(H_removed)];
    
    mu_H0 = param_phy.mur_im;
    x2 = H_new(end);
    y2 = mu_im_new(end);
    dy2 = (mu_im_new(end)-mu_im_new(end-1))/(H_new(end)-H_new(end-1));
    
    mat = [3*x2^2,2*x2;x2^3,x2^2];
    b = [dy2;y2-mu_H0];
    coeff = mat\b;
    
    syms x
    fct=@(x) coeff(1)*x.^3+coeff(2)*x.^2+mu_H0;
    mu_im_out = [mu_im_new;fct(H_removed)];
end


