function [] = computePreisach4p(densite)

global param_phy
mu0_const = 4*pi*1e-7;

%Optimization of a, b and s parameters and convert Wh to s parameter
%for the Preisach model

mur_flag = isfield(param_phy, 'mur_max_list'); 

for i = 1:length(param_phy.Br_list)
    Br = param_phy.Br_list(i);
    Bsat = param_phy.Bsat_list(i);
    Hc = param_phy.Hc_list(i);
    Wh = param_phy.Wh_list(i);
    
    
    b = @(x) x(1).*(x(2)+sqrt((x(4)-x(3))./x(3)));
    Fplus = @(x)  (x(4)-x(3))*(x(5)./b(x)).*(1+(x(5)./b(x)).^(x(2)+1)).^(-1/(x(2)+1));
    Gplus = @(x) x(3)-x(3)*(1+(x(5)./x(1)).^(x(2)+2)).^(-1);
    
    if mur_flag      
        mur_max = param_phy.mur_max_list(i);
        eq1 = @(x) mu0_const*x(5)+Fplus(x)+2*Gplus(x)-x(3);
        eq2 = @(x) 4*pi*x(1)*x(3)/(densite*(x(2)+2)*sin(pi/(x(2)+2)))-Wh;
        eq3 = @(x) 1+((x(4)-x(3))/mu0_const).*((1./b(x)).*(1+(x(5)./b(x)).^(x(2)+1)).^(-1./(x(2)+1))+(x(5)./(b(x)))...
            .*(-1./(x(2)+1)).*(1+(x(5)./b(x)).^(x(2)+1)).^((-1./(x(2)+1))-1).*((x(2)+1)./b(x)).*(x(5)./b(x)).^x(2))...
            + (2*x(3)/mu0_const)*(1+(x(5)./x(1)).^(x(2)+2)).^(-2).*((x(2)+2)./x(1)).*(x(5)./x(1)).^(x(2)+1)-mur_max;
        eq4 = @(x) x(3)-Br;
        eq5 = @(x) x(4)-Bsat;
        eq6 = @(x) x(5)-Hc;
        
        options = optimoptions('lsqnonlin','Display','iter');
        options.FunctionTolerance = 1e-18;
        options.OptimalityTolerance = 1e-18;
        options.StepTolerance = 1e-18;
        options.MaxFunctionEvaluations=5000;
        F = @(x) [eq1(x),eq2(x),eq3(x),eq4(x),eq5(x),eq6(x)];
        lb = [0,0,0,0,0];
        ub = [Inf,10,Inf,Inf,Inf];
        x0 = [1e4,5,Br,Bsat,Hc];
        
        [x,fval] = lsqnonlin(F,x0,lb,ub,options);
        param_phy.a_list(i) = x(1);
        param_phy.s_list(i) = x(2);
        param_phy.b_list(i) = x(1).*(x(2)+sqrt((Bsat-Br)./Br));
        param_phy.Br_list(i) = x(3);
        param_phy.Bsat_list(i) = x(4);
        param_phy.Hc_list(i) = x(5);
        
    else 
        eq1 = @(x) mu0_const*x(5)+Fplus(x)+2*Gplus(x)-x(3);
        eq2 = @(x) 4*pi*x(1)*x(3)/(densite*(x(2)+2)*sin(pi/(x(2)+2)))-Wh;
        eq3 = @(x) x(3)-Br;
        eq4 = @(x) x(4)-Bsat;
        eq5 = @(x) x(5)-Hc;
        
        options = optimoptions('lsqnonlin','Display','iter');
        options.FunctionTolerance = 1e-18;
        options.OptimalityTolerance = 1e-18;
        options.StepTolerance = 1e-18;
        options.MaxFunctionEvaluations=5000;
        F = @(x) [eq1(x),eq2(x),eq3(x),eq4(x),eq5(x)];
        lb = [0,0,0,0,0];
        ub = [Inf,10,Inf,Inf,Inf];
        x0 = [1e4,5,Br,Bsat,Hc];
        
        [x,fval] = lsqnonlin(F,x0,lb,ub,options);
        param_phy.a_list(i) = x(1);
        param_phy.s_list(i) = x(2);
        param_phy.b_list(i) = x(1).*(x(2)+sqrt((Bsat-Br)./Br));
        param_phy.Br_list(i) = x(3);
        param_phy.Bsat_list(i) = x(4);
        param_phy.Hc_list(i) = x(5);
        param_phy.mur_max_list(i) = 1+((x(4)-x(3))/mu0_const).*((1./param_phy.b_list(i))...
            .*(1+(x(5)./param_phy.b_list(i)).^(x(2)+1)).^(-1./(x(2)+1))+(x(5)./(param_phy.b_list(i)))...
            .*(-1./(x(2)+1)).*(1+(x(5)./param_phy.b_list(i)).^(x(2)+1)).^((-1./(x(2)+1))-1).*...
            ((x(2)+1)./param_phy.b_list(i)).*(x(5)./param_phy.b_list(i)).^x(2))...
            + (2*x(3)/mu0_const)*(1+(x(5)./x(1)).^(x(2)+2)).^(-2).*((x(2)+2)./x(1)).*(x(5)./x(1)).^(x(2)+1);
    end
    
end