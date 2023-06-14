clc
clear all
close all

mu0_const = 4*pi*1e-7;
currentFolder = pwd;

%% Input parameters
freq          = 10e3;                                  %Hz (No influence on the permeability curves)
L             = 2.5e-3;                                %m (adapt L to have a long enough domain with penetration depth)
rho           = 2.5e-7;                                %Ohm/m (No influence on the permeability curves)
H0_list       = [1:4,5:5:100]*1e3;                     %A/m
mattype       = 7;                                     %1-7
output        = '\Results\';                           %Output folder

global param_phy

%% Material type
% 1 : constant relative permeability
%     B = mu0*mur*H
if mattype == 1
    param_phy.mur_list     = [0.10e+03];
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 2 : non-linear atan model
%     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
elseif mattype == 2
    param_phy.murmax_list  = [6.0e+02];
    param_phy.bsat_list    = [1.8e+00];
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 3 : Preisach 3n coefficients non linear hysteretic model
%     B+ = Sum ai*atan((H+ci)/bi)
%     B- = Sum ai*atan((H-ci)/bi)
elseif mattype == 3
    param_phy.natan_list   = [1];
    param_phy.ai_list = [3.6/pi];
    param_phy.ci_list = [1.7e3];
    param_phy.bi_list = param_phy.ci_list./(tan(1.4./param_phy.ai_list));
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 4 : High field limit case (not used very much)
%     B = mu0*H+Bsat if H > Hsat
%     B = mu0*H-Bsat if H < -Hsat
%     3rd order polynomial between -Hsat and Hsat
elseif mattype == 4
    param_phy.Bsat_list = [1];
    param_phy.Hsat_list = [20];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 5: Preisach Model defining F and G functions from an elliptic major cycle
%    defined with a constant complex permeability value (not used very much)
elseif mattype == 5
    param_phy.mu_real_list = [100];
    param_phy.mu_imag_list = [-50];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 6: Preisach Model with elliptic hysteretic curves (minor and major, 
%    not optimal...)
elseif mattype == 6
    param_phy.mu_real_list = [100];
    param_phy.mu_imag_list = [-5];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 7 : Preisach 4 parameters non linear hysteretic model
elseif mattype == 7
    %Usage of 4340 AISI data in MATLAB Figure
    fig = open('AnalyseThermoMagnetique_4340ferrite.fig');
    axObjs = fig.Children;
    
    %Wh
    dataObjs =  axObjs(1).Children;
    Temp1 = dataObjs(2).XData;
    Whdata = dataObjs(2).YData;
    
    %Hc
    dataObjs =  axObjs(2).Children;
    Temp2 = dataObjs(2).XData;
    Hcdata = dataObjs(2).YData;
    
    %Br
    dataObjs =  axObjs(3).Children;
    Temp3 = dataObjs(2).XData;
    Brdata = dataObjs(2).YData;
    
    %Bsat
    dataObjs =  axObjs(4).Children;
    Temp4 = dataObjs(2).XData;
    Bsatdata = dataObjs(2).YData;
    
    %murmax
    dataObjs =  axObjs(5).Children;
    Temp5 = dataObjs(2).XData;
    murmaxdata = dataObjs(2).YData;
    
    minlen = min([length(Temp1),length(Temp2),length(Temp3),length(Temp4),length(Temp5)]);

    
    param_phy.Br_list = Brdata(1:minlen);    %T
    param_phy.Bsat_list = Bsatdata(1:minlen);  %T
    param_phy.Hc_list = Hcdata(1:minlen); %A/m
    param_phy.Wh_list = Whdata(1:minlen);
    param_phy.mur_max_list = murmaxdata(1:minlen);
    param_phy.a_list = zeros(size(param_phy.Br_list));
    param_phy.s_list = zeros(size(param_phy.Br_list));
    param_phy.b_list = zeros(size(param_phy.Br_list));
    
    %Optimization of a, b and s parameters and convert Wh to s parameter
    %for the Preisach model
    for i = 1:length(param_phy.Br_list)
        densite = 7850;
        Br = param_phy.Br_list(i);
        Bsat = param_phy.Bsat_list(i);
        Hc = param_phy.Hc_list(i)*1000;
        Wh = param_phy.Wh_list(i)*1000/densite;
        mur_max = param_phy.mur_max_list(i);
        
        b = @(x) x(1).*(x(2)+sqrt((x(4)-x(3))./x(3)));
        Fplus = @(x)  (x(4)-x(3))*(x(5)./b(x)).*(1+(x(5)./b(x)).^(x(2)+1)).^(-1/(x(2)+1));
        Gplus = @(x) x(3)-x(3)*(1+(x(5)./x(1)).^(x(2)+2)).^(-1);
        
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
        
    end
    param_phy.Temp_list = Temp1(1:minlen);        %celsius (informative only)
    
else
    error('Material type not implemented.')
end

addpath('FEM_transient_fortran', 'FEM_PEM_1D')

%% 1-D Transient losses
slabProblem_IO(freq, L, rho, H0_list, mattype, currentFolder, output);

%% PEM
powerEquivalentModel(freq, L, rho, H0_list, mattype, currentFolder, output);
