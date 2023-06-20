clc
clear all
close all

mu0_const = 4*pi*1e-7;
currentFolder = pwd;
addpath('FEM_transient_fortran', 'FEM_PEM_1D', 'Material_computation')

%% Input parameters
freq          = 10e3;                                  %Hz (No influence on the permeability curves)
L             = 2.5e-3;                                %m (adapt L to have a long enough domain with penetration depth)
rho           = 2.5e-7;                                %Ohm/m (No influence on the permeability curves)
H0_list       = [1:4,5:5:100]*1e3;                     %A/m
mattype       = 7;                                     %1-7
output        = '\Results\';                           %Output folder
flag_real     = false;                                 %Flag if we want to approximate with only a real permeability
                                                       %(i.e. Pj+Ph -> Pj)
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
    param_phy.natan_list   = [1,2];
    param_phy.ai_list = {[3.6/pi],[3.6/pi,0.36/pi]};
    param_phy.ci_list = {[1.7e3],[1.7e3,1.7e3]};
    param_phy.bi_list = {[618.74],[618.74,618.74]};
    
    param_phy.Temp_list = [25,50];        %celsius (informative only)
    
% 4 : High field limit case (not used very much)
%     B = mu0*H+Bsat if H > Hsat
%     B = mu0*H-Bsat if H < -Hsat
%     3rd order polynomial between -Hsat and Hsat
elseif mattype == 4
    param_phy.Bsat_list = [1];
    param_phy.Hsat_list = [100];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 5: Preisach Model defining F and G functions from an elliptic major cycle
%    defined with a constant complex permeability value (not used very much)
elseif mattype == 5
    param_phy.mu_real_list = [100];
    param_phy.mu_imag_list = [-10];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 6: Preisach Model with elliptic hysteretic curves (minor and major, 
%    not optimal...)
elseif mattype == 6
    param_phy.mu_real_list = [100];
    param_phy.mu_imag_list = [-5];
    
    param_phy.Temp_list = [25];        %celsius (informative only)
    
% 7 : Preisach 4 parameters non linear hysteretic model
elseif mattype == 7
    
    densite = 7850;
    param_phy.Br_list = [0.933108696857510,...
        0.920305377456978,0.930197751762866,...
        0.926456606409189,0.907901897554444,...
        0.834068184565533,0.821744739128845,...
        0.770564680189743,0.751389871375226,...
        0.679519143367053,0.600291623357859,...
        0.546451764614761,0.447966383479637,...
        0.294911109267253,0.117203883430224]; %T
    
    param_phy.Bsat_list = [1.96514555000160,...
        1.89816938124426,1.80278279035741,...
        1.72179361157548,1.63349133258213,...
        1.54758208394500,1.50581060201378,...
        1.42710185702372,1.37893617396983,...
        1.26059864901288,1.13165570989077,...
        1.03681069040786,0.923484496007370,...
        0.779667342725369,0.569298099828659];  %T
    
    param_phy.Hc_list = [1.95434793735610,...
        1.99679351082121,2.00604344074037,...
        1.77397723667928,1.49198530528908,...
        1.24745259010741,1.17284459110538,...
        1.04743752016454,0.975494190801733,...
        0.828481584863761,0.681070834335290,...
        0.589536030141750,0.507417304279865,...
        0.420054615271398,0.309345454205599]*1e3; %A/m
    
    param_phy.Wh_list = [10.6002733081384,...
        10.6712910496742,10.5273882482361,...
        8.45674877296000,6.82871746205513,...
        5.44235966335567,5.01539958455657,...
        4.35243578913929,3.97342776473797,...
        3.29549446907246,2.63188834671565,...
        2.24470941769569,1.87330007080518,...
        1.50684313735103,1.15671952792303]*1e3/densite; %J/kg;
    
    param_phy.mur_max_list = [670.569648615508,...
        687.872892200044,725.034216103888,...
        824.954654313757,911.258096055197,...
        968.896684679610,970.918365865108,...
        949.275856912466,930.437384409310,...
        900.758765048889,889.189340732500,...
        859.018984882164,774.678574845519,...
        615.496262753714,295.014147745148]; %not necessary, will be calculated if missing
    

    param_phy.Temp_list = [24.7297379595883,...
        113.390619575718,246.013103219732,...
        409.207259467707,461.194556251831,...
        506.975677300830,522.429199176500,...
        547.161505490999,560.895854473397,...
        584.962374821997,607.172498517265,...
        620.125559533445,632.554505459268,...
        645.927268534192,663.035642782262]; %celsius (informative only)
    
    computePreisach4p(densite);
    
else
    error('Material type not implemented.')
end

%% 1-D Transient losses
slabProblem_IO(freq, L, rho, H0_list, mattype, currentFolder, output);

%% PEM
powerEquivalentModel(freq, L, rho, H0_list, mattype, currentFolder, output, flag_real);
