function [] = powerEquivalentModel(freq, L, rho, H0_list, mattype, currentFolder, output, flag_real)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% powerEquivalentModel.m
%
% Main function to solve the Power-Equivalent Model in 1-D and determine mu
% from power data
%
% Steps:  1) Solve P = (rho/2)*[H*H'' + (H')^2] in fem_1D.m
%         2) Compute phi' = -(1/H)*sqrt((2/rho)*P_Joule - (H')^2) and phi''
%         3) Compute de mu_real and mu_imag
%
% July 2020
% Gregory Giard, Polytechnique Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Numerical parameters
global param_num

param_num.nel           = 200; %number of elements
param_num.nbIter        = 50; %maximum number of iteration (Newton)
param_num.critere       = 1e-4; %stop criteria (Newton) |delta|/|Hk|
param_num.deg           = 2; %Basis functions order (Lagrange degre 1, 2 et 3 implemented)
param_num.nquad         = 4; %Number of Gauss points
param_num.phi0          = pi/2; %phase at x = 0
param_num.sparsity      = false;
param_num.amortissement = .1; %damping factor (slower but steady convergence: default = 1)
param_num.correction    = true; %Apply correction on mu (correction and smoothing for low fields)
display                 = true; %Display results (usually turn on if we have only 1 equivalent curve to compute)
CF                      = [1,2]; %Boundary conditions on both edge (1 = Dirichlet, 2 = Neumann)


%% Physical parameters
mu0 = 4*pi*1E-7;
global param_phy
param_phy.rho           = rho;
param_phy.mattype       = mattype; 
param_phy.longueur      = L; %domain length
param_phy.omega         = 2*pi*freq; %frequency

%Loop over H0
for H0 = H0_list
    param_phy.H0 = H0;
    
    %Loop over the temperatures (magntetic measurments)
    for i = 1:length(param_phy.Temp_list)
        
        close all
        param_phy.Temp = param_phy.Temp_list(i);
        
        % 1 : constant relative permeability
        %     B = mu0*mur*H
        if param_phy.mattype == 1
            typename = 'Linear';
            param_phy.mur     = param_phy.mur_list(i);
            
        % 2 : non-linear atan model
        %     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
        elseif param_phy.mattype == 2
            typename = 'Anhyst';
            param_phy.murmax  = param_phy.murmax_list(i);
            param_phy.Bsat = param_phy.bsat_list(i);
            
        % 3 : Preisach 3n coefficients non linear hysteretic model
        %     B+ = Sum ai*atan((H+ci)/bi)
        %     B- = Sum ai*atan((H-ci)/bi)
        elseif param_phy.mattype == 3
            typename = 'Hysteresis';
            param_phy.natan   = param_phy.natan_list(i);
            param_phy.ai = param_phy.ai_list{i};
            param_phy.bi = param_phy.bi_list{i};
            param_phy.ci = param_phy.ci_list{i};
            
        % 4 : High field limit case (not used very much)
        %     B = mu0*H+Bsat if H > Hsat
        %     B = mu0*H-Bsat if H < -Hsat
        %     3rd order polynomial between -Hsat and Hsat
        elseif param_phy.mattype == 4
            typename = 'Limit';
            param_phy.Bsat = param_phy.Bsat_list(i);
            param_phy.Hsat = param_phy.Hsat_list(i);
            
        % 5: Preisach Model defining F and G functions from an elliptic major cycle
        %    defined with a constant complex permeability value (not used very much)
        elseif param_phy.mattype == 5
            typename = 'Ellipse';
            param_phy.mur_real     = param_phy.mur_real_list(i);
            param_phy.mur_im       = param_phy.mur_im_list(i);
        
        % 6: Preisach Model with elliptic hysteretic curves (minor and major, 
        %    not optimal...)
        elseif param_phy.mattype == 6
            typename = 'Ellipse2';
            param_phy.mur_real     = param_phy.mur_real_list(i);
            param_phy.mur_im       = param_phy.mur_im_list(i);
        
        % 7 : Preisach 4 parameters non linear hysteretic model
        elseif param_phy.mattype == 7
            typename = 'Hysteresis2';
            param_phy.Hc = param_phy.Hc_list(i);
            
        else
            error('Material type not implemented.')
        end
        
        %% Get P_Joule and P_Hyst
        
        %addpath('Resultats_transitoires_fortran','Resultats_transitoires_COMSOL','Fonctions_FEM','Fonctions_phi')
        addpath([currentFolder output '\Transient_results'],[currentFolder '\FEM_PEM_1D\Functions_FEM'],[currentFolder '\FEM_PEM_1D\Functions_phi'])
        
        
        % Get data in files
        [xvec_init, P_Joule, P_Hyst] = readLosses(freq, currentFolder, output);
        if flag_real %We replace the joule loss by total loss to make the permeability purely real
            P_Joule = P_Joule+P_Hyst;
            P_Hyst = 0.*P_Hyst;
        end
        P_exp = P_Joule+P_Hyst;
        
        % Initial approximation of the solution (must respect the BC)
        syms x
        Hk=@(x) param_phy.H0+0*x;
        dHk=@(x) 0*x;
        ddHk=@(x) 0*x;
        
        Hk = Hk(xvec_init);
        dHk = dHk(xvec_init);
        ddHk = ddHk(xvec_init);
        
        
        %% Finite element problem to obtain H (Step 1)
        
        % Solving the finite element problem
        [H, dH, ddH, xvec] = fem1D( CF, xvec_init, Hk, dHk, P_exp );
        P = (param_phy.rho/2)*((H).*ddH+dH.^2);
        
        %% Compute the phi equation (Step 2)
        dphi = (-1./H).*sqrt((2/param_phy.rho)*spline(xvec_init,P_Joule,xvec)-dH.^2);
        ddphi = diffO2(xvec,dphi);
        
        %% Compute mu' and mu'' (Step 3)
        mu_real = (1/mu0)*(param_phy.rho/param_phy.omega)*((2./H).*dH.*dphi+ddphi);
        mu_im = (1/mu0)*(-2*spline(xvec_init,P_Hyst,xvec)/param_phy.omega)./(H.^2);
        if (param_num.correction) [mu_real,mu_im] = correctMu(mu_real, mu_im, H, P); end %Correction de mu à faible champ
        mu = mu0*(mu_real+1i*mu_im);
        
        %% Verification (FoucaultSolve.m)
        PQ = struct('rho',param_phy.rho,'omega',param_phy.omega,'x',xvec,'dx',xvec(2)-xvec(1),'nx',length(xvec),'Hmax', param_phy.H0);
        H2 = foucaultSolve( PQ, transpose(mu) );
        [Ptot_approx,Pjoule_approx,Physt_approx] = powerLosses(PQ,H2,mu);
        
        %% Display results
        if display
            displayResults(xvec_init, P_Joule, P_Hyst, Hk, xvec, H, ...
                dH, ddH, P_exp, P, mu_real, mu_im, ...
                Pjoule_approx, Physt_approx, Ptot_approx)
        end
        
        %% Save results
        filename = strcat('H0',{' '},string(param_phy.H0),'_mu',typename,'_',num2str(param_phy.Temp),'deg_',num2str(freq/1e3),'kHz');
        if param_num.correction filename = strcat(filename,'_corrected'); end
        if flag_real filename = strcat(filename,'_real'); end
        filename = strcat(filename,'.txt');
        if flag_real
            dlmwrite(strcat(currentFolder, output, 'Mu_results\', filename), [H mu_real], '\t');
        else
            dlmwrite(strcat(currentFolder, output, 'Mu_results\', filename), [H mu_real mu_im], '\t');
        end
    end
end
