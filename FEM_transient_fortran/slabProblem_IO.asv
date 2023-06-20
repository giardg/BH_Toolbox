function [] = slabProblem_IO(freq, L, rho, H0_list, mattype, currentFolder, output)

% -------------------------------------------------------------------------
% MATLAB function that lauches the FEM executable compiled with Fortran for
% the resolution of the slab 1-D problem
% This program handles the Inputs/Outputs from the text files in the folder
% with the slab problem executable
%
% Input files
% Parameters.txt    : Contains the physical parameters
% Dirichelet.txt    : Contains the Dirichelet conditions (H(x=0,t))
%
% Output files
% Results.txt       : Contains H(x,t), B(x,t), J(x,t), Pj(x,t) and Ph(x,t)
%
% -------------------------------------------------------------------------

global param_phy

%Loop over H0
for H0 = H0_list
    
    
    %% FEM parameters that can be modified by a user
    
    % Path to ther fortran executable and I/O folder
    path =strcat(currentFolder,'\FEM_transient_fortran\SlabProblem\');
    
    % Number of elements
    ne          = 200;
    % Cycle duration (s)
    T           = 1/freq;
    % Number of time steps per period
    ntpc        = 2000;
    
    % Loop over the temperatures
    for k = 1:length(param_phy.Temp_list)
        
        % 1 : constant relative permeability
        %     B = mu0*mur*H
        if mattype == 1
            mur     = param_phy.mur_list(k);
            nomFich = ['LossesLinear_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
            % 2 : non-linear atan model
            %     B = mu0*H + 2*bsat/pi*atan(0.5*pi*mu0*(murmax-1)/bsat*H)
        elseif mattype == 2
            murmax  = param_phy.murmax_list(k);
            bsat    = param_phy.bsat_list(k);
            nomFich = ['LossesAnhyst_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
            % 3 : Preisach 3n coefficients non linear hysteretic model
            %     B+ = Sum ai*atan((H+ci)/bi)
            %     B- = Sum ai*atan((H-ci)/bi)
        elseif mattype == 3
            natan   = param_phy.natan_list(k);
            ai = param_phy.ai_list{k};
            ci = param_phy.ci_list{k};
            bi = param_phy.bi_list{k};
            nomFich = ['LossesHysteresis_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
            % 4 : High field limit case (not used very much)
            %     B = mu0*H+Bsat if H > Hsat
            %     B = mu0*H-Bsat if H < -Hsat
            %     3rd order polynomial between -Hsat and Hsat
        elseif mattype == 4
            Bsat = param_phy.Bsat_list(k);
            Hsat = param_phy.Hsat_list(k);
            nomFich = ['LossesLimit_' num2str(L*1E3) 'mm.txt'];
            
            % 5: Preisach Model defining F and G functions from an elliptic major cycle
            %    defined with a constant complex permeability value (not used very much)
        elseif mattype == 5
            mu_real = param_phy.mu_real_list(k);
            mu_imag = param_phy.mu_imag_list(k);
            nomFich = ['LossesEllipse_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
            % 6: Preisach Model with elliptic hysteretic curves (minor and major,
            %    not optimal...)
        elseif mattype == 6
            mu_real = param_phy.mu_real_list(k);
            mu_imag = param_phy.mu_imag_list(k);
            nomFich = ['LossesEllipse2_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
            % 7 : Preisach 4 parameters non linear hysteretic model
        elseif mattype == 7
            Br = param_phy.Br_list(k);
            Bsat = param_phy.Bsat_list(k);
            Hc = param_phy.Hc_list(k);
            s = param_phy.s_list(k);
            a = param_phy.a_list(k);
            b = param_phy.b_list(k);
            mur_max = param_phy.mur_max_list(k);
            nomFich = ['LossesHysteresis2_' num2str(L*1E3) 'mm_' num2str(param_phy.Temp_list(k)) 'deg_' num2str(freq/1e3) 'kHz.txt'];
            
        else
            error('Material type not implemented.')
        end
        
        nomFich = ['H0 ' num2str(H0) '_' nomFich];
        
        if ~isfile([currentFolder output 'Transient_results\' nomFich]) %check if data already exists
            
            % Applied magnetic field (Dirichelet condition) :
            % Amplitude 1 (A/m)
            %hamp1       = sqrt(1/(A-(B^2)/(4*C)));
            hamp1       = H0;
            
            % Frequency 1 (Hz)
            freq1       = 1/T;
            % Amplitude 2 (A/m) Only for SDF applications
            hamp2       = 0.0e+03;
            % Frequency 2 (Hz) Only for SDF applications
            freq2       = 190.0e+03;
            % Time vector (s)
            t = (0:ntpc-1)*T/ntpc;
            
            % Sin
            H_0_t = hamp1*sin(2*pi*freq1*t) + hamp2*sin(2*pi*freq2*t);
            
            %% Input writing :
            % Open file Parameters.txt
            fio = fopen(strcat(path,'Parameters.txt'),'wt');
            
            % Writing
            fprintf(fio,'Duree d''un cycle (s)\n');
            fprintf(fio,'%1.15e\n',T);
            fprintf(fio,'Nombre de pas de temps par cycle\n');
            fprintf(fio,'%u\n',ntpc);
            fprintf(fio,'Taille du domaine de calcul (m)\n');
            fprintf(fio,'%1.15e\n',L);
            fprintf(fio,'Nombre d''elements\n');
            fprintf(fio,'%u\n',ne);
            fprintf(fio,'Resistivite (Ohms/m)\n');
            fprintf(fio,'%1.15e\n',rho);
            fprintf(fio,'Type de courbe BH\n');
            fprintf(fio,'%u\n',mattype);
            fprintf(fio,'Parametres de la courbe BH\n');
            if mattype == 1
                fprintf(fio,'%1.15e\n',mur);
            elseif mattype == 2
                fprintf(fio,'%1.15e\n',murmax);
                fprintf(fio,'%1.15e\n',bsat);
            elseif mattype == 3
                fprintf(fio,'%1.15e\n',natan);
                for i = 1:natan
                    fprintf(fio,'%1.15e\n',ai(i));
                end
                for i = 1:natan
                    fprintf(fio,'%1.15e\n',bi(i));
                end
                for i = 1:natan
                    fprintf(fio,'%1.15e\n',ci(i));
                end
            elseif mattype == 4
                fprintf(fio,'%1.15e\n',Bsat);
                fprintf(fio,'%1.15e\n',Hsat);
            elseif mattype == 5
                fprintf(fio,'%1.15e\n',mu_real);
                fprintf(fio,'%1.15e\n',mu_imag);
                fprintf(fio,'%1.15e\n',hamp1);
            elseif mattype == 6
                fprintf(fio,'%1.15e\n',mu_real);
                fprintf(fio,'%1.15e\n',mu_imag);
                fprintf(fio,'%1.15e\n',hamp1);
            elseif mattype == 7
                fprintf(fio,'%1.15e\n',Br);
                fprintf(fio,'%1.15e\n',Bsat);
                fprintf(fio,'%1.15e\n',Hc);
                fprintf(fio,'%1.15e\n',s);
                fprintf(fio,'%1.15e\n',a);
                fprintf(fio,'%1.15e\n',b);
                fprintf(fio,'%1.15e\n',mur_max);
            else
                error('Material type not implemented.')
            end
            
            % Close file
            fclose(fio);
            
            % Open file Dirichelet.txt
            fio = fopen(strcat(path,'Dirichelet.txt'),'wt');
            
            % Writing
            for n = 1:ntpc
                fprintf(fio,'%1.15e\n',H_0_t(n));
            end
            
            % Close file
            fclose(fio);
            
            %% Launch time-transient simulation
            % system(char(strcat(path,{'SlabProblem '},path)));
            pathSys = [path '/SlabProblem/'];
            cd(pathSys)
            system('SlabProblem.exe '+string(path))
            % system('SlabProblem_tour.exe')
            
            %% Reading results
            % Open files
            fio = fopen(strcat(path,'Results.txt'),'rt');
            
            % Read buffer
            fgetl(fio); fgetl(fio); fgetl(fio);
            
            % Dimensions
            nt = fscanf(fio,'Nombre de pas de temps : %i\n',1);
            nx = fscanf(fio,'Nombre de points en x  : %i\n',1);
            
            % Vector initialization
            t  = zeros(1,nt);
            x  = zeros(nx,1);
            H  = zeros(nx,nt);
            B  = zeros(nx,nt);
            J  = zeros(nx,nt);
            Pj = zeros(nx,nt);
            Ph = zeros(nx,nt);
            
            % Read data
            for n = 1:nt
                
                % Read time step buffer
                fgetl(fio);
                t(n) = fscanf(fio,'Time = %e\n',1);
                fgetl(fio);
                
                % Read solution at time step n
                data = fscanf(fio,'%e',[6,nx]);
                fgetl(fio); fgetl(fio);
                
                % Copy data in structures
                x(1:nx)    = data(1,1:nx)';
                H(1:nx,n)  = data(2,1:nx)';
                B(1:nx,n)  = data(3,1:nx)';
                J(1:nx,n)  = data(4,1:nx)';
                Pj(1:nx,n) = data(5,1:nx)';
                Ph(1:nx,n) = data(6,1:nx)';
                
            end
            
            % Close files
            fclose(fio);
            
            Pj_x = mean(Pj.').';
            Ph_x = mean(Ph.').';
            
            dlmwrite([currentFolder output 'Transient_results\' nomFich], [x Pj_x Ph_x], '\t')
            
            delete([path 'Results.txt'])
            
        else
            disp(['Transient results for ' nomFich ' already exists, skipping to the next'])
        end
    end
end






