clc;
clearvars
close all
fclose all;

x0=400;
y0=400;
width=700;
height=550;

%% Parametres
mu0 = 4*pi*1E-7;
Tempvec = [24.7297,113.3906,246.0131,409.2073,461.1946,506.9757,522.4292,547.1615,560.8959,584.9624,607.1725,620.1256,632.5545,645.9273,663.0356, 775]; %liste des fichiers disponible
%Tempvec = [25];
H0_list = [1:4,5:5:100]*1e3; %liste des fichiers disponible
freq = 10*1e3;
longueur = 2.5e-3;
correction = true;
mattype = 7;

H0_grid = [];
H_grid = [];
mu_real_grid = [];
mu_im_grid = [];

H0_export = [];
H_export = [];
mu_real_export = [];
mu_im_export = [];

if mattype == 1
    typename = 'Linear';
elseif mattype == 2
    typename = 'Foucault';
elseif mattype == 3
    typename = 'Hysteresis';
elseif mattype == 4
    typename = 'Limites';
elseif mattype == 7
    typename = 'Hysteresis2';
else
    error('Le type de materiau choisi n''est pas implemente.')
end


mat_real = [];
mat_imag = [];

%% Ouverture de tous les fichiers
for H0 = H0_list
    for Temp = Tempvec
        
        if Temp ~= Tempvec(end) && ~any(Temp == Tempvec(end))
            filename = strcat('H0',{' '},string(H0),'_mu',typename,'_',num2str(Temp),'deg_',string(freq/1e3),'kHz');
            if correction filename = strcat(filename,'_corrected'); end
            filename = strcat(filename,'.txt');
            
            if isfile(filename)
                
                fid = fopen(filename,'r');
                data = textscan(fid, '%f%f%f', 'CollectOutput',1);
                data = cell2mat(data);
                H = real(data(:,1));
                mu_real = (data(:,2));
                mu_im = (data(:,3));
                
                new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real];
                new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im];
                mat_real = [mat_real; new_data_real];
                mat_imag = [mat_imag; new_data_imag];
                
                
            else
                error('Le fichier %s n''existe pas',filename)
            end
            
            %     if rem(H0*1e-3,5) ~= 0
            %         cpt = cpt+1;
            %         continue
            %     end
            
%         elseif Temp == Tempvec(end-6)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0.8+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0.8];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
%             
%         elseif Temp == Tempvec(end-5)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0.7+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0.7];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
%             
%         elseif Temp == Tempvec(end-4)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0.6+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0.6];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
%             
%             
%         elseif Temp == Tempvec(end-3)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0.5+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0.5];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
%             
%         elseif Temp == Tempvec(end-2)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0.4+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0.4];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
%             
%         elseif Temp == Tempvec(end-1)
%             new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0+1];
%             new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0];
%             mat_real = [mat_real; new_data_real];
%             mat_imag = [mat_imag; new_data_imag];
        else
            new_data_real = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_real.*0+1];
            new_data_imag = [H, H0*ones(size(data,1),1),Temp*ones(size(data,1),1),mu_im.*0];
            mat_real = [mat_real; new_data_real];
            mat_imag = [mat_imag; new_data_imag];
        end
        
    end
end

[~,idx] = sort(mat_real(:,1).^2+mat_real(:,2).^2+mat_real(:,3).^2);
mat_real = mat_real(idx,:);
[~,idx] = sort(mat_imag(:,1).^2+mat_imag(:,2).^2+mat_imag(:,3).^2);
mat_imag = mat_imag(idx,:);
mat_real = unique(mat_real,'rows');
mat_imag = unique(mat_imag,'rows');

%% Figure

dlmwrite(strcat('table_mu_real_temp_test.txt'), mat_real, '\t');
dlmwrite(strcat('table_mu_imag_temp_test.txt'), mat_imag, '\t');