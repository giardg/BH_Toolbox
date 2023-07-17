function [xvec, P_Joule, P_Hyst] = readLosses(freq, currentFolder, output)

global param_phy

mattype = param_phy.mattype;
longueur = param_phy.longueur; 
H0 = param_phy.H0;
Temp = param_phy.Temp;

switch mattype
    case 1
        pertesType = 'LossesLinear';
    case 2
        pertesType = 'LossesAnhyst';
    case 3
        pertesType = 'LossesHysteresis';
    case 4
        pertesType = 'LossesLimit';
    case 5
        pertesType = 'LossesEllipse';
    case 6
        pertesType = 'LossesEllipse2';
    case 7
        pertesType = 'LossesHysteresis2';
    otherwise
        error('Material not implemented\n')
end

filename = strcat('H0',{' '},string(H0),'_',pertesType,'_',string(longueur*1e3),'mm_',num2str(Temp),'deg_', num2str(freq/1e3), 'kHz.txt');

if isfile(strcat(currentFolder, output, '\Transient_results\', filename))
    A = load(strcat(currentFolder, output, '\Transient_results\', filename));
    xvec = A(:,1);
    P_Joule = A(:,2);
    P_Hyst = A(:,3);
else
     error('The file %s does not exist, it can be created by running FEM_transient_fortran',filename)
end