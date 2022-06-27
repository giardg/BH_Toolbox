function [xvec, P_Joule, P_Hyst] = readPuissance(freq, currentFolder, output)

global param_phy

mattype = param_phy.mattype;
longueur = param_phy.longueur; 
H0 = param_phy.H0;
Temp = param_phy.Temp;

switch mattype
    case 1
        pertesType = 'PertesLinear';
    case 2
        pertesType = 'PertesFoucault';
    case 3
        pertesType = 'PertesHysteresis';
    case 4
        pertesType = 'PertesLimites';
    case 5
        pertesType = 'PertesEllipse';
    case 6
        pertesType = 'PertesEllipse2';
    case 7
        pertesType = 'PertesHysteresis2';
    otherwise
        error('Materiau pas implemente\n')
end

filename = strcat('H0',{' '},string(H0),'_',pertesType,'_',string(longueur*1e3),'mm_',num2str(Temp),'deg_', num2str(freq/1e3), 'kHz.txt');

if isfile(strcat(currentFolder, output, '\Resultats_transitoire_fortran\', filename))
    A = load(strcat(currentFolder, output, '\Resultats_transitoire_fortran\', filename));
    xvec = A(:,1);
    P_Joule = A(:,2);
    P_Hyst = A(:,3);
else
     error('Le fichier %s n''existe pas, vous pouvez le creer a partir du programme dans FEM_transitoire_fortran',filename)
end