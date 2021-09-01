function [xvec, P_Joule, P_Hyst] = read_puissance()

global param_phy

mattype = param_phy.mattype;
longueur = param_phy.longueur; 
H0 = param_phy.H0;

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
    otherwise
        error('Materiau pas implemente\n')
end

filename = strcat('H0',{' '},string(H0),'_',pertesType,'_',string(longueur*1e3),'mm.txt');

if isfile('Resultats_transitoires_fortran/'+filename)
    A = load(filename);
    xvec = A(:,1);
    P_Joule = A(:,2);
    P_Hyst = A(:,3);
else
     error('Le fichier %s n''existe pas, vous pouvez le creer a partir du programme dans FEM_transitoire_fortran',filename)
end