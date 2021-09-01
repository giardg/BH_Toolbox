function [xvec_Hyst, P_Hyst, xvec_Joule, P_Joules] = donnees_puissance(file_Hyst, file_Joule)

if isfile('Resultats_transitoires_COMSOL/'+string(file_Hyst))
    A = load(file_Hyst);
    xvec_Hyst = A(:,1);
    P_Hyst = A(:,2);
else
    error('Le fichier %s n''existe pas, vous pouvez le creer a partir du programme dans CalibrationProblem1D',file_Hyst)
end

if isfile('Resultats_transitoires_COMSOL/'+string(file_Joule))
    B = load(file_Joule);
    xvec_Joule = B(:,1);
    P_Joules = B(:,2);
else
    error('Le fichier %s n''existe pas, vous pouvez le creer a partir du programme dans CalibrationProblem1D',file_Joule)
end
