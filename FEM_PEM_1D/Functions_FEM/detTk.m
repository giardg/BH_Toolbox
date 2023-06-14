function jac = detTk(intervalle)
% Fonction qui retourne le jacobien (1D) pour un certain intervalle

jac = (intervalle(2)-intervalle(1))/2; %rapport de longueur de l'element sur l'element de reference