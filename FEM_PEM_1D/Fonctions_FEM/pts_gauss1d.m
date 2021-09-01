function [pts,poids] = pts_gauss1d(npts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pts_gauss1d
% 
% Donne les points d'evaluation et les poids pour la quadrature de gauss en
% 1D avec npts points (issu du cours MTH8207 A2019)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (npts == 1)    % degre de precision = 1

   pts   = [ 0 ];

   poids = [ 2 ];

elseif (npts == 2) % degre de precision = 3

   pts   = [ - 1/sqrt(3) ;
               1/sqrt(3) ];

   poids = [ 1 ;
             1 ];

elseif (npts == 3) % degre de precision = 5


   pts   = [ -0.774596669241483 ;
              0.0;
              0.774596669241483 ];

   poids = [  0.5555555555555556 ;
              0.8888888888888889 ;
              0.5555555555555556  ];

elseif (npts == 4) % degre de precision = 7


   pts   = [ -0.861136311594052 ;
             -0.339981043584856 ;
              0.339981043584856 ;
              0.861136311594052];

   poids = [ 0.347854845137454 ;
             0.652145154862545 ;
             0.652145154862545 ;
             0.347854845137454 ];

else
   printf('Erreur dans pts_gauss1d(): quadrature avec %d points pas encore implementee.\n',npts);
end
