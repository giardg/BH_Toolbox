module ElementMod
  ! -------------------------------------------------------------------------------------------------------------------
  ! Mai 2018
  ! Maxime Tousignant
  ! maxime.tousignant@polymtl.ca
  ! -------------------------------------------------------------------------------------------------------------------
  ! Description :
  ! Quadrature de Gauss-Legendre de degree 4 pour integrer parfaitement des polynomes de Legendre d'ordre 2.
  ! Notes :
  ! La modification de l'ordre de la quadrature de Gauss ou des fonctions de bases va entrainer des erreurs.
  ! -------------------------------------------------------------------------------------------------------------------
  
  ! Quadrature de Gauss d'ordre 4 (hard coded)
  integer, parameter                                        :: nGauss     = 4
  double precision, dimension(nGauss), parameter                                                                        &
    :: xiGauss  = (/ -0.3399810435848563d0,  0.3399810435848563d0, -0.8611363115940526d0,  0.8611363115940526d0 /)
  double precision, dimension(nGauss), parameter                                                                        &
    :: wGauss   = (/  0.6521451548625461d0,  0.6521451548625461d0,  0.3478548451374538d0,  0.3478548451374538d0	/) 
    
  ! Fonctions de bases aux pts de gauss (hard coded)
  integer, parameter                                        :: degPsi     = 2
  double precision, dimension(nGauss)                       :: Psi1       = 0.5d0*xiGauss*(xiGauss-1.0d0)
  double precision, dimension(nGauss)                       :: Psi2       = 1.0d0-xiGauss*xiGauss
  double precision, dimension(nGauss)                       :: Psi3       = 0.5d0*xiGauss*(xiGauss+1.0d0)
  double precision, dimension(nGauss)                       :: dPsi1dxi   = xiGauss-0.5d0
  double precision, dimension(nGauss)                       :: dPsi2dxi   = -2.0d0*xiGauss
  double precision, dimension(nGauss)                       :: dPsi3dxi   = xiGauss+0.5d0
  
  ! Derivees des fonctions de base aux noeuds
  double precision, dimension(1+degPsi)                     :: dPsi1dxi_node  = (/-1.5d0 ,-0.5d0 , 0.5d0 /)
  double precision, dimension(1+degPsi)                     :: dPsi2dxi_node  = (/ 2.0d0 , 0.0d0 ,-2.0d0 /)
  double precision, dimension(1+degPsi)                     :: dPsi3dxi_node  = (/-0.5d0 , 0.5d0 , 1.5d0 /)
  
end module ElementMod
  