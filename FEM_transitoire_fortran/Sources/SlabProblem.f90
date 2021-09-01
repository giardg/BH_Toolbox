module SlabProblemMod
  ! -------------------------------------------------------------------------------------------------------------------
  ! Mai 2018
  ! Maxime Tousignant
  ! maxime.tousignant@polymtl.ca
  ! -------------------------------------------------------------------------------------------------------------------
  ! Description :
  ! Module de calcul du slab problem. Ce programme permet la resolution transitoire de l'equation 1D non lineaire de
  ! diffusion des courants de Foucault avec un formulation en H.
  !
  ! Equation resolue :
  ! rho * d^2H/dx^2 = - dB(H)/dt
  !        H(x=0,t) = H0(t)
  !    dH(x=L,t)/dx = 0
  !  
  ! -------------------------------------------------------------------------------------------------------------------
  
  ! Modules
  use ElementMod
  use ErrorMod
  use MaterialMod
  
  ! Parametres et variables
  implicit none
  
  ! Constantes
  ! Pi
  double precision, parameter, private                      :: pi         = 3.1415926535897932d+00
  ! Nombre d'or
  double precision, parameter, private                      :: golden     = 1.6180339887498949d+00

  ! Parametres des algorithmes
  ! Tolerance sur l'atteinte du regime permanent
  integer, parameter, private                               :: nCmax      = 50
  integer, parameter, private                               :: nCramp     = 2
  double precision, parameter, private                      :: tolSteady  = 1.0d-03
  ! Tolerance Newton-Raphson
  integer, parameter, private                               :: iterNRmax  = 100
  double precision, parameter, private                      :: tolNewton  = 1.0d-12
  ! Tolerance algorithme de recherche du facteur de relaxation
  integer, parameter, private                               :: iterRelMax = 8
  
  ! Valeur caracteristique du residu
  double precision, private                                 :: Rcara
  
  ! Format d'ecriture
  character(len=8)                                          :: f1dble = "(E23.15)"
  character(len=43)                                         :: f6dble = "(E23.15,E23.15,E23.15,E23.15,E23.15,E23.15)"
  
  ! Routine et fonctions privees
  ! Transitoire
  private                                                   :: GetPastTimes
  private                                                   :: PredictInitialSolution
  private                                                   :: EvalSteadyState
  ! Non lineaire
  private                                                   :: NonLinSimulation
  private                                                   :: UpdateSolution
  ! Systeme lineaire
  private                                                   :: Assemblage
  private                                                   :: Solve
  ! Post-pro
  private                                                   :: CalculCurrentDensity
  private                                                   :: CalculLosses
  ! General
  private                                                   :: Norm2
  
  contains
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine ReadParameters (path , ne,ndof,ntpc,dx,dt,rho,bh,err)
    ! Description :
    ! Lecture des parametres du probleme dans HOME/SlabProblem/Parameters.txt.
    implicit none
    ! Inputs :
    ! Path jusqu'au repertoire I/O
    character(len=maxStrLen), intent(in)                    :: path
    ! Outputs :
    ! Discretisation temps et espace
    integer, intent(out)                                    :: ne, ndof, ntpc 
    double precision, intent(out)                           :: dx, dt
    ! Proprietes du materiau                                
    double precision, intent(out)                           :: rho
    double precision, dimension(maxParam), intent(out)      :: bh
    ! Indicateur d'erreur                                 
    integer, intent(out)                                    :: err
    
    ! Variables
    integer                                                 :: i
    ! Fichiers
    character(len=maxStrLen)                                :: filename
    ! Descriptifs
    character(len=maxStrLen)                                :: junk
    ! Duree d'un cycle
    double precision                                        :: period
    ! Taille du domaine
    double precision                                        :: L
    ! Type de courbe BH
    integer                                                 :: mattype
    ! Erreur                                                
    character(len=maxStrLen)                                :: from
    character(len=maxStrLen)                                :: message
    logical :: file_exists
    ! Initialisation
    err = 0
    from = "ReadParameters"
    message = " "
    
    ! Nom du fichier
    filename = trim(path) // "Parameters.txt"
    
    ! Ouverture du fichier contenant les resultats pour une periode
    open (unit=1, file=trim(filename), iostat=err,  status="old", action="read")
    if ( err /= 0 ) then
      message = "Error opening file Parameters.txt"
      call Error (from,message )
      return
    end if
    
    ! Lecture du fichier
    ! Duree d'un cycle (s)
    read(1,*) junk
    read(1,*) period
    
    ! Nombre de pas de temps par cycle
    read(1,*) junk
    read(1,*) ntpc
    
    ! Pas de temps
    dt = period/dble(ntpc)
    
    ! Taille du domaine de calcul (m)
    read(1,*) junk
    read(1,*) L
    
    ! Nombre d'elements
    read(1,*) junk
    read(1,*) ne
    
    ! Calcul du nombre de degres de liberte (elements du second ordre)
    ndof = degPsi*ne
    
    ! Calcul du dx
    dx = L/dble(ne)
    
    ! Resistivite (Ohms/m)
    read(1,*) junk
    read(1,*) rho
    
    ! Type de courbe BH
    read(1,*) junk
    read(1,*) mattype
    bh(1) = dble(mattype)
    
    ! Parametres de la courbe BH
    read(1,*) junk
    ! 1 : lineaire permeabilite relative constante
    if ( mattype == 1 ) then
      read(1,*) bh(2)
    ! 2 : non lineaire modele arctangente
    else if ( mattype == 2 ) then
      read(1,*) bh(2)
      read(1,*) bh(3)
    ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
    else if ( mattype == 3 ) then
      read(1,*) bh(2)
      do i = 1, 3*nint(bh(2))
        read(1,*) bh(2+i)
      end do
    else if (mattype == 4) then
      read(1,*) bh(2)
      read(1,*) bh(3)
    else
      err = 1
      message = "Le type de materiau choisi n''est pas implemente."
      call Error (from,message )
    end if
      
    ! Fermeture du fichier
    close (1)  
    
    return
  end subroutine ReadParameters
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine ReadDirichelet (path,ntpc, H_0_t,err)
    ! Description :
    ! Lecture de la condition de Dirichelet a partir du fichier texte HOME/SlabProblem/Dirichelet.txt.
    implicit none
    ! Inputs :
    ! Path jusqu'au repertoire I/O
    character(len=maxStrLen), intent(in)                    :: path
    ! Discretisation
    integer, intent(in)                                     :: ntpc
    ! Outputs :
    ! Condition de Dirichelet
    double precision, dimension(ntpc), intent(out)          :: H_0_t
    ! Indicateur d'erreur
    integer, intent(out)                                    :: err
    
    ! Variables
    integer                                                 :: n
    ! Fichiers
    character(len=maxStrLen)                                :: filename
    ! Erreur                                                
    character(len=maxStrLen)                                :: from
    character(len=maxStrLen)                                :: message
    
    ! Initialisation
    err = 0
    from = "ReadDirichelet"
    message = " "
    
    ! Nom du fichier
    filename = trim(path) // "Dirichelet.txt"
    
    ! Ouverture du fichier contenant les resultats pour une periode
    open (unit=1, file=trim(filename), iostat=err,  status="old", action="read")
    if ( err /= 0 ) then
      message = "Error opening file Dirichelet.txt"
      call Error (from,message )
      return
    end if
    
    ! Lecture du fichier
    do n = 1, ntpc
      read(1,*) H_0_t(n)
    end do
    
    ! Fermeture du fichier
    close(1)
    
    return
  end subroutine ReadDirichelet
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine InitializeSlabProblem (path,ne,ndof,ntpc,dx,dt,rho, bh , H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t,err)
    ! Description :
    ! Initialisation des structures de donnees et calculs pre-resolution.
    implicit none
    ! Inputs :
    ! Path jusqu'au repertoire I/O
    character(len=maxStrLen), intent(in)                    :: path
    ! Discretisation temps et espace
    integer, intent(in)                                     :: ne, ndof, ntpc 
    double precision, intent(in)                            :: dx, dt
    ! Proprietes du materiau                               
    double precision, intent(in)                            :: rho
    ! Inouts :
    ! Proprietes du materiau                               
    double precision, dimension(maxParam), intent(inout)    :: bh
    ! Outputs :
    ! Champs temps et espace                                
    double precision, dimension(:,:), allocatable, intent(out)  :: H_x_t, B_x_t, J_x_t
    ! Pertes
    double precision, dimension(:,:), allocatable, intent(out)  :: Pj_x_t, Ph_x_t
    ! Indicateur d'erreur                                 
    integer, intent(out)                                    :: err
    
    ! Variables
    ! Erreur                                                
    character(len=maxStrLen)                                :: from
    character(len=maxStrLen)                                :: message
    
    ! Initialisation
    err = 0
    from = "InitializeSlabProblem"
    message = " "
    
    ! Allocation des structures de donnees
    allocate(H_x_t(ndof+1,ntpc) , stat=err)
    if ( err == 0 ) allocate(B_x_t(ndof+1,ntpc) , stat=err)
    if ( err == 0 ) allocate(J_x_t(ndof+1,ntpc) , stat=err)
    if ( err == 0 ) allocate(Pj_x_t(ndof+1,ntpc) , stat=err)
    if ( err == 0 ) allocate(Ph_x_t(ndof+1,ntpc) , stat=err)
    if ( err /= 0 ) then
      err = 1
      message = "Memoire insuffisante."
      call Error (from,message )
    end if
    
    ! Initialization a zero
    H_x_t = 0.0d0
    B_x_t = 0.0d0
    J_x_t = 0.0d0
    Pj_x_t = 0.0d0
    Ph_x_t = 0.0d0
    
    ! Recuperation des conditions de Dirichelet
    call ReadDirichelet (path,ntpc , H_x_t(ndof+1,1:ntpc),err)
    if ( err /= 0 ) return
    
    ! Valeur carateristique du residu
    Rcara = rho/dx * 0.5d0*( maxval(H_x_t(ndof+1,1:ntpc)) - minval(H_x_t(ndof+1,1:ntpc)) )
    if ( err /= 0 ) then
      err = 1
      message = "La valeur caracteristique du residu doit etre >0."
      call Error (from,message )
    end if
    
    ! Initialisation du materiau
    call InitializeMaterial (ne , bh , err)
    if ( err /= 0 ) return
    
    return
  end subroutine InitializeSlabProblem
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine FinalizeSlabProblem ( H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t )
    ! Description :
    ! Libere la memoire allouee.
    implicit none
    ! Inouts :
    ! Donnees physiques
    double precision, dimension(:,:), allocatable, intent(inout)  :: H_x_t, B_x_t, J_x_t, Pj_x_t, Ph_x_t
    
    if ( allocated(H_x_t) ) deallocate(H_x_t)
    if ( allocated(B_x_t) ) deallocate(B_x_t)
    if ( allocated(J_x_t) ) deallocate(J_x_t)
    if ( allocated(Pj_x_t) ) deallocate(Pj_x_t)
    if ( allocated(Ph_x_t) ) deallocate(Ph_x_t)
    
    call FinalizeMaterial ()
    
    return
  end subroutine FinalizeSlabProblem
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculSlabProblem (ne,ndof,ntpc,dx,dt,rho,bh , H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t , err)
    ! Description :
    ! Boucle transitoire. La routine calcule une nouvelle periode jusqu'a ce qu'un critere definissant le regime
    ! permanent soit atteint.
    implicit none
    ! Inputs :
    ! Discretisation temps et espace
    integer, intent(in)                                     :: ne, ndof, ntpc 
    double precision, intent(in)                            :: dx, dt
    ! Proprietes du materiau                                
    double precision, intent(in)                            :: rho
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Inpouts :                               
    ! Champs temps et espace                                
    double precision, dimension(ndof+1,ntpc), intent(inout) :: H_x_t, B_x_t, J_x_t
    ! Pertes
    double precision, dimension(ndof+1,ntpc), intent(inout) :: Pj_x_t, Ph_x_t
    ! Outputs :
    ! Indicateur d'erreur                                 
    integer, intent(out)                                    :: err
    
    
    ! Variables :
    ! Iterations
    integer                                                 :: nC, n
    ! Erreur
    character(len=maxStrLen)                                :: from, message
    ! Indices des pas precedents                            
    integer                                                 :: n_1, n_2, n_3
    ! Atteinte regime permanent                           
    double precision, dimension(ndof+1,ntpc)                :: B_x_t_
    double precision                                        :: steadyState
    ! Condition de Dirichelet
    double precision, dimension(ntpc)                       :: H_0_t
    
    ! Initialisation
    err = 0
    from = "TransientSimulation"
    message = " "
    steadyState = 1.0d0
    nC = 0
    
    ! Copie de la condition de Dirichelet
    H_0_t = H_x_t(ndof+1,1:ntpc)
     
    ! Calcul du regime permanent
    do while ( steadyState >= tolSteady .and. nC < nCmax )
      
      ! Nouveau cycle
      nC = nC + 1
      
      ! Affichage
      write(*,*)
      write(*,*) "Cycle no.   ",nC
      write(*,*) "------------------------------------------------------"
      
      ! Copie de l'induction au cycle precedent
      B_x_t_ = B_x_t
      
      ! Calcul d'une periode
      do n = 1, ntpc
        
        ! Affichage
        write(*,*)
        write(*,*) "Pas no.     ", n
        
        ! Indices des pas precedents
        call GetPastTimes (ntpc,n , n_1,n_2,n_3)
        
        ! Ramp-up de la condition de Dirichelet pour la premiere periode
        ! La condition de Dirichelet est sockee a la fin du vecteur solution
        if ( nC <= nCramp ) then
          H_x_t(ndof+1,n) = 0.5d0*( 1.0d0 - cos(pi*dble((nC-1)*ntpc+n-1)/dble(nCramp*ntpc)) )*H_0_t(n)
        else
          H_x_t(ndof+1,n) = H_0_t(n)
        end if
        
        ! Prediction de la solution
        call PredictInitialSolution (ndof,H_x_t(1:ndof,n_1),H_x_t(1:ndof,n_2),H_x_t(1:ndof,n_3) , H_x_t(1:ndof,n))

        ! Resolution non lineaire
        call NonLinSimulation (ne,ndof,dx,dt,rho,bh,B_x_t(1:ndof+1,n_1),B_x_t(1:ndof+1,n_2) , H_x_t(1:ndof+1,n) , err)
        if ( err /= 0 ) return
        
        ! Calcul de l'induction apres resolution et gestion de l'historique (si hysteresis)
        call UpdateMaterialBH (ne,ndof,bh,H_x_t(1:ndof+1,n) , B_x_t(1:ndof+1,n))
      
      enddo
      
      ! Evalution de l'atteinte du regime permanent
      call EvalSteadyState (ndof,ntpc,B_x_t_(1:ndof,1:ntpc),B_x_t(1:ndof,1:ntpc) , steadyState)
      
      ! Affichage
      write(*,*) 
      write(*,*) "SteadyState               = ",steadyState
      
    end do
    
    if ( steadyState >= tolSteady) then
      message = "Regime permanent non atteint."
      call Warning (from,message )
    end if
    
    ! Post-processing
    call CalculCurrentDensity (ne,ndof,ntpc,dx,H_x_t , J_x_t)
    
    ! Calcul des pertes
    call CalculLosses (ndof,ntpc,dx,dt,rho,H_x_t,B_x_t,J_x_t , Pj_x_t,Ph_x_t)
    
    return
  end subroutine CalculSlabProblem
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine GetPastTimes (ntpc,n , n_1,n_2,n_3)
    ! Description :
    ! Recupere les indices des pas de temps precedents.
    implicit none
    ! Inputs :
    ! Nombre de pas de temps par periode
    integer, intent(in)                                     :: ntpc
    ! Indice du pas de temps courant
    integer, intent(in)                                     :: n
    ! Outputs :
    ! Indices des pas de temps precedents
    integer, intent(out)                                    :: n_1, n_2, n_3
    
    ! Modulo ntpc
    if ( n == 1 ) then
      n_1 = ntpc
      n_2 = ntpc-1
      n_3 = ntpc-2
      
    else if ( n == 2 ) then
      n_1 = 1
      n_2 = ntpc
      n_3 = ntpc-1
      
    else if ( n == 3 ) then
      n_1 = 2
      n_2 = 1
      n_3 = ntpc
      
    else
      n_1 = n-1
      n_2 = n-2
      n_3 = n-3
      
    end if
    
    return
  end subroutine GetPastTimes
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine PredictInitialSolution (ndof,H_1,H_2,H_3 , H_n)
    ! Description :
    ! Extrapolation d'ordre 2 pour initialiser la solution.
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ndof
    ! Champ aux pas precedents                              
    double precision, dimension(ndof), intent(in)           :: H_1, H_2, H_3
    ! Outputs :                                           
    ! Prediction                                            
    double precision, dimension(ndof), intent(out)          :: H_n
    
    ! Prediction de la solution au pas n
    H_n = 0.5d0*( 5.0d0*H_1 - 4.0d0*H_2 + H_3 )
    
    return
  end subroutine PredictInitialSolution
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine EvalSteadyState (ndof,ntpc,U_x_t_,U_x_t , steadyState)
    ! Description :
    ! Evalutation de la variation entre la solution pour la nouvelle periode et l'ancienne.
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ndof, ntpc
    ! Solutions
    double precision, dimension(ndof,ntpc), intent(in)      :: U_x_t_, U_x_t
    ! Output :
    double precision, intent(out)                           :: steadyState
    
    ! Variables :
    integer                                                 :: n
    double precision                                        :: u, du
    
    ! Initialisation
    steadyState = 0.0d0
    
    ! Boucle sur les pas de temps
    do n = 1, ntpc
      
      ! Norme de la solution
      u = Norm2(ndof,U_x_t(1:ndof,n))
      
      ! Norme de la difference
      du = Norm2(ndof,U_x_t_(1:ndof,n)-U_x_t(1:ndof,n))
      
      if ( u > epsilon(1.0d0) ) then
        steadyState = steadyState + du/u
      endif
      
    enddo
    
    steadyState = steadyState/dble(ntpc)
      
    return
  end subroutine EvalSteadyState
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine NonLinSimulation (ne,ndof,dx,dt,rho,bh,B_1,B_2 , H , err)
    ! Description :
    ! Algorithme de Newton-Raphson pour resoudre le systeme non lineaire.
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ne, ndof
    double precision, intent(in)                            :: dx, dt
    ! Proprietes du materiau                                
    double precision, intent(in)                            :: rho
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Indution aux pas precedents                         
    double precision, dimension(ndof+1), intent(inout)      :: B_1, B_2
    ! InOut :                                             
    ! Prediction et condition de Dirichelet / Solution      
    double precision, dimension(ndof+1), intent(inout)      :: H
    ! Output :                                              
    ! Indicateur d'erreur                                 
    integer, intent(out)                                    :: err
                                                            
    ! Variables :                                           
    ! Erreur                                                
    character(len=maxStrLen)                                :: from, message
    ! Matrice Jacobienne                                    
    double precision, dimension(ndof,ndof)                  :: M
    ! Residu (avec signe -)                                              
    double precision, dimension(ndof)                       :: R
    ! Solution                                              
    double precision, dimension(ndof)                       :: dH
    ! Newton-Raphson                                        
    integer                                                 :: iterNR
    double precision                                        :: critNewton
    ! Facteur de relaxation                               
    double precision                                        :: alpha
    
    ! Initialisation
    err = 0
    from = "NonLinSimulation"
    message = " "
    iterNR = 1
    alpha  = 1.0d0
    
    ! Assemblage initial
    call Assemblage (ne,ndof,dx,dt,rho,bh,B_1,B_2,H , M,R)
    
    ! Evaluation de la convergence
    critNewton = Norm2(ndof,R)/Rcara
    
    ! Affichage
    write(*,*) "Residu initial            : ",critNewton
    
    ! Solve
    call Solve (ndof , M,R , dH,err)
    if ( err /= 0 ) return
    
    ! Mise a jour de la solution (alorithme de relaxation)
    call UpdateSolution (ne,ndof,dx,dt,rho,bh,B_1,B_2,critNewton,dH, alpha,H , M,R)
    
    ! Evaluation de la convergence
    critNewton = Norm2(ndof,R)/Rcara
    
    ! Affichage
    write(*,*) "Residu iter ",iterNR," : ",critNewton
    
    ! Boucle de Newton-Raphson
    do while ( critNewton >= tolNewton .and. iterNR < iterNRmax )
      
      ! Icrementation
      iterNR = iterNR + 1
      
      ! Solve
      call Solve (ndof , M,R , dH,err)
      if ( err /= 0 ) return
      
      ! Mise a jour de la solution (alorithme de relaxation)
      call UpdateSolution (ne,ndof,dx,dt,rho,bh,B_1,B_2,critNewton,dH, alpha,H , M,R)
      
       ! Affichage
       if ( alpha < 1.0d0 ) write(*,*) "  Relaxation              : ",alpha
      
      ! Evaluation de la convergence
      critNewton = Norm2(ndof,R)/Rcara
      
      ! Affichage
      write(*,*) "Residu iter ",iterNR," : ",critNewton
      
    end do
    
    ! Verification du critere de convergence
    if (critNewton >= tolNewton ) then
      err = 1
      message = "Critere de convergence de Newton-Raphson non atteint."
      call Error (from,message )
    
    end if
      
    return
  end subroutine NonLinSimulation
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine UpdateSolution (ne,ndof,dx,dt,rho,bh,B_1,B_2,crit0,dH , alpha,H , M,R)
    ! Description :
    ! Calcul la nouvelle solution H_k+1 = H_k + alpha * dH_k
    ! Le facteur de relaxation alpha est choisi pour etre le plus grand possible tout en faisant diminuer le residu.
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ne, ndof
    double precision, intent(in)                            :: dx, dt
    ! Proprietes du materiau                                
    double precision, intent(in)                            :: rho
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Indution aux pas precedents                         
    double precision, dimension(ndof+1), intent(inout)      :: B_1, B_2
    ! Critere de convergence                                
    double precision, intent(in)                            :: crit0
    ! Correction a la solution                              
    double precision, dimension(ndof), intent(in)           :: dH
    ! Inout :                                             
    ! Facteur de relaxation                               
    double precision, intent(inout)                         :: alpha
    ! Solution
    double precision, dimension(ndof+1), intent(inout)      :: H
    ! Outputs :
    ! Systeme lineaire
    double precision, dimension(ndof,ndof), intent(out)     :: M
    double precision, dimension(ndof), intent(out)          :: R
    
    ! Variables :
    double precision                                        :: crit1
    double precision, dimension(ndof+1)                     :: H0

    ! Parametres de l'algoritme de relaxation
    double precision                                        :: alphaMin = 1.0d0/golden**(dble(iterRelMax)-0.9d0)
    double precision                                        :: alphaOut = 1.0d0/golden**4.0d0 
    
    ! Copie de la solution initiale
    H0 = H
    
    ! On tente d'abord d'augmenter le facteur de relaxation
    if ( alpha >= 1.0d0/golden ) then
      alpha = 1.0d0
    else
      alpha = alpha*golden
    end if
    
    ! Mise-a-jour de la solution
    H(1:ndof) = H0(1:ndof) + alpha*dH
    
    ! Evaluation du systeme
    call Assemblage (ne,ndof,dx,dt,rho,bh,B_1,B_2,H , M,R)
    
    ! Evaluation de la convergence
    crit1 = Norm2(ndof,R)/Rcara
    
    ! Recherche du facteur de relaxation maximal qui fait diminuer le residu
    do while ( crit1 >= crit0 .and. alpha > alphaMin  )
      
      ! Reduction du facteur de relaxation
      alpha = alpha/golden
      
      ! Mise-a-jour de la solution
      H(1:ndof) = H0(1:ndof) + alpha*dH
      
      ! Evaluation du systeme
      call Assemblage (ne,ndof,dx,dt,rho,bh,B_1,B_2,H , M,R)
    
      ! Evaluation de la convergence
      crit1 = Norm2(ndof,R)/Rcara
      
    end do
    
    ! On tente de sortir du minimum local si c'est le cas
    if ( crit1 >= crit0 ) then
      
      ! Facteur de relaxation moyen
      alpha = alphaOut
      
      ! Mise-a-jour de la solution
      H(1:ndof) = H0(1:ndof) + alpha*dH
      
      ! Evaluation du systeme
      call Assemblage (ne,ndof,dx,dt,rho,bh,B_1,B_2,H , M,R)
      
    end if
      
    return
  end subroutine UpdateSolution
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine Assemblage (ne,ndof,dx,dt,rho,bh,B_1,B_2,H , M,R)
    ! Description :
    ! Intergration et assemblage du systeme d'equation pour le slab problem.
    implicit none
    ! Inputs :
    ! Discretisation                                        
    integer, intent(in)                                     :: ne, ndof
    double precision, intent(in)                            :: dx, dt
    ! Proprietes du materiau                                
    double precision, intent(in)                            :: rho
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champs                                                
    double precision, dimension(ndof+1), intent(in)         :: B_1, B_2
    double precision, dimension(ndof+1), intent(in)         :: H
    ! Outputs :
    ! Matrice Jacobienne
    double precision, dimension(ndof,ndof), intent(out)     :: M
    ! Residu (avec signe -)   
    double precision, dimension(ndof), intent(out)          :: R
    
    ! Variables :
    ! Iterations
    integer                                                 :: e, i, k
    ! Numero du point de calcul
    integer                                                 :: noPt
    ! Valeurs physiques aux pts de calcul
    double precision                                        :: Hk, dHkdx, Bk, Bk_1, bk_2, dBdHk
    
    ! Initialisation
    M = 0.0d0
    R = 0.0d0
    
    ! Traitement du premier element avec la condition de Dirichelet
    ! Les degrees de libertes connus sont stocke a la fin du vecteur
    i = 0
    
    ! Evaluation des quantites physiques aux pts de calcul
    do k = 1, nGauss
      
      ! Numero identifiant le pt de calcul
      noPt = GetNoPt(1,k,2 )
      
      ! Champ mag et sa derivee
      Hk = H(ndof+1)*Psi1(k)+H(i+1)*Psi2(k)+H(i+2)*Psi3(k)
      dHkdx = 2.0d0/dx*(H(ndof+1)*dPsi1dxi(k)+H(i+1)*dPsi2dxi(k)+H(i+2)*dPsi3dxi(k))
      
      ! Induction aux pas precedents
      Bk_1 = B_1(ndof+1)*Psi1(k)+B_1(i+1)*Psi2(k)+B_1(i+2)*Psi3(k)
      Bk_2 = B_2(ndof+1)*Psi1(k)+B_2(i+1)*Psi2(k)+B_2(i+2)*Psi3(k)
      
      ! Induction au pas present et sa derivee
      call CalculLocalBH (noPt,bh,Hk , Bk,dBdHk)
    
      ! Integration
      ! Construction de la matrice elementaire
      M(i+1,i+1) = M(i+1,i+1) + wGauss(k)*( 2.0d0/dx*rho*dPsi2dxi(k)*dPsi2dxi(k) +                                      &
                                            0.75d0*dx/dt*dBdHk*Psi2(k)*Psi2(k) )
      M(i+1,i+2) = M(i+1,i+2) + wGauss(k)*( 2.0d0/dx*rho*dPsi2dxi(k)*dPsi3dxi(k) +                                      &
                                            0.75d0*dx/dt*dBdHk*Psi2(k)*Psi3(k) )
      M(i+2,i+2) = M(i+2,i+2) + wGauss(k)*( 2.0d0/dx*rho*dPsi3dxi(k)*dPsi3dxi(k) +                                      &
                                            0.75d0*dx/dt*dBdHk*Psi3(k)*Psi3(k) )
    
      ! Construction du residu
      R(i+1) = R(i+1) - wGauss(k)*( rho*dHkdx*dPsi2dxi(k) + 0.25d0*dx/dt*(3.0d0*Bk-4.0d0*Bk_1+Bk_2)*Psi2(k) )
      R(i+2) = R(i+2) - wGauss(k)*( rho*dHkdx*dPsi3dxi(k) + 0.25d0*dx/dt*(3.0d0*Bk-4.0d0*Bk_1+Bk_2)*Psi3(k) )

    end do
    
    ! Symetrie
    M(i+2,i+1) = M(i+1,i+2)
    
    ! Boucle sur les autres elements
    do e = 2, ne
      
      ! Pointeur pour les degrees de liberte
      i = i + degPsi
      
      ! Evaluation des quantites physiques aux pts de calcul
      do k = 1, nGauss
        
        ! Numero identifiant le pt de calcul
        noPt = GetNoPt(e,k,2 )
      
        ! Champ mag et sa derivee
        Hk = H(i)*Psi1(k)+H(i+1)*Psi2(k)+H(i+2)*Psi3(k)
        dHkdx = 2.0d0/dx*(H(i)*dPsi1dxi(k)+H(i+1)*dPsi2dxi(k)+H(i+2)*dPsi3dxi(k))
      
        ! Induction aux pas precedents
        Bk_1 = B_1(i)*Psi1(k)+B_1(i+1)*Psi2(k)+B_1(i+2)*Psi3(k)
        Bk_2 = B_2(i)*Psi1(k)+B_2(i+1)*Psi2(k)+B_2(i+2)*Psi3(k)
      
        ! Induction au pas present et sa derivee
        call CalculLocalBH (noPt,bh,Hk , Bk,dBdHk)
        
        ! Integration
        ! Construction de la matrice elementaire
        M(i,i)     = M(i,i)     + wGauss(k)*( 2.0d0/dx*rho*dPsi1dxi(k)*dPsi1dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi1(k)*Psi1(k) )
        M(i,i+1)   = M(i,i+1)   + wGauss(k)*( 2.0d0/dx*rho*dPsi1dxi(k)*dPsi2dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi1(k)*Psi2(k) )
        M(i,i+2)   = M(i,i+2)   + wGauss(k)*( 2.0d0/dx*rho*dPsi1dxi(k)*dPsi3dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi1(k)*Psi3(k) )
        M(i+1,i+1) = M(i+1,i+1) + wGauss(k)*( 2.0d0/dx*rho*dPsi2dxi(k)*dPsi2dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi2(k)*Psi2(k) )
        M(i+1,i+2) = M(i+1,i+2) + wGauss(k)*( 2.0d0/dx*rho*dPsi2dxi(k)*dPsi3dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi2(k)*Psi3(k) )
        M(i+2,i+2) = M(i+2,i+2) + wGauss(k)*( 2.0d0/dx*rho*dPsi3dxi(k)*dPsi3dxi(k) +                                    &
                                              0.75d0*dx/dt*dBdHk*Psi3(k)*Psi3(k) )
        
        ! Construction du residu
        R(i)        = R(i)   - wGauss(k)*( rho*dHkdx*dPsi1dxi(k) + 0.25d0*dx/dt*(3.0d0*Bk-4.0d0*Bk_1+Bk_2)*Psi1(k) )
        R(i+1)      = R(i+1) - wGauss(k)*( rho*dHkdx*dPsi2dxi(k) + 0.25d0*dx/dt*(3.0d0*Bk-4.0d0*Bk_1+Bk_2)*Psi2(k) )
        R(i+2)      = R(i+2) - wGauss(k)*( rho*dHkdx*dPsi3dxi(k) + 0.25d0*dx/dt*(3.0d0*Bk-4.0d0*Bk_1+Bk_2)*Psi3(k) )

      end do
      
      ! Symetrie
      M(i+1,i)   = M(i,i+1)
      M(i+2,i)   = M(i,i+2)
      M(i+2,i+1) = M(i+1,i+2)
       
    end do
    
    return
  end subroutine Assemblage
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine Solve (n , A,b , x,err)
    ! Description :
    ! Computes the solution of an n by n linear system A*x = b
    ! Code by John Burkardt
    ! This code is distributed under the GNU LGPL license.
    implicit none 
    ! Inputs :
    ! Dimension du systeme
    integer, intent(in)                                     :: n
    ! Inouts :
    ! Matrice
    double precision, dimension(n,n), intent(inout)         :: A
    ! Membre de droite
    double precision, dimension(n), intent(inout)           :: b
    ! Outputs :
    ! Solution
    double precision, dimension(n), intent(out)             :: x
    ! Indicateur d'erreur
    !    0, no error detected
    !    1, consistent singularity
    !    2, inconsistent singularity
    integer, intent(out)                                    :: err
    
    ! Variables :
    ! Erreur                                                
    character(len=maxStrLen)                                :: from, message
    ! Iterateurs
    integer                                                 :: i, j, k
    integer                                                 :: imax
    integer, dimension(n)                                   :: ipiv
    double precision                                        :: amax

    ! Intialisation
    from = "Solve"
    message = " "
    err = 0
    ipiv = 0
    x = 0.0d+00
    
    ! Process the matrix.
    do k = 1, n
  
      !  In column k :
      !  Seek the row imax with the properties that:
      !  imax has not already been used as a pivot;
      !  A(imax,k) is larger in magnitude than any other candidate.
      amax = 0.0d+00
      imax = 0
      do i = 1, n
        if ( ipiv(i) == 0 ) then
          if ( amax < abs(a(i,k)) ) then
            imax = i
            amax = abs(a(i,k))
          end if
        end if
      end do

      !  If you found a pivot row imax, then,
      !  eliminate the k-th entry in all rows that have not been used for pivoting.
      if ( imax /= 0 ) then
        ipiv(imax) = k
        a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
        b(imax) = b(imax) / a(imax,k)
        a(imax,k) = 1.0d+00

        do i = 1, n
          if ( ipiv(i) == 0 ) then
            a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
            b(i) = b(i) - a(i,k) * b(imax)
            a(i,k) = 0.0d+00
          end if
        end do

      end if

    end do

    !  Now, every row with nonzero IPIV begins with a 1, and all other rows are all zero.
    !  Begin solution.
    do j = n, 1, -1

      imax = 0
      do k = 1, n
        if ( ipiv(k) == j ) then
          imax = k
        end if
      end do

      if ( imax == 0 ) then

        x(j) = 0.0d+00

        if ( b(j) == 0.0d+00 ) then
          err = 1
          message = "System is consistent, but singular."
          call Error (from,message )
          return
        else
          err = 1
          message = "System is inconsistent and singular."
          call Error (from,message )
          return
        end if

      else

        x(j) = b(imax)

        do i = 1, n
          if ( i /= imax ) then
            b(i) = b(i) - a(i,j) * x(j)
          end if
        end do

      end if

    end do

    return
  end subroutine Solve
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculCurrentDensity  (ne,ndof,ntpc,dx,H_x_t , J_x_t)
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ne, ndof, ntpc
    double precision, intent(in)                            :: dx
    ! Solution aux dof
    double precision, dimension(ndof+1,ntpc), intent(in)    :: H_x_t
    ! Outputs :
    ! Solution aux pts en x
    double precision, dimension(ndof+1,ntpc), intent(out)   :: J_x_t
    
    ! Variables :
    integer                                                 :: e, i
    
    ! initialisation
    J_x_t = 0.0d0
    
    ! Traitement du 1er element
    J_x_t(ndof+1,1:ntpc)  = 2.0d0/dx*( dPsi1dxi_node(1)*H_x_t(ndof+1,1:ntpc)  + dPsi2dxi_node(1)*H_x_t(1,1:ntpc)    +   &
                                       dPsi3dxi_node(1)*H_x_t(2,1:ntpc)       )
    J_x_t(1,1:ntpc)       = 2.0d0/dx*( dPsi1dxi_node(2)*H_x_t(ndof+1,1:ntpc)  + dPsi2dxi_node(2)*H_x_t(1,1:ntpc)    +   &
                                       dPsi3dxi_node(2)*H_x_t(2,1:ntpc)       )
        
    ! Traitement du 2e element
    J_x_t(2,1:ntpc)       = 1.0d0/dx*( dPsi1dxi_node(3)*H_x_t(ndof+1,1:ntpc)  + dPsi2dxi_node(3)*H_x_t(1,1:ntpc)    +   &
                                       dPsi3dxi_node(3)*H_x_t(2,1:ntpc)       + dPsi1dxi_node(1)*H_x_t(2,1:ntpc)    +   &
                                       dPsi2dxi_node(1)*H_x_t(3,1:ntpc)       + dPsi3dxi_node(1)*H_x_t(4,1:ntpc)    )
    J_x_t(3,1:ntpc)       = 2.0d0/dx*( dPsi1dxi_node(2)*H_x_t(2,1:ntpc)       + dPsi2dxi_node(2)*H_x_t(3,1:ntpc)    +   &
                                       dPsi3dxi_node(2)*H_x_t(4,1:ntpc)       )
                                       
    ! Traitement des autres elements
    i = 2
    do e = 3, ne
      i = i + degPsi
      J_x_t(i,1:ntpc)     = 1.0d0/dx*( dPsi1dxi_node(3)*H_x_t(i-2,1:ntpc)     + dPsi2dxi_node(3)*H_x_t(i-1,1:ntpc)    +   &
                                       dPsi3dxi_node(3)*H_x_t(i,1:ntpc)       + dPsi1dxi_node(1)*H_x_t(i,1:ntpc)      +   &
                                       dPsi2dxi_node(1)*H_x_t(i+1,1:ntpc)     + dPsi3dxi_node(1)*H_x_t(i+2,1:ntpc)    )
      J_x_t(i+1,1:ntpc)   = 2.0d0/dx*( dPsi1dxi_node(2)*H_x_t(i,1:ntpc)       + dPsi2dxi_node(2)*H_x_t(i+1,1:ntpc)    +   &
                                       dPsi3dxi_node(2)*H_x_t(i+2,1:ntpc)   )
    end do
    
    return
  end subroutine CalculCurrentDensity 
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculLosses (ndof,ntpc,dx,dt,rho,H_x_t,B_x_t,J_x_t , Pj_x_t,Ph_x_t)
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ndof, ntpc
    double precision, intent(in)                            :: dx, dt
    ! Resistivite
    double precision, intent(in)                            :: rho
    ! Solution
    double precision, dimension(ndof+1,ntpc), intent(in)    :: H_x_t, B_x_t, J_x_t
    ! Outputs :
    ! Pertes
    double precision, dimension(ndof+1,ntpc), intent(out)   :: Pj_x_t, Ph_x_t
    
    ! Variables :
    integer                                                 :: n, nm1, np1
    
    ! initialisation
    Pj_x_t = 0.0d0
    Ph_x_t = 0.0d0
    
    ! Calcul des pertes
    do n = 1, ntpc
      
      ! Index pas avant et apres
      if ( n == 1 ) then
        nm1 = ntpc
        np1 = 2
      else if ( n == ntpc ) then
        nm1 = ntpc - 1
        np1 = 1
      else
        nm1 = n - 1
        np1 = n + 1
      end if
      
      ! Calcul des pertes Joules en fonction du temps
      Pj_x_t(1:ndof+1,n) = rho*J_x_t(1:ndof+1,n)*J_x_t(1:ndof+1,n)
      
      ! Calcul des pertes par hysteresis en fonction du temps
      Ph_x_t(1:ndof+1,n) = H_x_t(1:ndof+1,n)*0.5d0/dt*(B_x_t(1:ndof+1,np1)-B_x_t(1:ndof+1,nm1))
      
    end do
    
    return
  end subroutine CalculLosses
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine WriteResults (path,ndof,ntpc,dx,dt,H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t , err)
    implicit none
    ! Inputs :
    ! Path jusqu'au repertoire I/O
    character(len=maxStrLen), intent(in)                    :: path
    ! Discretisation
    integer, intent(in)                                     :: ndof, ntpc
    double precision, intent(in)                            :: dx, dt
    ! Champ mag, induction, densite de courant, pertes
    double precision, dimension(ndof+1,ntpc), intent(in)    :: H_x_t, B_x_t, J_x_t
    double precision, dimension(ndof+1,ntpc), intent(in)    :: Pj_x_t, Ph_x_t
    ! Output :
    ! Erreur
    integer, intent(out)                                    :: err
    
    ! Variables :
    ! Erreur
    character(len=maxStrLen)                                :: from, message
    ! Date et heure
    character(len=maxStrLen)                                :: date, time, zone
    integer, dimension(8)                                   :: values
    ! Fichiers
    character(len=maxStrLen)                                :: filename
    ! Format
    character(len=8)                                        :: f1dble
    character(len=43)                                       :: f6dble
    ! Iterateurs
    integer                                                 :: n, i
    
    ! Initialisation
    from = "WriteResults"
    f1dble = "(E23.15)"
    f6dble = "(E23.15,E23.15,E23.15,E23.15,E23.15,E23.15)"
    
    ! Nom du fichier
    filename = trim(path) // "Results.txt"
    
    ! Ouverture du fichier contenant les resultats pour une periode
    open (unit=1, file=trim(filename), iostat=err,  status="replace", action="write")
    if ( err /= 0 ) then
      message = "Error opening file Results.txt"
      call Error (from,message )
      return
    end if
    
    ! Recuperation de la date et de l'heure
    call Date_and_Time ( date,time,zone,values)
    write(1,"(a,a,a,a,a,a)") "Date  : ",date(7:8)," / ",date(5:6)," / ",date(1:4)
    write(1,"(a,a,a,a,a,a)") "Heure : ",time(1:2)," : ",time(3:4)," : ",time(5:6)
    
    ! Ecriture des dimensions
    write(1,*)
    write(1,"(a,i0)") "Nombre de pas de temps : ",ntpc
    write(1,"(a,i0)") "Nombre de points en x  : ",ndof+1
    
    ! Ecriture des pas de temps
    do n = 1, ntpc
      
      ! En-tete
      write(1,*)
      write(1,"(a,i0)")    "Pas no. ",n
      write(1,"(a,E23.15)") "Time =  ",dt*(n-1)
      write(1,"(a)")       "| x(m)                 | H(A/m)               | B(T)                 | J(A/m^2)             &
                           &| Pj(W/m^3)            | Ph(W/m^3)            "
      
      ! Ecriture des donnes pour le pas de temps courant
      write(1,f6dble) 0.0d0, H_x_t(ndof+1,n), B_x_t(ndof+1,n), J_x_t(ndof+1,n), Pj_x_t(ndof+1,n), Ph_x_t(ndof+1,n)
      do i = 1, ndof
        write(1,f6dble) 0.5d0*dx*dble(i), H_x_t(i,n), B_x_t(i,n), J_x_t(i,n), Pj_x_t(i,n), Ph_x_t(i,n)
      end do
      
    end do
    
    ! Fermeture du fichier
    close (1)
    
    return
  end subroutine WriteResults
  ! -------------------------------------------------------------------------------------------------------------------
  function Norm2 (n,x )
    ! Description :
    ! Calcul de la norme euclidienne d'un vecteur.
    implicit none
    ! Inputs :
    integer, intent(in)                                     :: n
    double precision, dimension(n), intent(in)              :: x
    ! Outputs :                                           
    double precision                                        :: norm2
                                                            
    ! Variables :                                           
    integer                                                 :: i
    
    ! Calcul de la norme Euclidienne
    norm2 = 0.0d0
    do i = 1, n
      norm2 = norm2 + x(i)*x(i)
    enddo
    norm2 = sqrt(norm2)/dble(n)
    
    return
  end function Norm2 
  ! -------------------------------------------------------------------------------------------------------------------
  
end module SlabProblemMod