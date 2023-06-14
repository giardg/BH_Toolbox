module MaterialMod
  ! -------------------------------------------------------------------------------------------------------------------
  ! Mai 2018
  ! Maxime Tousignant
  ! maxime.tousignant@polymtl.ca
  ! -------------------------------------------------------------------------------------------------------------------
  ! Description :
  ! Calcul de la loi de comportement des materiaux magnetiques et gestion de l'historique pour l'hysteresis.
  !
  ! Ce module permet de modeliser trois types de courbes BH :
  ! 1 : lineaire permeabilite relative constante
  !     bh(1) = mattype
  !     bh(2) = mur
  ! 2 : non lineaire modele arctangente
  !     bh(1) = mattype
  !     bh(2) = murmax
  !     bh(3) = bsat
  ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
  !     bh(1) = mattype
  !     bh(2) = natan
  !     bh(:) = ai
  !     bh(:) = bi
  !     bh(:) = ci 
  ! -------------------------------------------------------------------------------------------------------------------
  
  ! Modules
  use ElementMod
  use ErrorMod
  
  ! Constantes
  ! Pi
  double precision, parameter, private                      :: pi         = 3.1415926535897932d+00
  ! Permeabilite du vide
  double precision, parameter, private                      :: mu0        = 1.2566370614359172d-06
  
  ! Nombre maximal de parametres pour la courbe BH
  integer, parameter                                        :: maxParam   = 256
  integer, parameter                                        :: maxAtan    = 5
  
  ! Nombre de points total
  integer, private                                          :: nPtsTot
  
  ! Historique de Preisach
  integer, parameter, private                               :: maxHistLen = 256
  ! Nombre de points dans l'historique
  integer, dimension(:), allocatable, private               :: nHist
  ! Contient prevG et partialM
  double precision, dimension(:,:), allocatable, private    :: dHist
  ! Contient les alpha
  double precision, dimension(:,:), allocatable, private    :: dHistA
  ! Contient les beta
  double precision, dimension(:,:), allocatable, private    :: dHistB
  
  ! Routine et fonctions privees
  private                                                   :: GetPreisachParameters
  private                                                   :: CalculPreisachBH
  private                                                   :: CalculFG
  private                                                   :: CalculFG4param
  private                                                   :: Everett
  private                                                   :: dEverettda
  private                                                   :: dEverettdb
  private                                                   :: EverettEllipse
  private                                                   :: dEverettEllipseda
  private                                                   :: dEverettEllipsedb
  
  contains
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine InitializeMaterial (ne , bh , err)
    ! Description :
    ! Calcul de parametres supplementaires necessaires aux modeles.
    ! Initialisation de l'historique pour le modele de Preisach.
    implicit none
    ! Inputs :
    ! Nombre d'elements
    integer, intent(in)                                     :: ne
    ! Inout :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(inout)    :: bh
    ! Indicateur d'erreur
    integer, intent(out)                                    :: err
    
    ! Variables :
    Integer                                                 :: i, ptr
    integer                                                 :: mattype
    ! Modele de Preisach
    integer                                                 :: natan
    double precision, dimension(maxAtan)                    :: ai, bi, ci
    double precision                                        :: muFmax, bsat, brem

    !Ellipse
    double precision                                        :: a, b, alpha, Hmax

    ! Erreur
    character(len=maxStrLen)                                :: from, message
    
    
    ! initialisation
    from = "InitializeMaterial"
    message = " "
    
    ! Type de courbe BH
    mattype = nint(bh(1))
    
    ! Nombre de pts total
    nPtsTot = ne*(nGauss+degPsi) + 1
    
    ! Verification des parametres de base et calcul des parametres supplementaires
    ! 1 : lineaire permeabilite relative constante
    if ( mattype == 1 ) then
      if ( bh(2) < 1.0d0 ) then
        err = 1
        message = "mur doit etre > 1."
        call Error (from,message )
        return
      end if
    
    ! 2 : non lineaire modele arctangente
    else if ( mattype == 2 ) then
      if ( bh(2) < 1.0d0 ) then
        err = 1
        message = "murmax doit etre > 1."
        call Error (from,message )
        return
      end if
      if ( bh(3) < epsilon(1.0d0) ) then
        err = 1
        message = "Bsat doit etre > 0."
        call Error (from,message )
        return
      end if
      bh(4) = 2.0d0*bh(3)/pi
      bh(5) = 0.5d0*pi*mu0*(bh(2)-1.0d0)/bh(3)
    
    ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients  
    else if ( mattype == 3 ) then
      
      ! Recuperation de ntanh, ai, bi, ci
      call GetPreisachParameters (bh , natan,ai,bi,ci,muFmax,brem,bsat)
      
      brem = 0.0d0
      bsat = 0.0d0
      do i = 1, natan
        if ( bi(i) < epsilon(1.0d0) ) then
          err = 1
          message = "Les coefficients bi doivent etre > 0."
          call Error (from,message )
          return
        end if
        brem = brem + ai(i)*atan(ci(i)/bi(i))
        bsat = bsat + ai(i)*0.5d0*pi
      end do
      
      if ( brem < 0.0d0 ) then
        brem = -brem
        do i = 1, natan
          ci(i) = -ci(i)
        end do
      end if
      bh((2*natan+3):(3*natan+2)) = ci(1:natan)
      
      if ( bsat < brem ) then
        err = 1
        message = "Bsat doit etre > Brem."
        call Error (from,message )
        return
      end if
      
      ! Position initiale dans le vecteur bh
      ptr = 3*natan + 2
      
      ! Calcul de muFmax
      muFmax = 0.0d0
      do i = 1, natan
        muFmax = muFmax + ai(i)/bi(i)/(1.0d0+(ci(i)/bi(i))**2.0d0)
      end do
      ptr = ptr + 1
      bh(ptr) = muFmax
      
      ! Ajout de Brem et Bsat
      ptr = ptr + 1
      bh(ptr) = brem
      ptr = ptr + 1
      bh(ptr) = bsat
      
      ! Allocation de l'historique
      Allocate ( nHist(nPtsTot), stat=err )
      if ( err == 0 ) Allocate (dHist(2,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistA(maxHistLen,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistB(maxHistLen,nPtsTot) , stat=err)
      if ( err /= 0 ) then
        message = "Memoire insuffisante."
        call Error (from,message )
        return
      end if
      
      ! Initialisation de l'historique
      nHist = 1
      dHist = 0.0d0
      dHistA = 0.0d0
      dHistB = 0.0d0
    
    elseif ( mattype == 4 ) then
      if ( bh(2) < 0.0d0 ) then
        err = 1
        message = "Br doit etre > 0."
        call Error (from,message )
        return
      end if

    else if ( mattype == 5 ) then
      
      ! Recuperation de ntanh, ai, bi, ci
      Hmax = bh(4)
      call GetEllipseParameters (bh,Hmax , a,b,alpha,brem)
      
      
      ! Allocation de l'historique
      Allocate ( nHist(nPtsTot), stat=err )
      if ( err == 0 ) Allocate (dHist(2,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistA(maxHistLen,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistB(maxHistLen,nPtsTot) , stat=err)
      if ( err /= 0 ) then
        message = "Memoire insuffisante."
        call Error (from,message )
        return
      end if
      
      ! Initialisation de l'historique
      nHist = 1
      dHist = 0.0d0
      dHistA = 0.0d0
      dHistB = 0.0d0

    else if ( mattype == 6 ) then
      
      ! Recuperation de ntanh, ai, bi, ci
      Hmax = bh(4)
      call GetEllipseParameters (bh,Hmax, a,b,alpha,brem)
      
      
      ! Allocation de l'historique
      Allocate ( nHist(nPtsTot), stat=err )
      if ( err == 0 ) Allocate (dHist(2,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistA(maxHistLen,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistB(maxHistLen,nPtsTot) , stat=err)
      if ( err /= 0 ) then
        message = "Memoire insuffisante."
        call Error (from,message )
        return
      end if
      
      ! Initialisation de l'historique
      nHist = 1
      dHist = 0.0d0
      dHistA = 0.0d0
      dHistB = 0.0d0

    else if ( mattype == 7 ) then
      
      ! Recuperation de ntanh, ai, bi, ci
      brem = bh(2)
      bsat = bh(3)
      
      ! Allocation de l'historique
      Allocate ( nHist(nPtsTot), stat=err )
      if ( err == 0 ) Allocate (dHist(2,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistA(maxHistLen,nPtsTot) , stat=err)
      if ( err == 0 ) Allocate (dHistB(maxHistLen,nPtsTot) , stat=err)
      if ( err /= 0 ) then
        message = "Memoire insuffisante."
        call Error (from,message )
        return
      end if
      
      ! Initialisation de l'historique
      nHist = 1
      dHist = 0.0d0
      dHistA = 0.0d0
      dHistB = 0.0d0

    
      
    end if
    
    return
  end subroutine InitializeMaterial
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine FinalizeMaterial ()
    ! Description :
    ! Libere la memoire allouee pour l'historique
    implicit none

    if ( allocated(nHist) ) deallocate(nHist)
    if ( allocated(dHist) ) deallocate(dHist)
    if ( allocated(dHistA) ) deallocate(dHistA)
    if ( allocated(dHistB) ) deallocate(dHistB)
    
    return
  end subroutine FinalizeMaterial
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine UpdateMaterialBH (ne,ndof,bh,H , B)
    ! Description :
    ! Calcul de B en post traitement et mise-a-jour de l'historique sur tout le domaine.
    implicit none
    ! Inputs :
    ! Discretisation
    integer, intent(in)                                     :: ne, ndof
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, dimension(ndof+1), intent(in)         :: H
    ! Outputs :
    ! Induction
    double precision, dimension(ndof+1), intent(out)        :: B
    
    ! Variables :
    ! Iterations
    integer                                                 :: e, i, j, k, noPt
    double precision, dimension(nGauss)                     :: Hk
    double precision                                        :: Bk
    
    ! Initialisation
    i = 0
    
    ! Premier element
    
    ! Noeuds
    noPt = GetNoPt (1,1,1 )
    call UpdateLocalBH (noPt,bh,H(ndof+1) , B(ndof+1))
    
    noPt = GetNoPt (1,2,1 )
    call UpdateLocalBH (noPt,bh,H(1) , B(1))
    
    noPt = GetNoPt (1,3,1 )
    call UpdateLocalBH (noPt,bh,H(2) , B(2))
    
    ! Pts de calcul
    do k = 1, nGauss
      
      ! Numero identifiant le pt de calcul
      noPt = GetNoPt(1,k,2 )
      
      ! Champ mag
      Hk(k) = H(ndof+1)*Psi1(k)+H(i+1)*Psi2(k)+H(i+2)*Psi3(k)
      
      ! Mise-a-jour de l'historique
      call UpdateLocalBH (noPt,bh,Hk(k) , Bk)
      
    end do
    
    ! Boucle sur les autres elements
    do e = 2, ne
      
      ! Pointeur pour les degrees de liberte
      i = i + degPsi
      
      ! Noeuds
      noPt = GetNoPt (e,2,1 )
      call UpdateLocalBH (noPt,bh,H(i+1) , B(i+1))
      
      noPt = GetNoPt (e,3,1 )
      call UpdateLocalBH (noPt,bh,H(i+2) , B(i+2))
      
      ! Pts de calcul
      do k = 1, nGauss
      
        ! Numero identifiant le pt de calcul
        noPt = GetNoPt(e,k,2 )
      
        ! Champ mag
        Hk(k) = H(i)*Psi1(k)+H(i+1)*Psi2(k)+H(i+2)*Psi3(k)
      
        ! Mise-a-jour de l'historique
        call UpdateLocalBH (noPt,bh,Hk(k) , Bk)
      
      end do
      
    end do
    
    return
  end subroutine UpdateMaterialBH
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine UpdateLocalBH (noPt,bh,H , B)
    ! Description :
    ! Calcul de B en post traitement et mise-a-jour de l'historique en un point donne.
    implicit none
    ! Inputs :
    ! Numero identifiant le point
    integer, intent(in)                                     :: noPt
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    ! Outputs :
    ! Induction
    double precision, intent(out)                           :: B
    
    ! Variables :
    ! Numero du modele de courbe BH
    integer                                                 :: mattype
    double precision                                        :: dBdH
    
    ! Initialisation
    B = 0.0d0
    
    ! Modele du materiau
    mattype = nint(bh(1))
    
    ! 1 : lineaire permeabilite relative constante
    if ( mattype == 1 ) then
      B = mu0*bh(2)*H
    
    ! 2 : non lineaire modele arctangente
    else if ( mattype == 2 ) then
      B = mu0*H + bh(4)*atan(bh(5)*H)
    
    ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients 
    else if ( mattype == 3 ) then
      call CalculPreisachBH (bh,H , nHist(noPt),dHist(1,noPt),dHist(2,noPt),                                            &
                             dHistA(1:maxHistLen,noPt),dHistB(1:maxHistLen,noPt) , B,dBdh)

  else if ( mattype == 4 ) then
    
    if (H < -1.0*bh(3)) then
      B = mu0*H - bh(2)
    else if (H > bh(3)) then
      B = mu0*H + bh(2)
    else
      B = (-bh(2)/(2.0d0*bh(3)**3.0d0))*(H**3.0d0)+(mu0+(3.0d0*bh(2))/(2.0d0*bh(3)))*H
    end if

  else if ( mattype == 5 ) then
    call CalculPreisachBH (bh,H , nHist(noPt),dHist(1,noPt),dHist(2,noPt),                                            &
                             dHistA(1:maxHistLen,noPt),dHistB(1:maxHistLen,noPt) , B,dBdh)
  
  else if ( mattype == 6 ) then
    call CalculPreisachBHEllipse (bh,H , nHist(noPt),dHist(1,noPt),dHist(2,noPt),                                            &
                              dHistA(1:maxHistLen,noPt),dHistB(1:maxHistLen,noPt) , B,dBdh)

  else if ( mattype == 7 ) then
    call CalculPreisachBH (bh,H , nHist(noPt),dHist(1,noPt),dHist(2,noPt),                                            &
                              dHistA(1:maxHistLen,noPt),dHistB(1:maxHistLen,noPt) , B,dBdh)
      
  end if
    
  return
  end subroutine UpdateLocalBH
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculLocalBH (noPt,bh,H , B,dBdH)
    ! Description :
    ! Calcul de B et de sa derivee en un point donne.
    implicit none
    ! Inputs :
    ! Numero identifiant le point
    integer, intent(in)                                     :: noPt
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    ! Outputs :
    ! Induction et sa derivee
    double precision, intent(out)                           :: B
    double precision, intent(out)                           :: dBdH
    
    ! Variables :
    ! Numero du modele de courbe BH
    integer                                                 :: mattype
    ! Modele de Preisach
    integer                                                 :: n
    double precision                                        :: prevG, partialM, prevH
    double precision, dimension(maxHistLen)                 :: prevMax, prevMin
    
    ! Initialisation
    B = 0.0d0
    dBdH = 0.0d0
    
    ! Modele du materiau
    mattype = nint(bh(1))
    
    ! 1 : lineaire permeabilite relative constante
    if ( mattype == 1 ) then
      B = mu0*bh(2)*H
      dBdH = mu0*bh(2)
      
    ! 2 : non lineaire modele arctangente
    else if ( mattype == 2 ) then
      B = mu0*H + bh(4)*atan(bh(5)*H)
      dBdH = mu0 + mu0*(bh(2)-1.0d0)/(1.0d0 + (bh(5)*H)**2.0d0)
      
    ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients 
    else if ( mattype == 3 ) then
      n = nHist(noPt)
      prevG = dHist(1,noPt)
      partialM = dHist(2,noPt)
      prevMax = dHistA(1:maxHistLen,noPt)
      prevMin = dHistB(1:maxHistLen,noPt)
      call CalculPreisachBH (bh,H , n,prevG,partialM,prevMax,prevMin , B,dBdh)

    else if ( mattype == 4 ) then
      if (H < -1.0*bh(3)) then
        B = mu0*H - bh(2)
        dBdH = mu0
      else if (H > bh(3)) then
        B = mu0*H + bh(2)
        dBdH = mu0
      else
        B = (-bh(2)/(2.0d0*bh(3)**3.0d0))*(H**3.0d0)+(mu0+(3.0d0*bh(2))/(2.0d0*bh(3)))*H
        dBdH = (-3.0d0*bh(2)/(2.0d0*bh(3)**3.0d0))*(H**2.0d0)+(mu0+(3.0d0*bh(2))/(2.0d0*bh(3)))
      end if

    else if ( mattype == 5 ) then
      n = nHist(noPt)
      prevG = dHist(1,noPt)
      partialM = dHist(2,noPt)
      prevMax = dHistA(1:maxHistLen,noPt)
      prevMin = dHistB(1:maxHistLen,noPt)
      call CalculPreisachBH (bh,H , n,prevG,partialM,prevMax,prevMin , B,dBdh)

    else if ( mattype == 6 ) then
      n = nHist(noPt)
      prevH = dHist(1,noPt)
      partialM = dHist(2,noPt)
      prevMax = dHistA(1:maxHistLen,noPt)
      prevMin = dHistB(1:maxHistLen,noPt)
      call CalculPreisachBHEllipse (bh,H , n,prevH,partialM,prevMax,prevMin , B,dBdh)

    else if ( mattype == 7 ) then
      n = nHist(noPt)
      prevG = dHist(1,noPt)
      partialM = dHist(2,noPt)
      prevMax = dHistA(1:maxHistLen,noPt)
      prevMin = dHistB(1:maxHistLen,noPt)
      call CalculPreisachBH (bh,H , n,prevG,partialM,prevMax,prevMin , B,dBdh)
      
    end if
    
    return
  end subroutine CalculLocalBH
  ! -------------------------------------------------------------------------------------------------------------------
  function GetNoPt (e,k,type )
    ! Description :
    ! Recupere le numero unique identifiant un noeud ou un point de calcul.
    implicit none
    ! Inputs :
    ! Numero de l'element
    integer, intent(in)                                     :: e
    ! Numero du point (calcul ou neud)                      
    integer, intent(in)                                     :: k
    ! Type de point (1 pour noeud, 2 pour gauss)            
    integer, intent(in)                                     :: type
    ! Output :                                              
    integer                                                 :: getNoPt

    ! Calcul du numero du point
    ! Degres de liberte
    if ( type == 1 ) then
      getNoPt = 6*(e-1) + 3*(k-1) + 1
      
    ! Pts d'integration
    else
      getNoPt = 6*(e-1) + k + 1
      if ( k >= 3 ) getNoPt = getNoPt + 1
      
    end if
    
    return
  end function GetNoPt
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine GetPreisachParameters (bh , natan,ai,bi,ci,muFmax,brem,bsat)
    ! Description :
    ! Recuperation des parametres a partir du vecteur bh.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Outputs :
    ! Nombre de tangentes hyperboliques
    integer, intent(out)                                    :: natan
    ! Parametres des tangentes hyperboliques
    double precision, dimension(maxAtan), intent(out)       :: ai, bi, ci
    ! Parametres supplementaires
    double precision, intent(out)                           :: muFmax, brem, bsat
    
    ! Variables :
    integer                                                 :: i, ptr
    
    ! Initialisation
    ptr = 2
    
    ! Recuperation du nombre de tangentes hyperboliques
    natan = nint(bh(ptr))
    
    ! Recuperation des ai
    do i = 1, natan
      ptr = ptr + 1
      ai(i) = bh(ptr)
    end do
    
    ! Recuperation des bi
    do i = 1, natan
      ptr = ptr + 1
      bi(i) = bh(ptr)
    end do
    
    ! Recuperation des ci
    do i = 1, natan
      ptr = ptr + 1
      ci(i) = bh(ptr)
    end do
    
    ! Recuperation de muFmax, Brem et Bsat
    ptr = ptr + 1
    muFmax = bh(ptr)
    ptr = ptr + 1
    brem = bh(ptr)
    ptr = ptr + 1
    bsat = bh(ptr)
    
    return
  end subroutine GetPreisachParameters
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculPreisachBH (bh,H , n,prevG,partialM,prevMax,prevMin, B,dBdH)
    ! Description :
    ! Calcul de B et sa derivee avec le modele de Preisach a 3n coefficients. Le calcul modifie l'historique.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ magnetique
    double precision, intent(in)                            :: H
    ! Inouts :
    ! Nombre de points dans l'historique
    integer, intent(inout)                                  :: n
    ! Valeur de G au pas de temps precedent
    double precision, intent(inout)                         :: prevG
    ! Valeur partielle de la magnetisation
    double precision, intent(inout)                         :: partialM
    ! Maxima precedents de G
    double precision, dimension(maxHistLen), intent(inout)  :: prevMax
    ! Minima precedents de G
    double precision, dimension(maxHistLen), intent(inout)  :: prevMin
    ! Outputs :
    ! Induction et sa derivee
    double precision, intent(out)                           :: B
    double precision, intent(out)                           :: dBdH
    
    ! Variables :
    integer                                                 :: j
    ! Modele de Preisach
    integer                                                 :: natan, prevn
    double precision                                        :: brem, F, dFdH, G, dGdH, Mirr, dMirrdG
    double precision                                        :: Hmax,Bdown,dBdowndH

    ! Initialisation
    B = mu0*H
    dBdH = mu0
    Mirr = 0.0d0
    dMirrdG = 0.0d0
    
    ! Calcul des fonctions F(H) et G(H) et recuperation de brem
    if (nint(bh(1)) == 3) then
      call CalculFG (bh,H , brem,F,dFdH,G,dGdH)
    else if (nint(bh(1)) == 5) then
      Hmax = bh(4)
      call CalculFGEllipse (bh,H,Hmax , brem,F,dFdH,G,dGdH,Bdown,dBdowndH)
    else if (nint(bh(1)) == 7) then
      call CalculFG4param (bh,H , brem,F,dFdH,G,dGdH)
    end if

    ! Calcul du modele de Preisach si existence d'une partie irreversible
    if ( brem > epsilon(1.0d0) ) then
      
      ! Mise-a-jour de la frontiere de Preisach
      ! Copie du nombre de points dans l'historique
      prevn = n
 
      ! Champ croissant
      if ( G > prevG ) then
        do j = prevn, 0, -1
          ! Saturation
          if ( j == 0 ) then
            n = 1
            prevMax(1) =  abs(G)
            prevMin(1) = -abs(G)
          ! Wipe-out ou creation d'un nouvel extremum
          elseif ( G < prevMax(j) ) then
            n = j+1
            prevMax(n) = G
            prevMin(n) = prevMin(j)
            exit
          endif
        enddo

      ! Champ decroissant
      else if ( G < prevG  ) then
        do j =  prevn, 0, -1
          ! Saturation
          if ( j == 0 ) then
            n = 1
            prevMax(1) =  abs(G)
            prevMin(1) = -abs(G)
          ! Wipe-out ou creation d'un nouvel extremum
            else if ( G > prevMin(j) ) then
            n = j+1
            prevMax(n) = prevMax(j)
            prevMin(n) = G
            exit
          endif
        enddo
      
      endif
      
      ! Calcul de la somme partielle si necessaire
      if ( prevn /= n ) then
        partialM = 0.0d0
        do j = 2, n-1
          if( prevMin(j-1) /= prevMin(j) )then
            partialM = partialM + 2.0d0*( Everett(prevMax(j),prevMin(j-1) ) - Everett(prevMax(j),prevMin(j)) )
          endif
        enddo
      endif
    
      ! Calcul de la contribution irreversible a l'induction
      ! Termes essentiels
      Mirr = -Everett(prevMax(1),prevMin(1) ) + 2.0d0*Everett(G,prevMin(n) )

      ! Ajout des termes supplementaires
      if( n > 1 )then
        Mirr = Mirr + partialM
        if( prevMin(n-1) /= prevMin(n) )then
          Mirr = Mirr + 2.0d0*( Everett(prevMax(n),prevMin(n-1) ) - Everett(prevMax(n),prevMin(n) ) )
        end if
      end if
      
      ! Calcul de la derivee de la partie irreversible
      ! Champ croissant
      if ( G > prevG ) then
        dMirrdG =  2.0d0*dEverettda(G,prevMin(n) )
        
      ! Champ decroissant
      else if ( G < prevG ) then
        dMirrdG = -2.0d0*dEverettdb(prevMax(n),G )
         
      ! Champ constant
      else
        dMirrdG = dEverettda(G,prevMin(n) ) - dEverettdb(prevMax(n),G )
      
      end if
      
      ! Normalisation
      Mirr = Mirr / brem
      dMirrdG = dMirrdG / brem
      
      ! Mise-a-jour de la valeur precedente de G
      prevG = G
      
    end if
    
      
    ! Addition des composantes reversibles et irreversibles
    ! Induction
    B = B + F + Mirr
      
    ! Permeabilite differentielle
    dBdH = dBdH + dFdH + dMirrdG*dGdH

    !write(*,*) B, dBdH, brem
    
    return
  end subroutine CalculPreisachBH
! -------------------------------------------------------------------------------------------------------------------
  subroutine CalculFG (bh,H , brem,F,dFdH,G,dGdH)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    ! Outputs :
    ! Induction a remanence
    double precision, intent(out)                           :: brem
    ! Valeur de la fonction F et sa derivee
    double precision, intent(out)                           :: F, dFdH
    ! Valeur de la fonction G et sa derivee
    double precision, intent(out)                           :: G, dGdH
      
    ! Variables :
    integer                                                 :: i, natan
    double precision, dimension(maxAtan)                    :: ai, bi, ci
    double precision                                        :: muFmax, bsat
    double precision                                        :: xia, xid
    
    ! Recuperation des parametres
    call GetPreisachParameters (bh , natan,ai,bi,ci,muFmax,brem,bsat)
        
    ! Cas H = 0
    if ( abs(H) < epsilon(1.0d0) ) then
      F = 0.0d0
      dFdH = muFmax
      G = 0.0d0
      dGdH = 0.0d0
      
    ! Cas H < 0
    else if ( H < 0.0d0 ) then
      
      ! Initialisation
      F = brem
      dFdH = 0.0d0
      G = -brem
      dGdH = 0.0d0
      
      ! Boucle sur le nombre de tangentes hyperboliques
      do i = 1, natan
        
        ! Calcul de la fonction reversible
        xid = (-H+ci(i))/bi(i)
        F = F - ai(i)*atan(xid)
        dFdH = dFdH + ai(i)/bi(i)/(1.0d0+xid*xid)
        
        ! Calcul de la fonction irreversible
        xia = (-H-ci(i))/bi(i)
        G  = G + 0.5d0*ai(i)*(atan(xid)-atan(xia))
        dGdH = dGdH - 0.5d0*ai(i)/bi(i)*(1.0d0/(1.0d0+xid*xid)-1.0d0/(1.0d0+xia*xia))
        
      end do
      
    ! Cas H > 0
    else
      
      ! Initialisation
      F = -brem
      dFdH = 0.0d0
      G = brem
      dGdH = 0.0d0
      
      ! Boucle sur le nombre de tangentes hyperboliques
      do i = 1, natan
      
        ! Calcul de la fonction reversible
        xid = (H+ci(i))/bi(i)
        F = F + ai(i)*atan(xid)
        dFdH = dFdH + ai(i)/bi(i)/(1.0d0+xid*xid)
      
        ! Calcul de la fonction irreversible
        xia = (H-ci(i))/bi(i)
        G  = G - 0.5d0*ai(i)*(atan(xid)-atan(xia))
        dGdH = dGdH - 0.5d0*ai(i)/bi(i)*(1.0d0/(1.0d0+xid*xid)-1.0d0/(1.0d0+xia*xia))
      
      end do
        
    end if
      
    return
  end subroutine CalculFG
  ! -------------------------------------------------------------------------------------------------------------------
  function Everett(alpha,beta )
    ! Description :
    ! Calcul de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    ! Outputs :
    double precision                                        :: Everett
    
    Everett = -alpha*beta
    if ( Everett < 0.0d0 ) Everett = 0.0d0
    
    return
  end function Everett
  ! -------------------------------------------------------------------------------------------------------------------
  function dEverettda(alpha,beta )
    ! Description :
    ! Calcul de la derivee de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    ! Outputs :
    double precision                                        :: dEverettda
    
    if ( alpha*beta < 0.0d0 ) then
      dEverettda = -beta
    else
      dEverettda = 0.0d0
    end if
    
    return
  end function dEverettda
  ! -------------------------------------------------------------------------------------------------------------------
  function dEverettdb(alpha,beta )
    ! Description :
    ! Calcul de la derivee de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    ! Outputs :
    double precision                                        :: dEverettdb
    
    if ( alpha*beta < 0.0d0 ) then
      dEverettdb = -alpha
    else
      dEverettdb = 0.0d0
    end if
    
    return
  end function dEverettdb
  ! ------------------------------------------------------------------------------------------------------------------- 

  subroutine GetEllipseParameters (bh,Hmax_in , a,b,alpha,brem)
    ! Description :
    ! Recuperation des parametres de l'ellipse (a,b,alpha) a partir du vecteur bh.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    double precision, intent(in)                            :: Hmax_in
    ! Outputs :
    ! Parametres supplementaires
    double precision, intent(out)                           :: a, b, alpha, brem

    !Autres
    double precision                                        :: delta, Hmax, Bmax, Hc, bnum, bdenom

    Hmax = abs(Hmax_in)
    delta = atan(-bh(3)/bh(2))
    Bmax = mu0*sqrt(bh(2)**2.0d0 + bh(3)**2.0d0)*Hmax
    Hc = Hmax*sin(delta)
    alpha = atan(Bmax/Hmax)*cos(delta*0.9d0)

    bnum = ((Hmax*sin(alpha)-Bmax*cos(delta)*cos(alpha))**2.0d0-(Hmax*sin(alpha)+Bmax*cos(delta)*sin(alpha)*tan(alpha))**2.0d0)
    bdenom = (1-((Hmax+Bmax*cos(delta)*tan(alpha))/Hc)**2.0d0)
    b = sqrt(bnum/bdenom)
    a = sqrt(((Hc*cos(alpha))**2.0d0)/(1-(Hc*sin(alpha)/b)**2.0d0))
    brem = sqrt(1.0d0/((sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0))
    
    return
  end subroutine GetEllipseParameters
  ! ------------------------------------------------------------------------------------------------------------------- 

  subroutine CalculFGEllipse (bh,H,Hmax , brem,F,dFdH,G,dGdH,Bdown,dBdowndH)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    double precision, intent(in)                            :: Hmax
    ! Outputs :
    ! Induction a remanence
    double precision, intent(out)                           :: brem
    ! Valeur de la fonction F et sa derivee
    double precision, intent(out)                           :: F, dFdH
    ! Valeur de la fonction G et sa derivee
    double precision, intent(out)                           :: G, dGdH
    ! Valeur de la fonction B down et sa derivee
    double precision, intent(out)                           :: Bdown, dBdowndH
      
    ! Variables :
    double precision                                        :: a, b, alpha
    double precision                                        :: A2, C, D, dCdH, dDdH, Bup, dBupdH,discr
    ! Recuperation des parametres
    call GetEllipseParameters (bh,Hmax , a,b,alpha,brem)
    !write(*,*) a, b, alpha, brem
    
    if (abs(Hmax) < epsilon(1.0d0) ) then
      brem = 0.0d0
      F = 0.0d0
      dFdH = 0.0d0
      G = 0.0d0 
      dGdH = 0.0d0 
      Bdown = 0.0d0
      dBdowndH = 0.0d0
    else


      if ( H < 0.0d0 ) then
      
        A2 = (sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0
        C = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*abs(H)
        D = ((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(H**2.0d0) - 1.0d0
        dCdH = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))
        dDdH = 2.0d0*((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*abs(H)
  
        discr = C**2.0d0 - 4.0d0*A2*D


        Bdown = (-C + sqrt(discr))/(2.0d0*A2)
        Bup = (-C - sqrt(discr))/(2.0d0*A2)
        dBdowndH = (-dCdH + (2.0d0*dCdH*C - 4.0d0*A2*dDdH)/(2.0d0*sqrt(discr)))/(2.0d0*A2)
        dBupdH = (-dCdH - (2.0d0*dCdH*C - 4.0d0*A2*dDdH)/(2.0d0*sqrt(discr)))/(2.0d0*A2)
      
        F = -(Bdown - brem -mu0*abs(H))
        G = -(brem - (Bdown-Bup)/2.0d0)
        dFdH = (dBdowndH - mu0)
        dGdH = -(dBdowndH-dBupdH)/2.0d0
      

      ! Cas H > 0
      else
      
        A2 = (sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0
        C = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(H)
        D = ((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(H**2.0d0) - 1.0d0
        dCdH = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))
        dDdH = 2.0d0*((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(H)
      
        discr = C**2.0d0 - 4.0d0*A2*D


        Bdown = (-C + sqrt(discr))/(2.0d0*A2)
        Bup = (-C - sqrt(discr))/(2.0d0*A2)
        dBdowndH = (-dCdH + (2.0d0*dCdH*C - 4.0d0*A2*dDdH)/(2.0d0*sqrt(discr)))/(2.0d0*A2)
        dBupdH = (-dCdH - (2.0d0*dCdH*C - 4.0d0*A2*dDdH)/(2.0d0*sqrt(discr)))/(2.0d0*A2)
    

        F = Bdown - brem -mu0*(H)
        G = brem - (Bdown-Bup)/2.0d0
        dFdH = (dBdowndH - mu0)
        dGdH = -(dBdowndH-dBupdH)/2.0d0


        
      end if

    end if
    !write(*,*) F, G, dFdH, dGdH
    return
  end subroutine CalculFGEllipse
  ! -------------------------------------------------------------------------------------------------------------------


  subroutine CalculPreisachBHEllipse (bh,H , n,prevH,partialM,prevMax,prevMin, B,dBdH)
    ! Description :
    ! Calcul de B et sa derivee avec le modele de Preisach a 3n coefficients. Le calcul modifie l'historique.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ magnetique
    double precision, intent(in)                            :: H
    ! Inouts :
    ! Nombre de points dans l'historique
    integer, intent(inout)                                  :: n
    ! Valeur de G au pas de temps precedent
    double precision, intent(inout)                         :: prevH
    ! Valeur partielle de la magnetisation
    double precision, intent(inout)                         :: partialM
    ! Maxima precedents de G
    double precision, dimension(maxHistLen), intent(inout)  :: prevMax
    ! Minima precedents de G
    double precision, dimension(maxHistLen), intent(inout)  :: prevMin
    double precision                                        :: E1, E2, dE1, dE2
    ! Outputs :
    ! Induction et sa derivee
    double precision, intent(out)                           :: B
    double precision, intent(out)                           :: dBdH
    
    ! Variables :
    integer                                                 :: j
    ! Modele de Preisach
    integer                                                 :: natan, prevn
    double precision                                        :: brem, F, dFdH, G, dGdH, Mirr, dMirrdH
    double precision                                        :: Hmax,Bdown,dBdowndH

    ! Initialisation
    B = mu0*H
    dBdH = mu0
    Mirr = 0.0d0
    
    ! Calcul des fonctions F(H) et G(H) et recuperation de brem
    call CalculFGEllipse (bh,H,bh(4) , brem,F,dFdH,G,dGdH,Bdown,dBdowndH)
    

    ! Calcul du modele de Preisach si existence d'une partie irreversible
    if ( brem > epsilon(1.0d0) ) then
      
      ! Mise-a-jour de la frontiere de Preisach
      ! Copie du nombre de points dans l'historique
      prevn = n
 
      ! Champ croissant
      if ( H > prevH ) then
        do j = prevn, 0, -1
          ! Saturation
          if ( j == 0 ) then
            n = 1
            prevMax(1) =  abs(H)
            prevMin(1) = -abs(H)
          ! Wipe-out ou creation d'un nouvel extremum
          elseif ( H < prevMax(j) ) then
            n = j+1
            prevMax(n) = H
            prevMin(n) = prevMin(j)
            exit
          endif
        enddo

      ! Champ decroissant
      else if ( H < prevH  ) then
        do j =  prevn, 0, -1
          ! Saturation
          if ( j == 0 ) then
            n = 1
            prevMax(1) =  abs(H)
            prevMin(1) = -abs(H)
          ! Wipe-out ou creation d'un nouvel extremum
            else if ( H > prevMin(j) ) then
            n = j+1
            prevMax(n) = prevMax(j)
            prevMin(n) = H
            exit
          endif
        enddo
      
      endif
      
      prevH = H

      ! Calcul de la somme partielle si necessaire
      if ( prevn /= n ) then
        partialM = 0.0d0
        do j = 2, n-1
          if( prevMin(j-1) /= prevMin(j) )then
            call EverettEllipse(prevMax(j),prevMin(j-1),bh, E1 )
            call EverettEllipse(prevMax(j),prevMin(j),bh, E2)
            partialM = partialM + 2.0d0*( E1 - E2)
          endif
        enddo
      endif
    
      ! Calcul de la contribution irreversible a l'induction
      ! Termes essentiels
      call EverettEllipse(prevMax(1),prevMin(1),bh, E1  )
      call EverettEllipse(H,prevMin(n),bh, E2  )
      Mirr = -E1 + 2.0d0*E2

      ! Ajout des termes supplementaires
      if( n > 1 )then
        Mirr = Mirr + partialM
        if( prevMin(n-1) /= prevMin(n) )then
          call EverettEllipse(prevMax(n),prevMin(n-1) ,bh, E1 )
          call EverettEllipse(prevMax(n),prevMin(n),bh, E2 )
          Mirr = Mirr + 2.0d0*( E1 - E2  )
        end if
      end if
      
      ! Calcul de la derivee de la partie irreversible
      ! Champ croissant
      if ( H > prevH ) then
        call dEverettEllipseda(H,prevMin(n),bh, dE1 )
        dMirrdH =  2.0d0*dE1
        
      ! Champ decroissant
      else if ( H < prevH ) then
        call dEverettEllipsedb(prevMax(n),H,bh, dE1 )
        dMirrdH = -2.0d0*dE1
         
      ! Champ constant
      else
        call dEverettEllipseda(H,prevMin(n),bh, dE1 )
        call dEverettEllipsedb(prevMax(n),H,bh, dE2 )
        dMirrdH = dE1 - dE2
      
      end if
      
      ! Normalisation
      !Mirr = Mirr / brem
      !dMirrdH = dMirrdH / brem
      
      Mirr = Mirr
      dMirrdH = dMirrdH
      ! Mise-a-jour de la valeur precedente de G
      
      
    end if
    
      
    ! Addition des composantes reversibles et irreversibles
    ! Induction
    !B = B  + Mirr
    B = Mirr
      
    ! Permeabilite differentielle
    !dBdH = dBdH + dMirrdH
    dBdH =  dMirrdH

    !write(*,*) B, dBdH, brem
    
    return
  end subroutine CalculPreisachBHEllipse
! -------------------------------------------------------------------------------------------------------------------


  subroutine EverettEllipse(alpha,beta,bh ,EverettEllipse_out)
    ! Description :
    ! Calcul de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Outputs :
    double precision, intent(out)                           :: EverettEllipse_out
    double precision                                        :: EverettEllipse1, EverettEllipse2
    ! Variables:
    ! Induction a remanence
    double precision                                        :: brem,Bdown1,Bdown2
    double precision                                        :: A2, C, D, discr,a,b,angle, C2, D2, discr2
    double precision                                        :: beta2, alpha2


    if (alpha < -beta) then
      beta2 = -alpha
      alpha2 = -beta
    else
      alpha2 = alpha
      beta2 = beta
    end if
   
    call GetEllipseParameters (bh,alpha2 , a,b,angle,brem)
    if (abs(alpha2) < epsilon(1.0d0) ) then
      Bdown1 = 0.0d0
      Bdown2 = 0.0d0
    else
      A2 = (sin(angle)/a)**2.0d0 + (cos(angle)/b)**2.0d0
      C = 2.0d0*cos(angle)*sin(angle)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(alpha2)
      D = ((cos(angle)/a)**2.0d0 + (sin(angle)/b)**2.0d0)*(alpha2**2.0d0) - 1.0d0
      discr = C**2.0d0 - 4.0d0*A2*D
      Bdown1 = (-C + sqrt(discr))/(2.0d0*A2)

      C2 = 2.0d0*cos(angle)*sin(angle)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(beta2)
      D2 = ((cos(angle)/a)**2.0d0 + (sin(angle)/b)**2.0d0)*(beta2**2.0d0) - 1.0d0
      discr2 = C2**2.0d0 - 4.0d0*A2*D2
      Bdown2 = (-C2 + sqrt(discr2))/(2.0d0*A2)
        
    end if
    EverettEllipse1 = (Bdown1 - Bdown2 )/2.0d0


    !EverettEllipse2 = -4.979d-05*(beta2-alpha2)+4.351d-09*(beta2**2.0d0-alpha2**2.0d0)+1.392d-13*(beta2**3.0d0-alpha2**3.0d0)-&
    !6.408d-13*(beta2*alpha2**2.0d0-beta2**2.0d0 * alpha2)
    EverettEllipse2 = -5.823d-05*(beta2-alpha2)+6.94d-9*(beta2**2-alpha2**2)

    if (alpha > 1000.0d0) then
      !write(*,*) EverettEllipse1, EverettEllipse2
    end if


    EverettEllipse_out = EverettEllipse1
    return
  end subroutine EverettEllipse
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine dEverettEllipseda(alpha,beta,bh, dEverettEllipseda_out )
    ! Description :
    ! Calcul de la derivee de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Outputs :
    double precision, intent(out)                           :: dEverettEllipseda_out
    double precision                                        :: dEverettEllipseda1, dEverettEllipseda2 
    ! Variables:
    ! Induction a remanence
    double precision                                        :: dBdowndH1, dBdowndH2
    double precision                                        :: alpha2, beta2
    

    if (alpha < -beta) then
      beta2 = -alpha
      alpha2 = -beta

    else
      alpha2 = alpha
      beta2 = beta
    end if
    !dEverettEllipseda2 = -4.979d-05*(-1.0d0)+4.351d-09*(-2.0d0*alpha2)+1.392d-13*(-3.0d0*alpha2**2.0d0)-&
    !6.408d-13*(2.0d0*beta2*alpha2-beta2**2.0d0)

    dEverettEllipseda2 = -5.823d-05*(-1)+6.94d-9*(-alpha2*2)

    call CalculdBdowndHmax1(alpha2,alpha2,bh, dBdowndH1)
    call CalculdBdowndHmax1(alpha2,beta2,bh, dBdowndH2)
    dEverettEllipseda1 = (dBdowndH1 -dBdowndH2 )/2.0d0
    if (alpha > 1000.0d0) then
      !write(*,*) dEverettEllipseda1, dEverettEllipseda2
    end if

    dEverettEllipseda_out = dEverettEllipseda1
    
    return
  end subroutine dEverettEllipseda
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine dEverettEllipsedb(alpha,beta,bh, dEverettEllipsedb_out )
    ! Description :
    ! Calcul de la derivee de la fonction d'Everett bilineaire.
    implicit none
    ! Inputs :
    double precision, intent(in)                            :: alpha, beta
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Outputs :
    double precision, intent(out)                           :: dEverettEllipsedb_out
    double precision                                        :: dEverettEllipsedb1, dEverettEllipsedb2
    ! Variables:
    ! Induction a remanence
    double precision                                        :: dBdowndH1, dBdowndH2

    double precision                                        :: alpha2, beta2

    

    if (alpha < -beta) then
      beta2 = -alpha
      alpha2 = -beta
    else
      alpha2 = alpha
      beta2 = beta
    end if
    !dEverettEllipsedb2 = -4.979d-05*(1.0d0)+4.351d-09*(2*beta2)+1.392d-13*(3*beta2**2.0d0)-&
    !6.408d-13*(alpha2**2.0d0-2.0d0*beta2 * alpha2)
    dEverettEllipsedb2 = -5.823d-05+6.94d-9*(beta2*2)
  


    call CalculdBdowndH(alpha2,alpha2,bh, dBdowndH1)
    call CalculdBdowndH(alpha2,beta2,bh, dBdowndH2)
    dEverettEllipsedb1 = ( dBdowndH1-dBdowndH2 )/2.0d0
    if (alpha > 1000.0d0) then
      !write(*,*) dEverettEllipsedb1, dEverettEllipsedb2
    end if

    dEverettEllipsedb_out = dEverettEllipsedb1

    return
  end subroutine dEverettEllipsedb
  ! ------------------------------------------------------------------------------------------------------------------- 

  subroutine CalculdBdowndH (Hmax,H,bh, dBdowndH)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    double precision, intent(in)                            :: Hmax
    ! Outputs :
    ! Valeur de la fonction B down et sa derivee
    double precision, intent(out)                           :: dBdowndH
      
    ! Variables :
    double precision                                        :: a, b, alpha, Bdown, discr, brem
    double precision                                        :: num, denom, A2, C, D
    ! Recuperation des parametres
    

    call GetEllipseParameters (bh,Hmax , a,b,alpha,brem)
    if (abs(Hmax) < epsilon(1.0d0) ) then
      Bdown = 0.0d0
      dBdowndH = 0.0d0
    else

      A2 = (sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0
      C = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(H)
      D = ((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(H**2.0d0) - 1.0d0
      discr = C**2.0d0 - 4.0d0*A2*D
      Bdown = (-C + sqrt(discr))/(2.0d0*A2)

      num = ((H*cos(alpha)+Bdown*sin(alpha))*cos(alpha))/a**2.0d0 + ((H*sin(alpha)-Bdown*cos(alpha))*sin(alpha))/b**2.0d0
      denom = (-(H*cos(alpha)+Bdown*sin(alpha))*sin(alpha))/a**2.0d0 + ((H*sin(alpha)-Bdown*cos(alpha))*cos(alpha))/b**2.0d0
      dBdowndH = num/denom

    end if

 

    

    return
  end subroutine CalculdBdowndH

  subroutine CalculdBdowndHmax1 (Hmax,H,bh, dBdowndHmax)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    double precision, intent(in)                            :: Hmax
    ! Outputs :
    ! Valeur de la fonction B down et sa derivee
    double precision, intent(out)                           :: dBdowndHmax
      
    ! Variables :
    double precision                                        :: a, b, alpha, Bdown, discr, brem
    double precision                                        :: num, denom, A2, C, D
    ! Recuperation des parametres


    call GetEllipseParameters (bh,Hmax , a,b,alpha,brem)


    if (abs(Hmax) < epsilon(1.0d0) ) then
      Bdown = 0.0d0
      dBdowndHmax = 0.0d0
    else
      A2 = (sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0
      C = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(H)
      D = ((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(H**2.0d0) - 1.0d0
      discr = C**2.0d0 - 4.0d0*A2*D
      Bdown = (-C + sqrt(discr))/(2.0d0*A2)

      num = ((H*sin(alpha)-Bdown*cos(alpha))**2.0d0)/(Hmax*b**2.0d0) + ((H*cos(alpha)+Bdown*sin(alpha))**2.0d0)/(Hmax*a**2.0d0)
      denom = ((H*cos(alpha)+Bdown*sin(alpha))*sin(alpha))/(a**2.0d0) - ((H*sin(alpha)-Bdown*cos(alpha))*cos(alpha))/(b**2.0d0)
      dBdowndHmax = num/denom
    end if
 


    return
  end subroutine CalculdBdowndHmax1

  subroutine CalculdBdowndHmax2 (Hmax,bh, dBdowndHmax)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: Hmax
    ! Outputs :
    ! Valeur de la fonction B down et sa derivee
    double precision, intent(out)                           :: dBdowndHmax
      
    ! Variables :
    double precision                                        :: a, b, alpha, Bdown, discr, brem
    double precision                                        :: num, denom, A2, C, D
    ! Recuperation des parametres


    call GetEllipseParameters (bh,Hmax , a,b,alpha,brem)


    if (abs(Hmax) < epsilon(1.0d0) ) then
      Bdown = 0.0d0
      dBdowndHmax = 0.0d0
    else
      A2 = (sin(alpha)/a)**2.0d0 + (cos(alpha)/b)**2.0d0
      C = 2.0d0*cos(alpha)*sin(alpha)*(1.0d0/(a**2.0d0) - 1.0d0/(b**2.0d0))*(Hmax)
      D = ((cos(alpha)/a)**2.0d0 + (sin(alpha)/b)**2.0d0)*(Hmax**2.0d0) - 1.0d0
      discr = C**2.0d0 - 4.0d0*A2*D
      Bdown = (-C + sqrt(discr))/(2.0d0*A2)

      num = ((Hmax*sin(alpha)-Bdown*cos(alpha))**2.0d0)/(Hmax*b**2.0d0) + &
      ((Hmax*cos(alpha)+Bdown*sin(alpha))**2.0d0)/(Hmax*a**2.0d0)- &
      (Hmax*cos(alpha)+Bdown*sin(alpha))*cos(alpha)/a**2.0d0 - (Hmax*sin(alpha)-Bdown*cos(alpha))*sin(alpha)/b**2.0d0

      denom = ((Hmax*cos(alpha)+Bdown*sin(alpha))*sin(alpha))/(a**2.0d0) - &
      ((Hmax*sin(alpha)-Bdown*cos(alpha))*cos(alpha))/(b**2.0d0)
      dBdowndHmax = num/denom
    end if
 


    return
  end subroutine CalculdBdowndHmax2


  subroutine CalculFG4param (bh,H , brem,F,dFdH,G,dGdH)
    ! Description :
    ! Calcul des fonction F(H) et G(H) pour le modele de Preisach. Recuperation de l'induction a remanence.
    implicit none
    ! Inputs :
    ! Parametres de la courbe BH
    double precision, dimension(maxParam), intent(in)       :: bh
    ! Champ mag
    double precision, intent(in)                            :: H
    ! Outputs :
    ! Induction a remanence
    double precision, intent(out)                           :: brem
    ! Valeur de la fonction F et sa derivee
    double precision, intent(out)                           :: F, dFdH
    ! Valeur de la fonction G et sa derivee
    double precision, intent(out)                           :: G, dGdH
      
    ! Variables :
    integer                                                 :: i, natan
    double precision, dimension(maxAtan)                    :: ai, bi, ci
    double precision                                        :: muFmax, bsat
    double precision                                        :: xia, xid
    
    
    brem = bh(2)
    muFmax = bh(8)*mu0

    ! Cas H = 0
    if ( abs(H) < epsilon(1.0d0) ) then
      F = 0.0d0
      dFdH = muFmax
      G = 0.0d0
      dGdH = 0.0d0
      
    ! Cas H < 0
    else if ( H < 0.0d0 ) then
      
      F = -(bh(3)-brem)*(-H/bh(7))*(1.0d0+(-H/bh(7))**(bh(5)+1.0d0))**(-1.0d0/(bh(5)+1.0d0))
      dFdH = ((bh(3)-brem)*(1.0d0+(-H/bh(7))**(bh(5)+1.0d0))**(-1.0d0/(bh(5)+1.0d0)))/(-H*(-H/bh(7))**(bh(5))+bh(7))
      G = -brem+brem*(1.0d0+(-H/bh(6))**(bh(5)+2.0d0))**(-1.0d0)
      dGdH = (brem*(bh(5)+2.0d0)*(-H)*(-H/bh(6))**(bh(5)))/((bh(6)**2.0d0)*(1.0d0+(-H/bh(6))**(bh(5)+2.0d0))**2.0d0)

      
    ! Cas H > 0
    else
      
      F = (bh(3)-brem)*(H/bh(7))*(1.0d0+(H/bh(7))**(bh(5)+1.0d0))**(-1.0d0/(bh(5)+1.0d0))
      dFdH = ((bh(3)-brem)*(1.0d0+(H/bh(7))**(bh(5)+1.0d0))**(-1.0d0/(bh(5)+1.0d0)))/(H*(H/bh(7))**(bh(5))+bh(7))
      G = brem-brem*(1.0d0+(H/bh(6))**(bh(5)+2.0d0))**(-1.0d0)
      dGdH = (brem*(bh(5)+2.0d0)*H*(H/bh(6))**(bh(5)))/((bh(6)**2.0d0)*(1.0d0+(H/bh(6))**(bh(5)+2.0d0))**2.0d0)
        
    end if
      
    return
  end subroutine CalculFG4param
end module