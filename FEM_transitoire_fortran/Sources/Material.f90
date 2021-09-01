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
  private                                                   :: Everett
  private                                                   :: dEverettda
  private                                                   :: dEverettdb
  
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
      

    else
      err = 1
      message = "Le type de materiau choisi n''est pas implemente."
      call Error (from,message )
      return
      
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
    double precision                                        :: prevG, partialM
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
    
    ! Initialisation
    B = mu0*H
    dBdH = mu0
    Mirr = 0.0d0
    dMirrdG = 0.0d0
    
    ! Calcul des fonctions F(H) et G(H) et recuperation de brem
    call CalculFG (bh,H , brem,F,dFdH,G,dGdH)

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
end module MaterialMod
    