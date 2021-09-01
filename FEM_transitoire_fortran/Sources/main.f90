program main
  ! -------------------------------------------------------------------------------------------------------------------
  ! Programme principal pour la resolution du slab problem 1D par FEM 1D
  ! Utilise le modele de Preisach pour le calcul de l'hysteresis
  ! Mai 2018
  ! Maxime Tousignant
  ! maxime.tousignant@polymtl.ca
  ! -------------------------------------------------------------------------------------------------------------------
  ! Description :
  !
  ! Le programme effectue une resolution transitoire des equations de Maxwell en regime basse frequence avec la 
  ! formulation en H. La resolution s'arrete lorsque le regime permanent est atteint.
  ! 
  ! Trois types de materiau magnetiques sont implementes.
  ! 1 : lineaire permeabilite relative constante
  ! 2 : non lineaire modele arctangente
  ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
  !
  ! Le programme gere les I/O via des fichiers text genere par le programme matlab SlabProblem_IO.m.
  !
  ! Parameters.txt contient :
  ! T         : Duree d'un cycle (s)
  ! ntpc      : Nombre de pas de temps par cycle
  ! L         : Taille du domaine de calcul (m)
  ! ne        : Nombre d'elements
  ! rho       : Resistivite (Ohms/m)
  ! mattype   : Type de courbe BH
  ! bh        : Parametres de la courbe BH
  !
  ! Dirichelet.txt contient :
  ! H_0_t     : Le champ magnetique en surface pour une periode en regime permanent
  !
  ! Results.txt contient : 
  ! H_x_t     : Champ magnetique H en fonction de la position x du temps pour une periode
  ! J_x_t     : Densite de courant J en fonction de la position x du temps pour une periode
  ! B_x_t     : Induction magnetique B en fonction de la position x du temps pour une periode
  ! Pj_x_t    : Pertes Joule en fonction de la position x du temps pour une periode
  ! Ph_x_t    : Pertes par hysteresis en fonction de la position x du temps pour une periode
  !
  ! -------------------------------------------------------------------------------------------------------------------
  ! Utilisation :
  !
  ! 1)  Telecharger gfortran a l'adresse suivante :
  !     https://gcc.gnu.org/wiki/GFortranBinaries
  !
  ! 2)  Compiler avec la commande suivante :
  !     gfortran Error.f90 Element.f90 Material.f90 SlabProblem.f90 main.f90 -o SlabProblem
  !
  ! 3)  Utiliser le programme matlab SlabProblem_IO.m pour lancer la simulation avec les parametres choisis et
  !     visualiser les resultats.
  !
  ! -------------------------------------------------------------------------------------------------------------------

  ! Modules
  use MaterialMod
  use SlabProblemMod
  
  ! Parametres et variables
  implicit none
  
  ! Path jusqu'au repertoire IO
  character(len=maxStrLen)                                  :: path
  
  ! Discretisation
  ! Nombre d'elements
  integer                                                   :: ne
  ! Nombre de degrees de liberte
  integer                                                   :: ndof
  ! Nombre de pas de temps par periode
  integer                                                   :: ntpc
  ! Pas de discretisation en espace
  double precision                                          :: dx
  ! Pas de temps
  double precision                                          :: dt
  
  ! Proprietes du materiau
  ! Resistivite electrique
  double precision                                          :: rho
  ! Parametres de la courbe BH
  ! Type de materiau magnetique
  ! 1 : lineaire permeabilite relative constante
  ! 2 : non lineaire modele arctangente
  ! 3 : non lineaire hysteretique modele de Preisach a 3n coefficients
  double precision, dimension(maxParam)                     :: bh
  

  ! Quantites physiques
  ! Champ magnetique
  double precision, dimension(:,:), allocatable             :: H_x_t
  ! Induction magnetique
  double precision, dimension(:,:), allocatable             :: B_x_t
  ! Densite de courant
  double precision, dimension(:,:), allocatable             :: J_x_t
  ! Pertes par courants de Foucault
  double precision, dimension(:,:), allocatable             :: Pj_x_t
  ! Pertes par hysteresis
  double precision, dimension(:,:), allocatable             :: Ph_x_t
  
  ! Erreur
  integer                                                   :: err
  
  ! -------------------------------------------------------------------------------------------------------------------
  
  ! Initialisation
  err = 0
  
  ! Recuperation du path
  call get_command_argument(1 , path)
  !path = "C:/Users/Gregory.Giard/Maitrise/Projet/mu_1D/FEM_transitoire_fortran/SlabProblem/"
  write(*,*) path
  ! Lecture des parametres a partir d'un fichier texte et verification
  call ReadParameters (path , ne,ndof,ntpc,dx,dt,rho,bh,err)

  ! Initialisation
  if (err == 0 ) call InitializeSlabProblem (path,ne,ndof,ntpc,dx,dt,rho, bh , H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t,err)
  
  ! Lancement de la simulation
  if (err == 0 ) call CalculSlabProblem (ne,ndof,ntpc,dx,dt,rho,bh , H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t , err)
  
  ! Ecriture des resultats
  if (err == 0 ) call WriteResults (path,ndof,ntpc,dx,dt,H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t , err)
  
  ! Desalocation des structures de donnees
  call FinalizeSlabProblem ( H_x_t,B_x_t,J_x_t,Pj_x_t,Ph_x_t )

  stop
end