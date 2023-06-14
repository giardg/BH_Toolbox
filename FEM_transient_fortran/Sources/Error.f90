module ErrorMod
  ! -------------------------------------------------------------------------------------------------------------------
  ! Mai 2018
  ! Maxime Tousignant
  ! maxime.tousignant@polymtl.ca
  ! -------------------------------------------------------------------------------------------------------------------
  ! Description :
  ! Gere l'ecriture des erreurs.
  ! -------------------------------------------------------------------------------------------------------------------
  
  ! Longueur maximale du message d'erreur
  integer, parameter                                        :: maxStrLen  = 256
  
  contains
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine Error (from,message )
    ! Description :
    ! Affichage des messages d'erreur.
    implicit none
    ! Inputs
    character(len=maxStrLen), intent(in)                    :: from
    character(len=maxStrLen), intent(in)                    :: message
  
    write(*,"(a,a,a)") "Error in ",trim(from)," :"
    write(*,*) trim(message)
  
    return
  end subroutine Error
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine Warning (from,message )
    ! Description :
    ! Affichage des messages d'erreur.
    implicit none
    ! Inputs
    character(len=maxStrLen), intent(in)                    :: from
    character(len=maxStrLen), intent(in)                    :: message
  
    write(*,"(a,a,a)") "Warning in ",trim(from)," :"
    write(*,*) trim(message)
  
    return
  end subroutine Warning
  ! -------------------------------------------------------------------------------------------------------------------
end module