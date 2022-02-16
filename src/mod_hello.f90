module mod_hello

use prog_type
implicit none
!
contains
!
!-----------------------------------------------------------------------
!
  subroutine hello_world()
    !
    ! ... arguments
    !
    write(*,*) 'Hello, this is the cont3Dslab startup program'
    !
    !
  end subroutine hello_world
  !
end module mod_hello
  
