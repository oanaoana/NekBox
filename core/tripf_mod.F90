!> cleaned
module tripf
  use kinds, only : DP
  use size_m
  implicit none
!     Arrays and parameters needed for tripf.f
      integer, parameter :: maxwalls = 12
      integer, parameter :: maxwallpar = 5, maxpwallpar = 1
      integer :: Nelx=16,Nely=16,Nelz=4
      real(DP) :: smth=0.2

      integer, parameter :: maxlxyz=max0(lx1,ly1,lz1)
 
      integer ::  seed,ntdt,trip
      Integer ::  nwalls,nnelx1x2(maxwallpar),kpts(maxwallpar), nwallpar,npwallpar
      real(DP) :: wallpar(maxwallpar), pwallpar(3,maxpwallpar,maxwalls)
      real(DP) :: znek(lelv*maxlxyz,maxwalls)
      real(DP) :: fzt1(lelv*maxlxyz,maxwalls)
      real(DP) :: fzt2(lelv*maxlxyz,maxwalls)
      real(DP) :: fzt3(lelv*maxlxyz,maxwalls)
      real(DP) :: fzt4(lelv*maxlxyz,maxwalls)
      real(DP) :: tripx(maxwalls),tripy(maxwalls),tripz(maxwalls)
      character(1):: direction(maxwalls)

!  contains

!  subroutine init_mesh()
!    use size_m, only : lelt, ldim
!    implicit none

!    ifsolv = .false.
!    allocate(vertex((2**ldim)*lelt))
!    allocate(ifdfrm(LELT), iffast(lelt))
!    iffast = .true.
    
!  end subroutine init_mesh


end module tripf
