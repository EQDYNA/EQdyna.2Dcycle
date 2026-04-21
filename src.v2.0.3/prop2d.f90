SUBROUTINE prop2d
  use globalvar
  implicit none
  !
  !### program to read and store properties for two-dimensional
  !        continuum elements
  !
  !        note: this routine is presently restricted to the
  !              isotropic linearly-elastic case
  !
  !              iopt = 0; plane stress
  !                   = 1; plane strain
  !		 see book pp.83
  !
  integer (kind=4) :: n,m
  real (kind = dp) :: amu2,alam
  do n=1,numat
    rdampk(n) = rdampk(n) * dt1
	
    !...set material constants for out-of-plane components
    amu2 = e(n)/(1.0d0 + pois(n))
    alam = amu2*pois(n)/(1.0d0 - 2.0d0 * pois(n))
    !
    c(1,4,n) = alam
    c(2,4,n) = alam
    c(3,4,n) = 0.0d0
    c(4,4,n) = alam + amu2
    !
    c(4,1,n) = c(1,4,n)
    c(4,2,n) = c(2,4,n)
    c(4,3,n) = c(3,4,n)
    !... set material constants for in-plane components
    if (iopt == 0) then
      alam = alam*amu2/(alam + amu2)
    endif
    !
    c(1,1,n) = alam + amu2
    c(1,2,n) = alam
    c(2,2,n) = c(1,1,n)
    c(1,3,n) = 0.0d0
    c(2,3,n) = 0.0d0
    c(3,3,n) = 0.5d0*amu2
    !
    c(2,1,n) = c(1,2,n)
    c(3,1,n) = c(1,3,n)
    c(3,2,n) = c(2,3,n)
  enddo

end SUBROUTINE prop2d
