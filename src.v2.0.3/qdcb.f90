SUBROUTINE qdcb(shg,b)
	use globalvar
	implicit none
	!
	!### program to set up the strain-displacement matrix "b" for
	!	2-D continuum elements.
	!	note: ibbar=0 (standard b-matrix) only
	!
	integer (kind=4) :: j,j2,j2m1
	real (kind = dp),dimension(nrowsh,nen) :: shg
	real (kind = dp),dimension(nrowb,nee) :: b
	!
	b = 0.0d0	!initialize
	!
	do j=1,nen
		!
		j2   = 2*j
		j2m1 = j2 - 1
		!
		b(1,j2m1) = shg(1,j)
		b(1,j2  ) = 0.0d0
		b(2,j2m1) = 0.0d0
		b(2,j2  ) = shg(2,j)
		b(3,j2m1) = shg(2,j)
		b(3,j2  ) = shg(1,j)
	enddo
  !
end SUBROUTINE qdcb
