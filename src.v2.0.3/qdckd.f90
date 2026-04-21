SUBROUTINE qdckd(bb,cc,dl,elresf,constk)
  use globalvar
  implicit none
  !
  !### program to form internal force ("-k*d") for a continuum element
  !        with "nen" nodes
  !	Avoid calling subroutines for matrix operations to speed up.
  !	B.D. 3/24/05
  !	Further speed up by taking into account zero in B and D 
  !		for 2-D, isotropic, plane-strain(stress) elasticity.
  !	B.D. 3/24/05 again.
  !	Rewrite again: strain-displacement matrix "b" was stored before
  !		time loops. B.D. 7/2/05
  !
  integer (kind=4) :: l,j,i
  real (kind = dp) :: temp,constk
  real (kind = dp),dimension(nee) :: elresf,work,dl
  real (kind = dp),dimension(nstr) :: strain,stress
  real (kind = dp),dimension(nrowb,nee) :: bb	!correspond to b
  real (kind = dp),dimension(nrowc,nrowc) :: cc	!correspond to c 
  !
  !*** loop on integration points ***
  !	at present, nint=1 is assumed. B.D. 7/2/05
  do l=1,nint
    temp = constk*w(l)
    !... calculate strains
    strain(1) = bb(1,1)*dl(1) + bb(1,3)*dl(3) + bb(1,5)*dl(5) + bb(1,7)*dl(7)
    strain(2) = bb(2,2)*dl(2) + bb(2,4)*dl(4) + bb(2,6)*dl(6) + bb(2,8)*dl(8)
    strain(3) = 0.0d0
    do j=1,nee
      strain(3) = strain(3) + bb(3,j)*dl(j)
    enddo
    !... calculate stresses
    stress(1) = cc(1,1)*strain(1) + cc(1,2)*strain(2)
    stress(2) = cc(2,1)*strain(1) + cc(2,2)*strain(2)
    stress(3) = cc(3,3)*strain(3)
    !... calculate element internal force
    do i=1,nstr
      stress(i) = temp * stress(i)
    enddo
    work(1) = bb(1,1)*stress(1) + bb(3,1)*stress(3)
    work(2) = bb(2,2)*stress(2) + bb(3,2)*stress(3)
    work(3) = bb(1,3)*stress(1) + bb(3,3)*stress(3)
    work(4) = bb(2,4)*stress(2) + bb(3,4)*stress(3)
    work(5) = bb(1,5)*stress(1) + bb(3,5)*stress(3)
    work(6) = bb(2,6)*stress(2) + bb(3,6)*stress(3)
    work(7) = bb(1,7)*stress(1) + bb(3,7)*stress(3)
    work(8) = bb(2,8)*stress(2) + bb(3,8)*stress(3)
    do i=1,nee
      elresf(i) = elresf(i) + work(i)
    enddo
  enddo
end SUBROUTINE qdckd
