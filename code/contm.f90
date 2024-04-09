SUBROUTINE contm(shg,det,elmass,constm)
	use globalvar
	implicit none     
	!=============================================================
	!### program to form mass matrix for a continuum element
	!        with "nen" nodes: lumped mass only here (imass=1)!
	!=============================================================
	integer (kind=4) :: l,j,n,k
	real (kind = dp) :: dsum,totmas,constm,temp1,temp2
	real (kind = dp),dimension(nint) :: det
	real (kind = dp),dimension(nee) :: elmass
	real (kind = dp),dimension(nen) :: work
	real (kind = dp),dimension(nrowsh,nen,nint) :: shg
	!...initialize
	dsum   = 0.0d0
	totmas = 0.0d0
	work   = 0.0d0
	!...calculate
	do l=1,nint
		temp1 = constm*w(l)*det(l)
		totmas = totmas + temp1
		do j=1,nen
			temp2 = temp1*shg(nrowsh,j,l)**2
			dsum = dsum + temp2
			work(j) = work(j) + temp2
		enddo
	enddo
	!...scale diagonal to conserve total mass
	temp1 = totmas/dsum
	!...store terms in a column
	do j=1,nen
		temp2 = temp1*work(j)
		n = (j - 1)*ned
		do k=1,ned
			elmass(n + k) = elmass(n + k) + temp2
			!it seems accumulation here does not make sense.
			!nevertherless, as elmass is intialized to be zero,
			!that does not matter in results. 
			!it seems two components(k=1,2) are same for one node.
			!B.D. 2/23/05
			!Accumulation makes sense if nint > 1. B.D. 7/26/05
		enddo
	enddo
	!
end SUBROUTINE contm
