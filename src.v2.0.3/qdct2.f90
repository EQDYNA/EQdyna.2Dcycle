SUBROUTINE qdct2
	use globalvar
	implicit none
	!
	!### program to calculate effective explicit,lumped mass matrix
	!	for 4-node quadrilateral,elastic continuum element
	!	and assemble into the global ledt-hand-side matrix
	!	B.D. 6/30/05
	!   Add in calculations of right-hand-side mass matrix and 
	!	element "b" matrix here, which were orignially in
	!	"qdct3.f90", to avoid recalculation or storage of 
	!	element shape function and avoid calculation of 
	!	time-independent varaibels in "qdct3.f90" for 
	!	speeding up. B.D. 7/3/05
	!   In addtion, right-hand-side mass is aslo needed to 
	!	assembled for "faulting" to use. B.D. 7/3/05
	!   Further, calculate and store arrays for hourglass control
	!	routine "hrglss" to use. B.D. 7/22/05
	!
	logical :: lquad
	integer (kind=4) :: nel,i,i1,j,m,k,k1
	real (kind = dp) :: constm,vol,ce,co
	real (kind = dp),dimension(nint) :: det
	real (kind = dp),dimension(nee) :: eleffm 	!lumped mass only 
	real (kind = dp),dimension(nesd,nen) :: xl
	real (kind = dp),dimension(2,2) :: xs
	real (kind = dp),dimension(nrowb,nee) :: b
	real (kind = dp),dimension(nrowsh,nen,nint) :: shg
  
  
	do nel=1,numel
		eleffm = 0.0d0	!initialize element mass array
		do i=1,nen		!localize coordinate array for the element
			i1 = ien(i,nel)
			do j=1,nesd
				xl(j,i) = x(j,i1)
			enddo
		enddo
		m = 1!mat(nel)	!material number
		lquad = .true.
		if (ien(3,nel) == ien(4,nel)) then
			lquad = .false.	!quadrilateral degenerate to triangular
		endif
		!
		!*** compute global shape function and derivatives ***
		!		at integration points
		call qdcshg(xl,det,shg,nel,lquad,xs)
		!
		!*** form mass matrix for left-hand-side***
		!
		constm = (1.0d0 + rdampm(m)*0.5d0*dt1)*th(m)*rho(m)
		if (constm /= 0.0d0) then
			call contm(shg,det,eleffm,constm)
		endif
		!
		!*** assemble element effective mass matrix into global ***
		!        left-hand-side matrix
		do i=1,nen
			do j=1,ned
				k=lm(j,i,nel)
				k1=j+(i-1)*ned
				if(k > 0) then
					alhs(k) = alhs(k) + eleffm(k1)
				endif
			enddo
		enddo
		!
		!*** following part to calculate and store global ***
		!	variables for time loop routines "qdct3.f90",
		!	"faulting.f90", and "hrglss.f90" to use for 
		!	computational efficiency. 
		! Note: this bases on an assumption that nodal coordinate
		!	does not need to update. This is a valid approximation
		!	for infinitesmal deformation. B.D. 7/22/05 
		!...... save determinant of Jacobi matrix for 1-point 
		!		Gaussian point only. B.D. 7/3/05
		eledet(nel) = det(1)
		!...... form and store mass matrix for right-hand-side ***
		!	note: different "constm" from above B.D. 7/3/05 
		constm = th(m)*rho(m)
		eleffm = 0.0d0	!initialize it again!!! B.D. 7/25/05
		call contm(shg,det,eleffm,constm)	!calculate
		do i=1,nee
			elemass(i,nel) = eleffm(i)	!store
		enddo
		!......assemble element mass to nodal mass 
		!	for "faulting" to use. B.D. 7/3/05
		!	note: ned degree of an element node have a same mass 
		!	  (see "contm"), so only use one of them to assemble.
		do i=1,nen
			k  = ien(i,nel)
			k1 = 1 + (i-1) *ned
			fnms(k) = fnms(k) + eleffm(k1)
		enddo
		!...... calculate strain-displacement matrix "b" for
		! 	1-Gaussian point case only and store.
		!	   note: nint=1 is assumed. B.D. 7/3/05
		
		call qdcb(shg,b)	!calculate
		
		do j=1,nee
			do i=1,nrowb
				eleb(i,j,nel) = b(i,j)	!store
			enddo
		enddo
		!...... calculate S matrix for hourglass control
		!  Note: isotropic, Posian material is assumed.
		!  	th(m) = 1 is explicitly used to save time.
		!  Note: entries of xs come from qdcshg routine.
		vol = 0.5d0 * ((xl(1,3)-xl(1,1)) * (xl(2,4)-xl(2,2))  &
				   + (xl(1,2)-xl(1,4)) * (xl(2,3)-xl(2,1)))
		!ce = 16.0d0* e(m) / 15.0d0
		ce = 4.0d0*amu*(lambda+amu)/(lambda+2.0d0*amu)
		co = ce * vol / 48.d0
		ss(1,nel) = co * (xs(2,2)*xs(2,2) + xs(1,2)*xs(1,2))
		ss(2,nel) = co * (xs(2,1)*xs(2,1) + xs(1,1)*xs(1,1))
		ss(3,nel) = - co * (xs(2,2)*xs(2,1) + xs(1,2)*xs(1,1))
		!...... calculate phi for hourglass control
		!  Based on Hughes (2000) pp.252-254.
		!  Note: only 1-point gaussian works so shg(:,:,1).
		do i=1,nen
			vol = 0.0d0		!temporary use
			do j=1,nen
				vol = vol + (-1.0d0)**j * (xl(1,j) * shg(1,i,1)  &
					 + xl(2,j) * shg(2,i,1))
			enddo
			phi(i,nel) = (-1.0d0)**i - vol 	!middle value
		enddo
		vol = 0.0d0
		do i=1,nen
			vol = vol + phi(i,nel) * phi(i,nel)
		enddo
		vol = sqrt(vol)
		do i=1,nen
			phi(i,nel) = 2.0d0 * phi(i,nel) / vol	!normalize
		enddo
		!  
	enddo
	!
end SUBROUTINE qdct2
