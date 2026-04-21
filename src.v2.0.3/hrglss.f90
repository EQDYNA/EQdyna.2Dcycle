SUBROUTINE hrglss
	use globalvar
	implicit none
	!
	!### program to calculate hourglass resistence and to add to
	!	residual force vector for the four-node quadrilateral
	!	element and assemble into the global right-hand-side 
	!	vector.
	!  Notice about some assumptions: 1-point Gaussian quadrature
	!	only; infinitesimal deformation so that coordinates x() 
	!	are not updated, so that lots of variables can be 
	!	clculated once and save for time step loops to save time.
	!  B.D. 7/22/05
	!
	logical :: lquad
	integer (kind=4) :: nel,i,j,k,m
	real (kind = dp),dimension(ned) :: phid
	real (kind = dp),dimension(ned,nen) :: dl,vl,fhr
	!
	!*** loop over elements ***
	!
	!$omp parallel default(shared) private(nel,i,j,k,m,lquad,dl,vl,fhr,phid)
	!$omp do schedule (static)
	do nel=1,numel
		m = 1
		!...if element degenerate
		lquad = .true.
		if(ien(3,nel) == ien(4,nel)) then
		  lquad = .false.
		endif
		!...only hourglass control for quadrilateral
		if(lquad) then
			!...localize dl and account for Rayleigh damping
			do i=1,nen
				k=ien(i,nel)
				do j=1,ned
					dl(j,i) = d(j,k)
					vl(j,i) = v(j,k)
					dl(j,i) = dl(j,i) + rdampk(m) * vl(j,i)
				enddo
			enddo  
			!... calculate sum(phi*dl)
			do i=1,ned
				phid(i) = 0.0d0
				do j=1,nen
					phid(i) = phid(i) + phi(j,nel) * dl(i,j)
				enddo
			enddo
			!... calculate hourglass resistence
			do i=1,nen
				fhr(1,i) = phi(i,nel) * (ss(1,nel)*phid(1) + ss(3,nel)*phid(2))
				fhr(2,i) = phi(i,nel) * (ss(3,nel)*phid(1) + ss(2,nel)*phid(2))
			enddo
			!... assemble to global right-hand-side force vector
			do i=1,nen
				do j=1,ned
					k=lm(j,i,nel)
					if(k > 0) then
                                                !$OMP ATOMIC
						brhs(k) = brhs(k) - fhr(j,i)
					endif
				enddo
			enddo
		endif
    !
	enddo
	!$omp end do nowait
	!$omp end parallel
    !
end SUBROUTINE hrglss
       
  
  
