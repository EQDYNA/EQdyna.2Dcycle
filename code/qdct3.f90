SUBROUTINE qdct3
	use globalvar
	implicit none
	!
	!### program to calculate element contributions to residual force
	!        for the four-node quadrilateral, elastic continuum element 
	!        and assemble into the global right-hand-side vector.
	!		Note: only for 1-Gaussian point case now! B.D. 7/2/05
	!
	!  Explicitly use central difference method: al() always zero if no
	!	Rayleigh dampling (rdampm = 0), no global a() is needed.
	!	B.D. 7/21/05
	!
	logical :: formkd,zerodl,formma,zeroal
	integer (kind=4) :: nel,m,i,j,ntemp,k,k1
	real (kind = dp) :: det,constk
	real (kind = dp),dimension(nee) :: elresf,eleffm
	real (kind = dp),dimension(ned,nen) :: dl,vl,al
	!real (kind = dp),dimension(nrowsh,nen) :: shg

	!
	!*** loop over elements ***
	!
	!!$omp parallel do default(shared) private(formkd,zerodl,formma,zeroal,nel,j,i,k,m,ntemp,k1,dl,vl,al,det,constk,elresf,eleffm)
	do nel=1,numel
		!
		formma = .false.
		formkd = .false.
		m = 1
		!...localize dl,vl,al
		do j=1,nen
			ntemp = ien(j,nel)
			do i=1,ned
				dl(i,j) = d(i,ntemp)
				vl(i,j) = v(i,ntemp)
				al(i,j) = 0.0d0
			enddo
		enddo
		!...compute effective dl accounting for Rayleigh damping
		do j=1,nen
			do i=1,ned
				dl(i,j) = dl(i,j) + rdampk(m)*vl(i,j)
				al(i,j) = al(i,j) + rdampm(m)*vl(i,j)
			enddo
		enddo
		!...determine if element makes inertial contribution
		zeroal = .true.
		outer1: do j=1,nen	!Giving names to control constructs
			inner1: do i=1,ned
						k=id(i,j)
						if(al(i,j) /= 0.0d0) then
							zeroal = .false.
							exit outer1
						endif
					enddo inner1
				enddo outer1
		if ( (.not.zeroal) .and. (imass /= 2) .and. (rho(m) /= 0.0d0) ) then
			formma = .true.
		endif
		!...determine if element makes stiffness contribution
		zerodl = .true.
		outer2: do j=1,nen
			inner2: do i=1,ned
						if(dl(i,j) /= 0.0d0) then
							zerodl = .false.
							exit outer2
						endif
					enddo inner2
				enddo outer2
		if (.not.zerodl) then
			formkd = .true.
		endif
		!...gravity effect and determine formma again
		!    at present, no gravity. zerog = .true. B.D. 3/26/05
		!	  zerog = .true.
		!     call ztest(grav,nesd,zerog)
		!  if ((.not.zerog) .and. (lfbody.ne.0) .and. (rho(m).ne.zero)	&
		! &   .and. (imass.ne.2)) then
		!     formma = .true.
		!     do 400 i=1,ned
		!     temp = grav(i)*g1(lfbody)
		!         do 300 j=1,nen
		!           al(i,j) = al(i,j) - temp
		!300    continue
		!400  continue
		!  endif
		!*** if either, start time-consuming computing ***
		if (formma .or. formkd) then
			!...initialize element right hand vector
			elresf = 0.0d0
			!...directly assign shg and det to avoid recalculation to speed up.
			! 	only for 1-point Gaussian case. 	B.D. 3/25/05
			!	   No shg needed in this new formation. B.D. 7/3/05
			!do i=1,nrowsh
			!	do j=1,nen
			!	  shg(i,j) = eleshap(i,j,nel)
			!	enddo
			!enddo
			det = eledet(nel)
			!*** form inertial and/or body force ***
			if (formma) then
				!......assign eleffm from previous stored. B.D. 3/26/05
				do i=1,nee
					eleffm(i)=elemass(i,nel)
				enddo
				!......call to compute elresf	
				call contma(eleffm,al,elresf)
			endif
			!*** form internal force: most time-consuming part ***
			if (formkd) then
				!...... form internal force
				constk = - det * th(m)
				call qdckd(eleb(1,1,nel),c(1,1,m),dl,elresf,constk)
			endif
			!*** only either formma or formkd, assemble ***
			do i=1,nen
				do j=1,ned
					k=lm(j,i,nel)
					k1=j+(i-1)*ned
					if(k > 0) then
						brhs(k) = brhs(k) + elresf(k1)
					endif
				enddo
			enddo
		endif
		!
	enddo
	!!$omp end parallel do
	!
	!*** form surface force ***
	!     note: assembly of surface loads is performed inside qdcsuf
	!	at present, no surface force applied! B.D. 7/2/05
	!if ( (nsurf.gt.0) .and. (lfsurf.gt.0)) then
	!   call qdcsuf(ielno,ien,x,xl,iside,mat,th,press,shear,elresf, &
	!             brhs,lm,g1(lfsurf),nsurf,nen,nsd,nesd,ned,nee,iopt)
	!endif
	!
end SUBROUTINE qdct3
