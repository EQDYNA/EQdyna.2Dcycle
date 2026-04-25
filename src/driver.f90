SUBROUTINE driver
	use globalvar
	implicit none
	!
	!### solution driver program
	! As shown above, only fault arrays are transferred.
	!
	! Modify to explictly use central difference method as DYNA3D,
	!	rather than alpha-method. B.D. 7/21/05
	!
	logical :: lstr,lstrf
	integer (kind=4) :: nsave,nsq,n,i,j,k,k1,l,nel,m, &
	   ndprt1,nsprt1,nhplt1,niter1,ifault
	integer (kind=8) :: tic, toc, clockrate
	integer :: dt_vals(8)
	character(len=23) :: ts
	!   print interval (in steps) for the in-loop status line.
	!   dt1=0.01 s, so 300 steps = every 3 simulated seconds.
	integer, parameter :: print_every = 300

    timeused(2) = timeused(2) + (time2 - time1) * 1.d-9
    !...reset per-cycle wall-clock accumulators at start of each dynamic-rupture run
    t_vd = 0.0d0; t_qdct3 = 0.0d0; t_hrglss = 0.0d0; t_faulting = 0.0d0; t_brhs = 0.0d0
    call system_clock(count_rate=clockrate)
    !
    quitdriver = .false.
    timedyna = 0.0d0
    do n=1,nstep1
		!
		timedyna = timedyna + dt1
		!...per-step progress line (one line per step by default; change print_every to thin out)
		if (mod(n,print_every) == 0) then
			call date_and_time(values=dt_vals)
			write(ts,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,".",I3.3)') &
				dt_vals(1),dt_vals(2),dt_vals(3),dt_vals(5),dt_vals(6),dt_vals(7),dt_vals(8)
			write(*,'(A,A,A,I6,A,F7.3,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A)') &
				'[',ts,']  n=',n,'  t=',timedyna,'s  vd=',t_vd, &
				' qdct3=',t_qdct3,' hrglss=',t_hrglss,' flt=',t_faulting,' brhs=',t_brhs,' s'
			flush(6)
		endif

		call system_clock(tic)
		!$omp parallel do default(shared) private(l,j,k) schedule(static)
		do l=1,numnp
			do j=1,ndof
				k=id(j,l)
				if(k > 0) then	!only non-boundary,update
					v(j,l) = v(j,l) + brhs(k) * dt1
					d(j,l) = d(j,l) + v(j,l) * dt1
				endif
			enddo
		enddo
		!$omp end parallel do
		call system_clock(toc); t_vd = t_vd + real(toc-tic, dp)/real(clockrate, dp)
	!if (timedyna < 1.0d0) then
	!	codetermination = .false.
	!else
	!	codetermination = .true.
	!	do l = 1, numnp
	!		do j = 1, ndof
	!			k = id(j,l)
	!			if (k>0) then
	!				if (abs(v(j,l)) > 0.02d0) then
	!					codetermination = .false.
	!				endif
	!			endif
	!		enddo
	!	enddo
	!endif
		! !*** store desired results at set time intervals ***
		! if(n == 1) then
			! lstr = .true.
			! nsave = 1
		! elseif (mod(n,nhplt1) == 0) then	
			! lstr = .true.	!to store faulting results.
			! nsave = nsave + 1
			! locplt = locplt + 1
		! else
			! lstr = .false.
		! endif
		! if(lstr) then
			! ftimestr(nsave) = timedyna	!save time info for output
		! endif
		! if (lstr) then
			! !...nodal output
			! if((ndout > 0) .and. (locplt > 1)) then
				! dout(1,locplt) = timedyna
				! do i=1,ndout
					! j = idhist(1,i)
					! k = idhist(2,i)
					! l = idhist(3,i)
					! if(l == 1) then
						! dout(i+1,locplt) = d(k,j)
					! elseif(l == 2) then
						! dout(i+1,locplt) = v(k,j)
					! elseif(l == 3) then
						! k1 = id(k,j)
						! dout(i+1,locplt) = brhs(k1)
					! endif
				! enddo
			! endif
			! !	
		! endif
		!
		!*** initialize for right hand force to use ***
		!
		brhs = 0.0d0	!initialize it every interation

		call system_clock(tic); call qdct3;     call system_clock(toc)
		t_qdct3    = t_qdct3    + real(toc-tic, dp)/real(clockrate, dp)

		call system_clock(tic); call hrglss;    call system_clock(toc)
		t_hrglss   = t_hrglss   + real(toc-tic, dp)/real(clockrate, dp)

		!*** call faulting to revise residual force ***
		call system_clock(tic); call faulting(n); call system_clock(toc)
		t_faulting = t_faulting + real(toc-tic, dp)/real(clockrate, dp)
		if (debug == 1) write(*,*) 'AFter faulting'

		call system_clock(tic)
		!$omp parallel do schedule(static)
		do i = 1, neq
			brhs(i) = brhs(i) * alhs_inv(i)
		enddo
		!$omp end parallel do
		call system_clock(toc); t_brhs = t_brhs + real(toc-tic, dp)/real(clockrate, dp)
		!
		if (quitdriver .eqv. .TRUE.) then
			exit
		endif
	enddo 	!end time step loop n
	!
end SUBROUTINE driver
