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
	
    timeused(2) = timeused(2) + (time2 - time1) * 1.d-9   
    !
    quitdriver = .false.
    timedyna = 0.0d0	
    do n=1,nstep1
		!
		timedyna = timedyna + dt1	 	
		!...print on screen for monitoring
		if (mod(n,2000) == 1) then
			write(*,*) '=                                                                   ='
			write(*,*) '=     Current time in dynamic rupture                               ='
			write(*,'(X,A,40X,f7.3,4X,A)') '=',timedyna, 's'
			write(*,*) '=     Your earthquake is rampaging ... Be patient ...               ='
		endif

		do l=1,numnp
			do j=1,ndof
				k=id(j,l)
				if(k > 0) then	!only non-boundary,update
					v(j,l) = v(j,l) + brhs(k) * dt1
					d(j,l) = d(j,l) + v(j,l) * dt1
				endif
			enddo
		enddo
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
		
		!      time1 = gethrtime()
		call qdct3
		!if (debug == 1) write(*,*) 'AFter qdct3'
		!      time2 = gethrtime()
		timeused(3) = timeused(3) + (time2 - time1) * 1.d-9 
		!      time1 = gethrtime()
		call hrglss
		!if (debug == 1) write(*,*) 'AFter hrglss'
		!      time2 = gethrtime()
		timeused(4) = timeused(4) + (time2 - time1) * 1.d-9 
		!
		!*** call faulting subrotuine to revise residual force ***
		!  This is main revision on general dynamic code to
		!  study earthquke rupture problems. The main purpose
		!  is to revise right-hand-side vector "brhs" by 
		!  faulting boundary. B.D. 7/3/05
		!
		! if(n == 1) then !note: 1 time step advance of faulting results
			! lstrf = .true.
		! elseif (mod(n+1,nhplt1) == 0) then 	
			! lstrf = .true.	!to store faulting results.
		! else
			! lstrf = .false.
		! endif
		!       time1 = gethrtime()
		
		call faulting(n)
		if (debug == 1) write(*,*) 'AFter faulting'
		!      time2 = gethrtime()
		timeused(5) = timeused(5) + (time2 - time1) * 1.d-9 

		do i = 1,neq
			brhs(i) = brhs(i)/alhs(i)
		enddo
		!
		if (quitdriver == .TRUE.) then
			exit
		endif
	enddo 	!end time step loop n
	!
end SUBROUTINE driver
