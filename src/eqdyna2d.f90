PROGRAM eqdyna2d
	use globalvar
	implicit none

	real (kind = dp) :: timebegin, timeover
	integer (kind=4) :: ic, n,i,j,k,m,alloc_err
	character (len = 30) :: mm
	write(*,*) '====================================================================='
	write(*,*) '================== Welcome to EQdyna 2D ' // EQDYNA_VERSION // ' =================='
	write(*,*) '= UT Austin (Inst. for Geophysics) & Texas A&M (Center for Tectono.) ='
	write(*,*) '=========== Contacts: dliu@ig.utexas.edu, bduan@tamu.edu ============'
	write(*,*) '=                                                                   ='
	write(*,*) '=   EQdyna 2D uses FEM to simulate multicycle earthquake dynamics   ='
	write(*,*) '=   dynamic ruptures on geometrically realistic fault systems.      ='
	write(*,*) '=                                                                   ='
	write(*,*) '=   Model and system related parameters can be adjusted in          ='
	write(*,*) '=       FE_Global.txt,                                              ='
	write(*,*) '=       FE_Model_Geometry.txt, and                                  ='
	write(*,*) '=       FE_Fault_Geometry.txt.                                      ='
	write(*,*) '=                                                                   ='
	write(*,*) '=   Additional files required for earthquake cycle are              ='
	write(*,*) '=       Rate_direction.txt   (C_mesh=2 loading), or                 ='
	write(*,*) '=       nsmpGeoPhys.txt      (C_mesh=3 loading), and                ='
	write(*,*) '=       binaryop if this is not the first earthquake cycle          ='
	write(*,*) '====================================================================='
	
	call readglobal
	call readmodelgeometry		
	allocate(ftcn(ntotft), nfnode(ntotft))
	ftcn   = 0
	nfnode = 0
	call readfaultgeometry
	write(*,*) '=                                                                   ='
	write(*,*) '=     3 input files has been read in                                ='	
	
	amu     =    vs**2*rou
	lambda  =    vp**2*rou - 2.0d0*amu
	youngs  =    amu*(3.0d0*lambda + 2.0d0*amu)/(lambda + amu)
	poisr   =    lambda/2.0d0/(lambda + amu)
	
	!if (debug==1) write(*,*) 'ftcn'
	!if (debug==1) write(*,*) (ftcn(i),i=1,ntotft)
	ftn = maxval(ftcn)
	!if (debug==1) write(*,*) 'ftn = ', ftn
	
	xnode0  =    0.0d0
	ien0    =    0
	nsmp0   =    0
	nsmpnv  =    0.0d0
	write(*,*) '=                                                                   ='
	write(*,*) '=     Building finite element mesh ...                              ='	
	!if (debug==1) write(*,*) 'before meshgen'
	if      (C_mesh == 1) then
		call meshgen
	elseif 	(C_mesh == 2) then
		call meshgen1
	elseif  (C_mesh == 3) then
		call loadGmshMesh
	else
		write(*,*) '= ERROR: Unknown C_mesh value:', C_mesh
		stop
	endif
	!if (debug==1) write(*,*) 'after meshgen'
	
	totftnode = sum(nfnode)
	maxftnode = maxval(nfnode)
	
	allocate(output4plot(5,totftnode), fistr(2,totftnode), x(nsd,numnp))

	if (C_mesh == 2) then
		allocate(rd(2,maxftnode*ntotft))
		open(3,file='Rate_direction.txt',form='formatted',status='old')
			read(3,*) (rd(1,i),rd(2,i),i=1,maxftnode*ntotft)
		close(3)
		write(*,*) '=                                                                   ='
		write(*,*) '=     Rate_direction.txt is loaded                                  ='
	else
		allocate(nsmpgp(9,maxftnode*ntotft))
		nsmpgp = 0.0d0
		open(3,file='nsmpGeoPhys.txt',form='formatted',status='old')
			read(3,*) ((nsmpgp(j,i),j=1,9),i=1,maxftnode*ntotft)
		close(3)
		write(*,*) '=                                                                   ='
		write(*,*) '=     nsmpGeoPhys.txt is loaded                                     ='
		! Populate nsmpnv for faulting.f90 (dynamic rupture force calculation)
		! faulting.f90 indexes nsmpnv with i=(ift-1)*maxftnode+jj (padded), same as nsmpgp
		! nsmpnv(1)=-ty, nsmpnv(2)=tx (fault normal), nsmpnv(3)=segment length (m)
		do j = 1, ntotft
			do i = 1, nfnode(j)
				k = (j-1)*maxftnode + i
				nsmpnv(1,k) = -nsmpgp(2,k)          ! -ty -> fault normal x
				nsmpnv(2,k) =  nsmpgp(1,k)          !  tx -> fault normal y
				nsmpnv(3,k) =  nsmpgp(3,k)*1.0d3    ! segment length km->m
			enddo
		enddo
		write(*,*) '=     nsmpnv populated from nsmpGeoPhys.txt                         ='
	endif

	! Write key input arrays for cross-comparison debugging
	open(98, file='debug_nsmp0.txt', form='formatted', status='unknown')
	open(99, file='debug_nsmpnv.txt', form='formatted', status='unknown')
	if (C_mesh == 2) then
		open(96, file='debug_rd.txt', form='formatted', status='unknown')
	elseif (C_mesh == 3) then
		open(97, file='debug_nsmpgp.txt', form='formatted', status='unknown')
	endif
	do j = 1, ntotft
		do i = 1, nfnode(j)
			k = (j-1)*maxftnode + i
			write(98,'(3i10)') k, nsmp0(1,k), nsmp0(2,k)
			write(99,'(i10,3e18.8)') k, nsmpnv(1,k), nsmpnv(2,k), nsmpnv(3,k)
			if (C_mesh == 2) then
				write(96,'(i10,2e18.8)') k, rd(1,k), rd(2,k)
			elseif (C_mesh == 3) then
				write(97,'(i10,9e18.8)') k, (nsmpgp(n,k), n = 1, 9)
			endif
		enddo
	enddo
	close(98)
	close(99)
	if (C_mesh == 2) close(96)
	if (C_mesh == 3) close(97)
	if (C_mesh == 2) then
		write(*,*) '=     debug_nsmp0.txt, debug_nsmpnv.txt, and debug_rd.txt written    ='
	elseif (C_mesh == 3) then
		write(*,*) '=     debug_nsmp0.txt, debug_nsmpnv.txt, and debug_nsmpgp.txt written ='
	else
		write(*,*) '=     debug_nsmp0.txt and debug_nsmpnv.txt written                   ='
	endif

	do i = 1,nsd
		do j = 1,numnp
			x(i,j) = xnode0(i,j)
		enddo 
	enddo 
	
	allocate(ien(nen,numel), mat(numel),lm(ned,nen,numel), stat=alloc_err)
	if(alloc_err /= 0) then
		write(*,*) 'Insufficient space to allocate array ien'
		stop
	endif	
	
	do i = 1,nen
		do j = 1, numel
			ien(i,j) = ien0(i,j)
		enddo 
	enddo 
	
	nint = 1
	allocate(shl(nrowsh,nen,nint),w(nint))
	!write(*,*) 'Before qdcshl'
	call qdcshl
	!write(*,*) 'After qdcshl'
	
	allocate(eleb(nrowb,nee,numel),eledet(numel),elemass(nee,numel), &
  	   ss(3,numel),phi(nen,numel))	
	allocate(rho(numat),rdampm(numat),rdampk(numat),th(numat), &
  	    c(nrowc,nrowc,numat),e(numat), pois(numat))
		
	rho(1) = rou
	rdampm(1) = 0.0d0 
	rdampk(1) = 0.1d0
	e(1) = youngs
	pois(1) = poisr
	th(1) = 1.0d0
	nstep1 = floor(term/dt1) + 1
	!write(*,*) 'Total time steps =', nstep1
	!write(*,*) 'Before prop2d'
	call prop2d
	!write(*,*) 'After prop2d'
	
	allocate(fnms(numnp))
	fnms = 0.0d0
	
	allocate(id(ndof,numnp),stat=alloc_err)
	if(alloc_err /= 0) then
		write(*,*) 'Insufficient space to allocate array id'
		stop
	endif
	neq = 0	!establish equation numbers after above input
	do n=1,numnp
		do i=1,ndof
			neq = neq + 1
			id(i,n) = neq	!overwrite id array with equation number
		enddo
	enddo 	
	
	call formlm	 
	
	write(*,*) '=                                                                   ='
	write(*,*) '=     Model Information                                             ='
	write(*,*) '=                                                                   ='	
	write(*,*) '=     Total node number                                             ='	
	write(*,'(X,A,40X,i8)') '=',numnp
	write(*,*) '=     Total element number                                          ='	
	write(*,'(X,A,40X,i8)') '=',numel
	write(*,*) '=     Total equation number                                         ='
	write(*,'(X,A,40X,i8)') '=',neq	
	write(*,*) '=     Model x & y range                                             ='	
	write(*,'(X,A,40X,f7.2,3X,f7.2,2X,2A)') '=',xmax/1.0d3, xmin/1.0d3,'km'
	write(*,'(X,A,40X,f7.2,3X,f7.2,2X,2A)') '=',ymax/1.0d3, ymin/1.0d3,'km'

	allocate(alhs(neq),brhs(neq))
	alhs = 0.0d0	
	brhs = 0.0d0

	allocate(d(ndof,numnp),v(ndof,numnp))
	d = 0.0d0  
	v = 0.0d0	

	timeused(1) = (time2 - time1) * 1.d-9	!time in second
	do i = 1, maxftnode*ntotft
		fric(1,i) = fric_fs
		fric(2,i) = fric_fd
		fric(3,i) = fric_fv
		fnft(i) = 1000.0d0
	enddo
	do i = 1,ntotft
		fric(1,(i-1)*maxftnode+1) = 1000.0d0
		fric(1,(i-1)*maxftnode+nfnode(i)) = 1000.0d0
	enddo 
	call qdct2
	!write(*,*) 'After qdct2'
	!*** precompute reciprocal of lumped effective mass (alhs never changes after qdct2) ***
	!    enables the driver hot loop to use multiply instead of divide on every timestep.
	allocate(alhs_inv(neq))
	do i = 1, neq
		if (alhs(i) /= 0.0d0) then
			alhs_inv(i) = 1.0d0 / alhs(i)
		else
			alhs_inv(i) = 0.0d0
		endif
	enddo
	
	write(mm,'(i6)') icstart
	mm=trim(adjustl(mm))	
	
	open(4,file='cyclelog.txt'//mm,form = 'formatted', status = 'unknown')
		write(4,*) icstart, icstart
	close(4)	
	if (icstart>1) then 
		open(5, file = 'binaryop', form = 'unformatted', status = 'unknown')
			read(5) ((output4plot(i,j), i = 1,5), j = 1, totftnode)
		close(5)	
		do i = 1,totftnode
			fistr(1,i) = output4plot(1,i)
			fistr(2,i) = output4plot(2,i)
		enddo 				
	endif
	
	write(*,*) '=                                                                   ='
	write(*,*) '=     Readt to simulate earthquake cycles ...                       ='		

	do ic = icstart, icend
		write(*,*) '=                                                                   ='
		write(*,*) '=     Entering earthquake cycle No.                                 ='
		write(*,'(X,A,40X,i4)') '=',ic
		write(*,*) '=     Calculating interseismic deformation ...                      ='	
		
		call interstress(ic)
		open(94, file='debug_fistr_init.txt', form='formatted', position='append', status='unknown')
		n = 0
		do j = 1, ntotft
			do i = 1, nfnode(j)
				k = (j-1)*maxftnode + i
				n = n + 1
				write(94,'(3i10,2e18.8)') ic, k, n, fistr(1,n), fistr(2,n)
			enddo
		enddo
		close(94)
		
		!write(*,*) 'After interstress and before driver'
		
		fnft = 1000.0d0
		v = 0.0d0
		d = 0.0d0 
		write(*,*) '=                                                                   ='	
		write(*,*) '=     Calculating dynamic rupture ...                               ='	
		
		call driver

		do i = 1,totftnode
			fistr(1,i) = output4plot(1,i)
			fistr(2,i) = output4plot(2,i)
		enddo
		open(95, file='debug_output4plot.txt', form='formatted', position='append', status='unknown')
		n = 0
		do j = 1, ntotft
			do i = 1, nfnode(j)
				k = (j-1)*maxftnode + i
				n = n + 1
				write(95,'(3i10,5e18.8)') ic, k, n, (output4plot(m,n), m = 1, 5)
			enddo
		enddo
		close(95)
		open(5, file = 'binaryop', form = 'unformatted', status = 'unknown')
			write(5) ((output4plot(i,j), i = 1,5), j = 1, totftnode)
		close(5)
		open(61, file = 'totalop.txt'//mm, position = 'append', status = 'unknown')
			write(61,'(5e18.7)') ((output4plot(i,j), i = 1,5), j = 1, totftnode)
		close(61)
		open(4,file='cyclelog.txt'//mm,form = 'formatted', status = 'unknown')
			write(4,*) icstart, ic
		close(4)

		write(*,*) '=                                                                   ='
		write(*,*) '=     Finishing the current earthquake cycle                        ='
	enddo 	
	

end PROGRAM eqdyna2d
