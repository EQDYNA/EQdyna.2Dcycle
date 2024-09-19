program eqdyna2d
    use globalvar
    implicit none

    real (kind = dp) :: timebegin, timeover
    integer (kind=4) :: ic, n,i,j,k,alloc_err
    character (len = 30) :: mm, cycleId
    write(*,*) '====================================================================='
    write(*,*) '================== Welcome to EQdyna 2D 2.0.3 ======================='
    write(*,*) '===== Product of Earthquake Modeling Lab @ Texas A&M University ====='
    write(*,*) '========== Website https://seismotamu.wixsite.com/emlam ============='
    write(*,*) '=========== Contacts: dliu@ig.utexas.edu, bduan@tamu.edu ============='
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
    write(*,*) '=       Rate_direction.txt, and                                     ='
    write(*,*) '=       binaryop if this is not the first earthquake cylce          ='
    write(*,*) '====================================================================='
    
    call readglobal
    call readmodelgeometry        
    allocate(ftcn(ntotft), nfnode(ntotft))
    ftcn = 0
    nfnode = 0
    
    call readfaultgeometry
    write(*,*) '=                                                                   ='
    write(*,*) '=     3 input files has been read in                                ='    
    
    amu     =    vs**2*rou
    lambda  =    vp**2*rou - 2.0d0*amu
    youngs  =    amu*(3.0d0*lambda + 2.0d0*amu)/(lambda + amu)
    poisr   =    lambda/2.0d0/(lambda + amu)
    
    if (debug==1) write(*,*) 'ftcn'
    if (debug==1) write(*,*) (ftcn(i),i=1,ntotft)
    ftn = maxval(ftcn)
    if (debug==1) write(*,*) 'ftn = ', ftn
    
    xnode0  =    0.0d0
    ien0    =    0
    nsmp0   =    0
    nsmpnv  =    0.0d0
    write(*,*) '=                                                                   ='
    write(*,*) '=     Building finite element mesh ...                              ='    
    if (debug==1) write(*,*) 'before meshgen'
    if      (C_mesh == 1) then 
        call meshgen
    elseif     (C_mesh == 2) then 
        call meshgen1
    elseif (C_mesh==3) then
        call loadGmshMesh
    endif 
    if (debug==1) write(*,*) 'Mesh is generated/loaded.'
    
    totftnode = sum(nfnode)
    maxftnode = maxval(nfnode)
    if (debug==1) write(*,*) 'Total num of ft nodes is', totftnode, 'max num per ft node is', maxftnode
    
    allocate(output4plot(5,totftnode), fistr(2,totftnode), x(nsd,numnp), &
        rd(2,maxftnode*ntotft))
        
    !open(3,file='Rate_direction.txt',form = 'formatted',status = 'old')
    !        read(3,*) (rd(1,i),rd(2,i),i=1,maxftnode*ntotft)
    !close(3)
    do i = 1, maxftnode*ntotft
        rd(1,i) = 300.d0
        rd(2,i) = 10.d0
    enddo 
    write(*,*) '=                                                                   ='
    write(*,*) '=     Rate_direction.txt is loaded                                  ='    
    
    do i = 1,nsd
        do j = 1,numnp
            x(i,j) = xnode0(i,j)
        enddo 
    enddo 
    
    allocate(ien(nen,numel), mat(numel),lm(ned,nen,numel), stat=alloc_err)
    if(alloc_err /= 0) then
        write(*,*) 'Insufficient space to allocate array ien'
         
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
         
    endif
    neq = 0    !establish equation numbers after above input
    do n=1,numnp
        do i=1,ndof
            neq = neq + 1
            id(i,n) = neq    !overwrite id array with equation number
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

    timeused(1) = (time2 - time1) * 1.d-9    !time in second
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
    
    write(mm,'(i6)') icstart
    mm=trim(adjustl(mm))    
    
    open(4, file='cyclelog.txt'//mm, form = 'formatted', status = 'unknown')
        write(4,*) icstart, icstart
    close(4)    
    if (icstart>1) then 
        open(5, file = 'restart.txt', form='formatted', status = 'unknown')
            read(5,*) ((output4plot(i,j), i = 1,5), j = 1, totftnode)
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
        
        !write(*,*) 'After interstress and before driver'
        
        fnft = 1000.0d0
        v = 0.0d0
        d = 0.0d0 
        write(*,*) '=                                                                   ='    
        write(*,*) '=     Calculating dynamic rupture ...                               ='    
        
        call driver    
        
        write(*,*) 'Exiting driver. Write results to files. '
        do i = 1,totftnode
            fistr(1,i) = output4plot(1,i)
            fistr(2,i) = output4plot(2,i)
        enddo         
        
        ! open(5, file = 'binaryop', form = 'unformatted', status = 'unknown')
            ! write(5) ((output4plot(i,j), i = 1,5), j = 1, totftnode)
        ! close(5)
        
        open(5, file = 'restart.txt', form = 'formatted', status = 'unknown')
            write(5,*) ((output4plot(i,j), i = 1,5), j = 1,totftnode)
        close(5)
        
        write(cycleId,'(i6)') ic
        cycleId=trim(adjustl(cycleId))  
    
        open(1001, file = 'totalop.txt'//cycleId, form='formatted', status = 'unknown')
            write(1001,'(5e18.7)') ((output4plot(i,j), i = 1,5), j = 1, totftnode)
        close(1001)   
        
        ! open(6, file = 'totalop.txt'//mm, position = 'append', status = 'unknown')
            ! write(6,'(5e18.7)') ((output4plot(i,j), i = 1,5), j = 1, totftnode)
        ! close(6)    
        open(4,file='cyclelog.txt'//mm,form = 'formatted', status = 'unknown')
            write(4,*) icstart, ic
        close(4)    
        
        write(*,*) '=                                                                   ='    
        write(*,*) '=     Finishing the current earthquake cycle                        ='
    enddo     
    

end program eqdyna2d
