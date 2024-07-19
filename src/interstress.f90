subroutine interstress(ic)

    use globalvar
    implicit none

    integer (kind=4) :: ic, ida(maxftnode*ntotft), ift, ntag, i, j, nuc, ndt
    real (kind = dp) :: shs(totftnode), ns(totftnode), ss0(totftnode), ns0(totftnode), &
        strengthexcess(totftnode) 
    real (kind = dp) :: eta, etaMin, dtInSeconds, dtInYears, tinter, tcon, rn, rs
    real (kind = dp) :: theta, minstrengthexcess
    integer (kind = 4):: nuci(300), nucntag = 0
    character (len = 30) :: m1
    logical :: quit 

    write(m1,'(i6)') icstart
    m1=trim(adjustl(m1))

    if (debug==1) write(*,*) 'INTERSTRESS: start. Cycle id is', ic

    if (ic == 1) then 
        ss0 = -ambientnorm * fric_fini
        ns0 = ambientnorm
    else
        do i = 1, totftnode
            ss0(i) = fistr(1,i) 
            ns0(i) = fistr(2,i) 
        enddo 
    endif
    ntag = 0  
    do ift = 1, ntotft
        do i = 1, nfnode(ift)
            ntag = ntag + 1
            ida(ntag) = (ift - 1)*maxftnode + i
        enddo              
    enddo

    if (debug==1) write(*,*) 'INTERSTRESS: initialization.'
    if (debug==1) write(*,*) 'INTERSTRESS: totftnode', totftnode
    if (debug==1) write(*,*) 'INTERSTRESS: nfnode', nfnode

    shs = 0.0d0
    ns = 0.0d0
    quit = .FALSE.
    dtInYears = 1.0d0 ! in years
    tcon = 365.0d0*24.0d0*3600.0d0
    if (debug==1) write(*,*) tcon
    dtInSeconds = dtInYears * tcon
    tinter = 0.0d0
    do ndt = 1, 100000
        if (quit .eqv. .TRUE.) then
            exit
        endif
        tinter = tinter + dtInSeconds
        if (mod(ndt,100) == 1) then
            write (*,*) 'Interseismic t = ', tinter/tcon, ' years'
        endif
        do i = 1, totftnode
            theta = atan(nsmpnv(1, ida(i))/nsmpnv(2,ida(i))) ! in RAD
            !theta = rd(2,ida(i))/180.0d0*pi
            !if (debug==1) write(*,*) 'INTERSTRESS: node id', i, 'theta', theta
            
            !if (theta >= 45.0d0/180.0d0*pi) then 
            !    theta = 45.0d0/180.0d0*pi
            !endif
            
            ! Viscosity eta is inversely proportional to the maximum shear strain rate (Liu et al., 2022).
            ! NOTE: don't scale for now.
            rd(1,ida(i)) = maxShearStrainLoadRate
            eta   = eta0*maxShearStrainLoadRate/rd(1,ida(i))      
            rs    = rd(1,ida(i))*cos(2.0d0*theta)*eta ! eq.(3)
            rn    = -rd(1,ida(i))*sin(2.0d0*theta)*eta ! eq.(3)
            ns(i) = (ns0(i) - ambientnorm - rn)*exp(-tinter*amu/eta)+rn+ambientnorm ! eq.(2) 
            shs(i) = (ss0(i) - rs) * exp(-tinter*amu/eta) + rs ! eq.(1)
            
            etaMin = -fric_fs*ambientnorm/maxShearStrainLoadRate
            if (eta<etaMin) then 
                write(*,*) 'Minimum viscosity is ', etaMin
                write(*,*) 'For ambient normal stress ', ambientnorm/1.d6, ' MPa'
                write(*,*) 'and static friction ', fric_fs
                write(*,*) 'maxShearStrainLoadRate is', maxShearStrainLoadRate
                write(*,*) 'Max Shr direction to strike angle is', theta/pi*180.d0
                write(*,*) 'Viscosity is ', eta
                write(*,*) 'Viscosity is lower than the minimum requirement. Exiting ...'
                stop
            endif
            
            
            strengthexcess(i) = (abs(ns(i))*fric_fs - shs(i))
            if (i==1.or.i==nfnode(1).or.i==1+nfnode(1).or. &
                            i==nfnode(1)+nfnode(2).or.i==(nfnode(1)+nfnode(2)+1).or. &
                            i==(nfnode(1)+nfnode(2)+nfnode(3))) then
                strengthexcess(i) = 100.0d6
            endif
            ! if (debug==1) then
                ! write(*,*) 'rs, rn', rs, rn
                ! write(*,*) 'ns(i), shs(i)', ns(i), shs(i)
                ! write(*,*) eta, amu
                ! write(*,*) strengthexcess(i)/1e6, abs(ns(i))/1e6, fric_fs, shs(i)/1e6
            ! endif
        enddo
        
        !if (debug==1) write(*,*) 'INTERSTRESS: finish assigning straining rate.'
        
        minstrengthexcess = minval(strengthexcess)
        loc = minloc(strengthexcess, dim = 1)
        
        !if (debug==1) write(*,*) 'INTERSTRESS: strength excess is', strengthexcess
        !if (debug==1) write(*,*) 'INTERSTRESS: loc is', loc, 'min strength is', minstrengthexcess
        
        nucntag = 0
        do i = 1, totftnode
            if ((i >1.and.i<nfnode(1)).or.(i>(1+nfnode(1)).and. &
                            i<(nfnode(1)+nfnode(2))).or.(i>(nfnode(1)+nfnode(2)+1).and. &
                            (i<(nfnode(1)+nfnode(2)+nfnode(3))))) then
                if (shs(i) > fric_fs * abs(ns(i))) then 
                    nucntag = nucntag + 1
                    nuci(nucntag) = i
                    quit  = .TRUE.
                endif
            endif
        enddo
    enddo

    do i = 1, totftnode
        fistr(1,i) = shs(i)
        fistr(2,i) = ns(i)
    enddo 

    if (debug==1) write(*,*) 'tcon', tcon
    if (debug==1) write(*,*) 'tinter is', tinter

    write(*,*) '=                                                                   ='
    write(*,*) '=     Interseismic takes                                            ='
    write(*,'(X,A,40X,f15.2,3X,5A)')  '=', tinter/tcon, 'years'
    write(*,*) '=     Nucleation occurs at on-fault node                            ='
    write(*,'(X,A,40X,i5)') '=',loc
    write(*,*) 'normal, shear, strength excess are', ns(loc)/1e6, shs(loc)/1e6, strengthexcess(i)/1e6

    open(2001,file='interval.txt'//m1, position='append')
        write(2001,'(1e21.14)') tinter/tcon
    close(2001)

    if (debug==1) write(*,*) 'INTERSTRESS: exiting.'
    ! if (nucntag > 1) then 
        ! write(*,*) 'Many nucleation points', nucntag
        ! write(*,*) (nuci(i), i = 1, nucntag) 
    ! endif
end subroutine interstress
