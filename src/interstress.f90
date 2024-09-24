subroutine interstress(ic)

    use globalvar
    implicit none

    integer (kind=4) :: ic, ftPairIDtoNodeID(maxftnode*ntotft), ift, ntag, i, j, nuc, ndt
    real (kind = dp) :: shearStress(totNumFtNode), normStress(totNumFtNode), ss0(totNumFtNode), ns0(totNumFtNode), &
        strengthExcess(totNumFtNode) 
    real (kind = dp) :: eta, etaMin, dtInSeconds, dtInYears, interseisElapsedTime, conYrToSec, rn, rs
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
        do i = 1, totNumFtNode
            ss0(i) = fistr(1,i) 
            ns0(i) = fistr(2,i) 
        enddo 
    endif
    ntag = 0  
    do ift = 1, ntotft
        do i = 1, nfnode(ift)
            ntag = ntag + 1
            ftPairIDtoNodeID(ntag) = (ift - 1)*maxftnode + i
        enddo              
    enddo

    if (debug==1) write(*,*) 'INTERSTRESS: initialization.'
    if (debug==1) write(*,*) 'INTERSTRESS: totNumFtNode', totNumFtNode
    if (debug==1) write(*,*) 'INTERSTRESS: nfnode', nfnode

    shearStress = 0.0d0
    normStress = 0.0d0
    quit = .FALSE.
    dtInYears = 1.0d0 ! in years
    conYrToSec = 365.0d0*24.0d0*3600.0d0
    if (debug==1) write(*,*) conYrToSec
    dtInSeconds = dtInYears * conYrToSec
    interseisElapsedTime = 0.0d0
    do ndt = 1, 100000
        if (quit .eqv. .TRUE.) then
            exit
        endif
        interseisElapsedTime = interseisElapsedTime + dtInSeconds
        if (mod(ndt,100) == 1) then
            write (*,*) 'Interseismic t = ', interseisElapsedTime/conYrToSec, ' years'
        endif
        do i = 1, totNumFtNode
            theta = atan(nsmpGeoPhys(1, ftPairIDtoNodeID(i))/nsmpGeoPhys(2,ftPairIDtoNodeID(i))) ! in RAD
            !theta = rd(2,ftPairIDtoNodeID(i))/180.0d0*pi
            !if (debug==1) write(*,*) 'INTERSTRESS: node id', i, 'theta', theta
            
            !if (theta >= 45.0d0/180.0d0*pi) then 
            !    theta = 45.0d0/180.0d0*pi
            !endif
            
            ! Viscosity eta is inversely proportional to the maximum shear strain rate (Liu et al., 2022).
            ! NOTE: don't scale for now.
            !rd(1,ftPairIDtoNodeID(i)) = maxShearStrainLoadRate
            !eta   = eta0*maxShearStrainLoadRate/rd(1,ftPairIDtoNodeID(i))
            !rs    = rd(1,ftPairIDtoNodeID(i))*cos(2.0d0*theta)*eta ! eq.(3)
            !rn    = -rd(1,ftPairIDtoNodeID(i))*sin(2.0d0*theta)*eta ! eq.(3)
            
            eta = nsmpGeoPhys(9,ftPairIDtoNodeID(i))     
            rs    = nsmpGeoPhys(6,ftPairIDtoNodeID(i)) *cos(2.0d0*theta)*eta ! eq.(3)
            rn    = -nsmpGeoPhys(6,ftPairIDtoNodeID(i)) *sin(2.0d0*theta)*eta ! eq.(3)
            normStress(i) = (ns0(i) - ambientnorm - rn)*exp(-interseisElapsedTime*amu/eta)+rn+ambientnorm ! eq.(2) 
            shearStress(i) = (ss0(i) - rs) * exp(-interseisElapsedTime*amu/eta) + rs ! eq.(1)
            
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
            
            
            strengthexcess(i) = (abs(normStress(i))*fric_fs - shearStress(i))
            if (i==1.or.i==nfnode(1).or.i==1+nfnode(1).or. &
                            i==nfnode(1)+nfnode(2).or.i==(nfnode(1)+nfnode(2)+1).or. &
                            i==(nfnode(1)+nfnode(2)+nfnode(3))) then
                strengthExcess(i) = 100.0d6
            endif
            ! if (debug==1) then
                ! write(*,*) 'rs, rn', rs, rn
                ! write(*,*) 'ns(i), shearStress(i)', ns(i), shearStress(i)
                ! write(*,*) eta, amu
                ! write(*,*) strengthexcess(i)/1e6, abs(ns(i))/1e6, fric_fs, shearStress(i)/1e6
            ! endif
        enddo
        
        !if (debug==1) write(*,*) 'INTERSTRESS: finish assigning straining rate.'
        
        minStrengthExcess = minval(strengthExcess)
        loc = minloc(strengthExcess, dim = 1)
        
        !if (debug==1) write(*,*) 'INTERSTRESS: strength excess is', strengthexcess
        !if (debug==1) write(*,*) 'INTERSTRESS: loc is', loc, 'min strength is', minstrengthexcess
        
        nucntag = 0
        do i = 1, totNumFtNode
            if ((i >1.and.i<nfnode(1)).or.(i>(1+nfnode(1)).and. &
                            i<(nfnode(1)+nfnode(2))).or.(i>(nfnode(1)+nfnode(2)+1).and. &
                            (i<(nfnode(1)+nfnode(2)+nfnode(3))))) then
                if (shearStress(i) > fric_fs * abs(normStress(i))) then 
                    nucntag = nucntag + 1
                    nuci(nucntag) = i
                    quit  = .TRUE.
                endif
            endif
        enddo
    enddo

    do i = 1, totNumFtNode
        fistr(1,i) = shearStress(i)
        fistr(2,i) = normStress(i)
    enddo 

    if (debug==1) write(*,*) 'conYrToSec', conYrToSec
    if (debug==1) write(*,*) 'interseisElapsedTime is', interseisElapsedTime

    write(*,*) '=                                                                   ='
    write(*,*) '=     Interseismic takes                                            ='
    write(*,'(X,A,40X,f15.2,3X,5A)')  '=', interseisElapsedTime/conYrToSec, 'years'
    write(*,*) '=     Nucleation occurs at on-fault node                            ='
    write(*,'(X,A,40X,i5)') '=',loc
    write(*,*) 'normal, shear, strength excess are', normStress(loc)/1e6, shearStress(loc)/1e6, strengthExcess(i)/1e6

    open(2001,file='interval.txt'//m1, position='append')
        write(2001,'(1e21.14)') interseisElapsedTime/conYrToSec
    close(2001)

    if (debug==1) write(*,*) 'INTERSTRESS: exiting.'
    ! if (nucntag > 1) then 
        ! write(*,*) 'Many nucleation points', nucntag
        ! write(*,*) (nuci(i), i = 1, nucntag) 
    ! endif
end subroutine interstress
