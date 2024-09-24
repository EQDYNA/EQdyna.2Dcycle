subroutine interstress(cycleID)

    use globalvar
    implicit none

    integer (kind=4) :: cycleID
    integer (kind=4) :: ftPairIDtoNodeID(maxftnode*ntotft)
    integer (kind=4) :: ift, ntag, i, timeStepID
    real (kind = dp) :: shearStress(totNumFtNode)
    real (kind = dp) :: normStress(totNumFtNode)
    real (kind = dp) :: initShear(totNumFtNode)
    real (kind = dp) :: initNorm(totNumFtNode)
    real (kind = dp) :: strengthExcess(totNumFtNode)
    real (kind = dp) :: eta, etaMin
    real (kind = dp) :: dtInSeconds, dtInYears
    real (kind = dp) :: interseisElapsedTime
    real (kind = dp) :: conYrToSec
    real (kind = dp) :: dNorm, dShear
    real (kind = dp) :: minusLocStrikeAngleToX
    real (kind = dp) :: locMaximumShearRate, locLoadAngleToX, locLoadWeight, locFtType, locFtDip
    real (kind = dp) :: minStrengthExcess
    integer (kind = 4):: nuci(300), nucntag = 0
    character (len = 30) :: m1
    logical :: quit 

    write(m1,'(i6)') icstart
    m1=trim(adjustl(m1))

    if (debug==1) write(*,*) 'INTERSTRESS: start. Cycle id is', cycleID

    if (cycleID == 1) then 
        initShear = -ambientNorm * fric_fini
        initNorm = ambientNorm
    else
        initShear = fistr(1,:)
        initNorm = fistr(2,:)
        !do i = 1, totNumFtNode
        !    initShear(i) = fistr(1,i) 
        !    initNorm(i) = fistr(2,i) 
        !enddo 
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
    do timeStepID = 1, 100000
        if (quit .eqv. .TRUE.) then
            exit
        endif
        interseisElapsedTime = interseisElapsedTime + dtInSeconds
        if (mod(timeStepID,100) == 1) then
            write (*,*) 'Interseismic t = ', interseisElapsedTime/conYrToSec, ' years'
        endif
        do i = 1, totNumFtNode
            ! NOTE: the minus sign is good for left-lateral loading and faults.
            minusLocStrikeAngleToX = atan(nsmpGeoPhys(1,ftPairIDtoNodeID(i))/nsmpGeoPhys(2,ftPairIDtoNodeID(i)))
            !minusLocStrikeAngleToX = rd(2,ftPairIDtoNodeID(i))/180.0d0*pi
            !if (debug==1) write(*,*) 'INTERSTRESS: node id', i, 'minusLocStrikeAngleToX', minusLocStrikeAngleToX
            
            !if (minusLocStrikeAngleToX >= 45.0d0/180.0d0*pi) then 
            !    minusLocStrikeAngleToX = 45.0d0/180.0d0*pi
            !endif
            
            ! Viscosity eta is inversely proportional to the maximum shear strain rate (Liu et al., 2022).
            ! NOTE: don't scale for now.
            !rd(1,ftPairIDtoNodeID(i)) = maxShearStrainLoadRate
            !eta   = eta0*maxShearStrainLoadRate/rd(1,ftPairIDtoNodeID(i))
            !rs    = rd(1,ftPairIDtoNodeID(i))*cos(2.0d0*minusLocStrikeAngleToX)*eta ! eq.(3)
            !rn    = -rd(1,ftPairIDtoNodeID(i))*sin(2.0d0*minusLocStrikeAngleToX)*eta ! eq.(3)
            locFtType = nsmpGeoPhys(4,ftPairIDtoNodeID(i))
            locFtDip = nsmpGeoPhys(5,ftPairIDtoNodeID(i))
            locMaximumShearRate = nsmpGeoPhys(6,ftPairIDtoNodeID(i))
            locLoadAngleToX = nsmpGeoPhys(7,ftPairIDtoNodeID(i))
            locLoadWeight = nsmpGeoPhys(8,ftPairIDtoNodeID(i))
            eta = nsmpGeoPhys(9,ftPairIDtoNodeID(i))     
            dShear = locMaximumShearRate*locLoadWeight*cos(2.0d0*minusLocStrikeAngleToX)*eta ! eq.(3)
            dNorm = -locMaximumShearRate*locLoadWeight*sin(2.0d0*minusLocStrikeAngleToX)*eta ! eq.(3)
            normStress(i) = (initNorm(i)-ambientNorm-dNorm)*exp(-interseisElapsedTime*amu/eta)+dNorm+ambientNorm ! eq.(2) 
            shearStress(i) = (initShear(i)-dShear) * exp(-interseisElapsedTime*amu/eta)+dShear ! eq.(1)
            
            etaMin = -fric_fs*ambientNorm/maxShearStrainLoadRate
            if (eta<etaMin) then 
                write(*,*) 'Minimum viscosity is ', etaMin
                write(*,*) 'For ambient normal stress ', ambientNorm/1.d6, ' MPa'
                write(*,*) 'and static friction ', fric_fs
                write(*,*) 'maxShearStrainLoadRate is', maxShearStrainLoadRate
                write(*,*) 'Max Shr direction to strike angle is', minusLocStrikeAngleToX/pi*180.d0
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
        enddo

        minStrengthExcess = minval(strengthExcess)
        loc = minloc(strengthExcess, dim = 1)
  
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
