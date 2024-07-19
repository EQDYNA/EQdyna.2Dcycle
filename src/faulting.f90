SUBROUTINE faulting(step)
    use globalvar
    implicit none

    logical :: lstr
    integer (kind=4)::ifault,i,j,k, ii, jj, kk, step
    real (kind = dp)::nx,ny,tx,ty,txy,mmast,mslav,mtotl,ttao,tnrm,taoc, &
        taox,taoy,ftix, ftiy, fnfault, ftfault, slip,sliprate,xmu,trupt,tr,temp1, maxslip, maxsliprate
    real (kind = dp),dimension(4,2,3)::fvd=0.0d0 
    real (kind = dp):: xcoor0, ycoor0, ift0, radi
    real (kind = dp),dimension(totftnode) :: maxslip_arr, maxsliprate_arr
    !
    if (loc < nfnode(1)+1) then
        ift0 = 1
        xcoor0 = x(1, nsmp0(1, loc))
        ycoor0 = x(2, nsmp0(1, loc))
    elseif ((loc < nfnode(1) + nfnode(2) + 1).and.(loc > nfnode(1))) then
        ift0 =2
        xcoor0 = x(1, nsmp0(1, loc-nfnode(1)+maxftnode))
        ycoor0 = x(2, nsmp0(1, loc-nfnode(1)+maxftnode))
    elseif ((loc < nfnode(1) + nfnode(2) + nfnode(3) + 1).and.(loc > nfnode(1) + nfnode(2))) then    
        ift0 =3
        xcoor0 = x(1,nsmp0(1,loc-nfnode(1)-nfnode(2)+ 2*maxftnode))
        ycoor0 = x(2,nsmp0(1,loc-nfnode(1)-nfnode(2)+ 2*maxftnode))
    elseif ((loc < nfnode(1) + nfnode(2) + nfnode(3) + nfnode(4) + 1).and.(loc > nfnode(1) + nfnode(2) + nfnode(3))) then    
        ift0 =4
        xcoor0 = x(1,nsmp0(1,loc-nfnode(1)-nfnode(2)-nfnode(3)+3*maxftnode))
        ycoor0 = x(2,nsmp0(1,loc-nfnode(1)-nfnode(2)-nfnode(3)+3*maxftnode))
    endif
    !write(*,*) ift0
    !*** loop over slave nodes ***
    kk = 0
    do ii=1,ntotft
        do jj = 1, nfnode(ii)
            i = (ii-1)*maxftnode + jj
            kk = kk + 1
            !...get unit normal, shear: use variables for easily coding. B.D. 8/25/06
            ! AN error: I imported nx,ny while Dr. Duan imported tx,ty.
            nx = -nsmpnv(1,i)
            ny = -nsmpnv(2,i)
            txy = nsmpnv(3,i)
            ftfault = -fistr(1,kk)
            fnfault = fistr(2,kk)
            tx = -ny ! The sign of tx,ty are reversed from those in Cajonpass Cyc2d_v3a.
            ty = nx
            !nx = ty
            !ny = -tx
            !
            !...get nodal force,velocity, and diplacement in x,y.
            !   B.D. 11/23/06
            !...aslo add Rayleigh stiffness damping before using d.
            !   assume stifness coefficient is first material: rdampk(1).
            !   B.D. 11/26/06
            do j=1,2  !1-slave, 2-master
                temp1=(-1)**(j-1)
                do k=1,2  !1-x comp, 2-y comp
                    fvd(k,j,1) = brhs(id(k,nsmp0(j,i))) !+ temp1*fistr(k,i)  !1-force
                    fvd(k,j,2) = v(k,nsmp0(j,i)) !2-vel
                    fvd(k,j,3) = d(k,nsmp0(j,i)) + rdampk(1)*fvd(k,j,2) !3-disp
                enddo
            enddo
            !...resolve x,y components onto normal and shear components.
            !   B.D. 11/23/06
            do j=1,3    !1-force,2-vel,3-disp
                do k=1,2  !1-slave,2-master
                    fvd(3,k,j) = fvd(1,k,j)*nx + fvd(2,k,j)*ny  !3-norm
                    fvd(4,k,j) = fvd(1,k,j)*tx + fvd(2,k,j)*ty  !4-shear
                enddo
            enddo
            !
            !...nodal mass. Mass of each element may not be distributed among its 
            ! nodes evenly. Instead, distribution is related to element shape. 
            !   Note: nodal mass should not be directly obtained from left-hand-side
            ! diagnoal mass matrix, because that's the effective mass, which takes 
            ! damping coefficient into accout. Instead, I computed nodal mass from 
            ! element mass and assembled in "qdct2.f90".B.D.7/3/05
            mslav = fnms(nsmp0(1,i))        
            mmast = fnms(nsmp0(2,i))
            mtotl = mslav + mmast
            mtotl = mtotl * txy
            !
            !...trial traction to enforce continuity. B.D. 11/23/06
            ttao = (mslav*mmast*(fvd(4,1,2)-fvd(4,2,2))/dt1 + mmast*fvd(4,1,1) &
                    - mslav*fvd(4,2,1)) / mtotl + ftfault
            tnrm = (mslav*mmast*((fvd(3,1,2)-fvd(3,2,2))+(fvd(3,1,3)-fvd(3,2,3))/dt1)/dt1 &
                    + mmast*fvd(3,1,1) - mslav*fvd(3,2,1)) / mtotl + fnfault       
            !
            !...friction law to determine friction coefficient
            !   B.D. 11/23/06
            ! for slip-weakening:
            slip = sqrt((fvd(1,2,3)-fvd(1,1,3))**2+(fvd(2,2,3)-fvd(2,1,3))**2) !slip mag
            sliprate = sqrt((fvd(1,2,2)-fvd(1,1,2))**2+(fvd(2,2,2)-fvd(2,1,2))**2) !sliprate mag

            !... based on choices, call corresponding friction laws.
            ! B.D. 8/19/06
            if(friclaw == 1) then
                call slip_weak(slip,fric(1,i),xmu)
            elseif(friclaw == 2) then
                trupt =  timedyna - fnft(i)
                call time_weak(trupt,fric(1,i),xmu)
            elseif(friclaw==3)then
                call slip_weak_healing(slip,fric(1,i),xmu)
            elseif(friclaw==4)then
                call slip_rate_weak(slip,sliprate,fric(1,i),xmu)
            elseif(friclaw==5)then
                call slip_rate_caltech(slip,sliprate,fric(1,i),xmu)
            endif

            if(ii == ift0 .and. xmu > fric(2,i).and.fric(1,i)<500.0d0) then
              !only nucleation fault and before finishing dropping, do...
              radi = ((x(1,nsmp0(1,i)) - xcoor0)**2 + (x(2,nsmp0(1,i)) - ycoor0)**2)**0.5
              if(radi <= radius) then !only within nucleation zone, do...
                tr = radi / vrupt0
                if(tr <= timedyna) then !only ready or already fail, do...
                  trupt = timedyna - tr
                  call time_weak(trupt,fric(1,i),xmu)
                endif
              endif
            endif                          
            if (ii == ift0.and.abs(x(1,nsmp0(1,i))-xcoor0)<0.1d0 .and. &
                abs(x(2,nsmp0(1,i))-ycoor0)<0.1d0 .and. &
                (step>1 .and. mod(step,2)==numOfStepsToPrint)) then
                write(*,*) '=                                                                   ='
                write(*,*) '=     Slipx at epicenter                                            ='
                write(*,'(X,A,40X,E15.7)') '=',(fvd(1,2,3)-fvd(1,1,3))
                write(*,*) '=     Sliprate at epicenter                                         ='
                write(*,'(X,A,40X,E15.7)') '=',sliprate
                write(*,*) fvd
                write(*,*) 'traction ttao, tnrm', ttao/1e6, tnrm/1e6
                write(*,*) mslav, mmast
                write(*,*) xcoor0/1e3, ycoor0/1e3
                write(*,*) 'ift, node id in a fault ', ii, jj
                write(*,*) 'total ft node id ', i
                write(*,*) 'nsmp node ids', nsmp0(1,i), nsmp0(2,i) 
            endif    
            !
            !...adjust ttao and tnrm based on jump conditions on fault.
            !   before calculate taoc, first adjust tnrm if needed. 
            !   after this, they are true (corrected) values. B.D. 11/23/06
            !   if(tnrm > 0) tnrm = 0   !norm must be <= 0, otherwise no adjust
            if (tnrm > minnorm) then 
                tnrm = minnorm
            !elseif (tnrm < -100.0d6) then
            !    tnrm = -100.0d6
            endif
            taoc = -xmu * tnrm
            if(abs(ttao) >= taoc) then
                if(ttao==0) then
                    write(*,*) 'ERROR: ttao=0!! at nfn=',i, ftfault, fnfault
                endif
                ttao=taoc*ttao/abs(ttao)   !otherwise, no adjust
                if(fnft(i)>900.d0 .and. sliprate>0.001d0) fnft(i)=timedyna  !failure time for time-weakening
            endif
            if(tnrm > 0.0d0) ttao = 0.0d0
            
            output4plot(1,kk) = -ttao
            output4plot(2,kk) = tnrm
            output4plot(3,kk) = slip
            output4plot(4,kk) = sliprate
            output4plot(5,kk) = fnft(i)
            maxslip_arr(kk) = slip
            maxsliprate_arr(kk) = sliprate
            !
            !...add the above fault boundary traction to the split nodes.
            !   first resolve normal and shear back to x-,y-. 
            !   then subtract them from slave, add to master as the above calculation
            !   based on this convention. see Day et al. (2005). B.D. 11/23/06
            taox = (tnrm*nx + ttao*tx) * txy
            taoy = (tnrm*ny + ttao*ty) * txy
            ftix = (fnfault*nx + ftfault*tx) * txy
            ftiy = (fnfault*ny + ftfault*ty) * txy
            brhs(id(1,nsmp0(1,i))) = brhs(id(1,nsmp0(1,i))) - taox + ftix
            brhs(id(2,nsmp0(1,i))) = brhs(id(2,nsmp0(1,i))) - taoy + ftiy
            brhs(id(1,nsmp0(2,i))) = brhs(id(1,nsmp0(2,i))) + taox - ftix
            brhs(id(2,nsmp0(2,i))) = brhs(id(2,nsmp0(2,i))) + taoy - ftiy
            !
            !...Store fault stresses and slip/slipvel for fault nodes at set time interval.
            ! note: need element dimension to calculate stress from force
            ! B.D. 7/4/05
            ! directly write out for monitoring. B.D. 8/20/06
            !...modify to be conisistent with the above change. B.D. 11/23/06
            !
        enddo    !ending loop on slave nodesa
    enddo
    maxslip = maxval(maxslip_arr)
    maxsliprate = maxval(maxsliprate_arr)
    if (maxsliprate < 0.001d0 .and. step*dt1 >5.0d0) then
        write(*,*) '=                                                                   ='
        write(*,*) '=     Time to quit dynamic rupture                                  ='
        write(*,*) '=                                                                   ='
        write(*,*) '=     Maximum sliprate, m/s                                         ='    
        write(*,'(X,A,40X,E15.7)') '=', maxsliprate
        write(*,*) '=     Maximum slip, m                                               ='
        write(*,'(X,A,40X,E15.7)') '=', maxslip    
        quitdriver = .TRUE.
    endif

end SUBROUTINE faulting     
