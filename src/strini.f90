program strini
! Initial stress for Pingding.04/11/2019
implicit none
integer (kind=4) :: totft = 3, nft(3), id(5307), maxftnode = 1769, ift, ntag, ntotnd,  ndt, nnode, i,nct,nuc
real (kind = 8) :: nv(3, 5307),rd(2,5307) 
real (kind=8) :: amu, ant0,ant, str, miu0, mius, tx,ty, dt,t,tcon,temp1,temp2,temp3, rn, rs
real (kind = 8) :: theta, j, minstrengthexcess, ss0
real (kind = 8),allocatable,dimension(:) :: ss, ns, strengthexcess,nv1, sf1s,sf1n,sf1s0,sf1n0,sf1na
real (kind=8),allocatable,dimension(:) :: sf2s,sf2n,sf2s0,sf2n0,sf2na
real (kind=8),allocatable,dimension(:) :: theta1,theta2,rn1,rt1,rn2,rt2
integer (kind = 4)::nuci(100), nucntag = 0
logical :: ny = .true.
logical :: quit 
amu = 3.20d10 ! Shear modulus
ant0 = 2.8d21*2.0d0 ! Viscosity for loading rate of 450 nrad/yr. 
str = 1.427d-14 ! Maximum shearing strain loading rate
!This loading rate equals to 450 nrad/yr.
miu0 = 0.0d0
mius = 0.4d0 
        open(1, file = './mesh/Mesh_general_info.txt', form = 'formatted', status = 'old') 
                read(1,*) temp1, temp2
                read(1,*) (nft(i), i = 1, totft)
        close(1)
        open(2, file = './mesh/nsmpnv.txt', form = 'formatted', status = 'old')
                read(2,*) (nv(1,i), nv(2,i),nv(3,i), i = 1, maxftnode*totft)
        close(2)
        open(3, file = 'Rate_direction.txt', form = 'formatted', status = 'old')
                read(3,*) (rd(1,i),rd(2,i),i=1,maxftnode*totft)
        close(3)
        ntotnd = sum(nft)
        write(*,*) 'Total fault nodes = ', ntotnd
        write(*,*) (nft(i), i = 1, totft )
        allocate(ss(ntotnd), ns(ntotnd), strengthexcess(ntotnd))
        ntag = 0  
        do ift = 1, totft
                do i = 1, nft(ift)
                        ntag = ntag + 1
                        id(ntag) = (ift - 1)* maxftnode + i
                enddo              
        enddo
        !initialize
        ss = 0.0d0
        ns = 0.0d0
quit = .FALSE.
dt = 1.0d0 ! in years
tcon = 365.0d0 * 24.0d0 *3600.0d0
dt = dt * tcon
t = 0.0d0
do ndt = 1, 10000000
        t = t + dt
        !write (*,*) 't = ', t/tcon, ' years'
        if (quit == .TRUE.) then
                exit
        endif
        do i = 1, ntotnd 
                theta = atan(nv(1, id(i))/nv(2,id(i))) ! in RAD
                !theta = acos(nv(2,id(i))) 
                theta = -theta + 11.71d0/180.0d0*3.1415926d0 - rd(2,id(i))/180.0d0*3.1415926d0
                if (theta>=45.0d0/180.0d0*3.1415926d0) then
                        theta = 45.0d0/180.0d0*3.1415926d0
                endif
                ant = ant0*450.0d0/rd(1,id(i))
		rs = str/450.0d0*rd(1,id(i)) * cos(2.0d0*theta) * ant
                rn = - str/450.0d0*rd(1,id(i)) * sin(2.0d0*theta) *ant
		
                ss0 = (rn + 100.0d6)*miu0
                ns(i) = - rn * exp(-t*amu/ant) + rn + -100.0d6 
                ss(i) = -rs * exp(-t*amu/ant) + rs
                strengthexcess(i) = (abs(ns(i))*mius - ss(i))
		if (i==1.or.i==295.or.i==296.or.i==473.or.i==474.or.i==2242) then
			strengthexcess = 100.0d6
		endif
        enddo
        minstrengthexcess = minval(strengthexcess)
        if (minstrengthexcess<1.0d6) then 
        !        dt = 0.05d0 * tcon
        endif
        do i = 1, ntotnd
		temp3 = mius
		if ((i > 1.or.i<nft(1)).or.(i>(1+nft(1)).and.i<(nft(1)+nft(2))).or.(i>(nft(1)+nft(2)+1).and.i<(nft(1)+nft(2)+nft(3)))) then
                	if (ss(i) >= mius * abs(ns(i))) then 
                        	nucntag = nucntag + 1
                        	nuci(nucntag) = i
                        	quit  = .TRUE.
                	endif
		endif
        enddo
        ! if (mod(ndt,10) == 1) then
                ! open(3, file = 'stress4plot.txt', form = 'formatted', position = 'append')
                        ! write(3,'(2e18.7e4)') (ss(j), ns(j), j = 1, ntotnd)  
        ! endif
enddo
  !assign values

write(*,*) 'mim strengthexcess', minstrengthexcess/1e6, 'in MPa'
  
open(4, file='finalstress.txt', form = 'formatted', status='unknown')
        write(4, '(2e18.7e4)') (ss(i), ns(i), i = 1, ntotnd)
close(4) 
  
open(5,file='interval.time', position='append')
        write(5,'(1e18.7e4)') t/tcon
close(5)

open(6, file = 'nucloc.txt', status = 'unknown')
        write(6, *) 1, (nuci(i), i =1, nucntag)
close(6)
if (nucntag > 1) then 
	write(*,*) 'Many nucleation points', nucntag
	write(*,*) 
endif
write(*,*) 't ==', t/tcon, 'years'
end program strini  
      
