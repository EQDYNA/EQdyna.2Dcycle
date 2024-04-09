program strbld
! Initial stress for Pingding.04/11/2019
implicit none
integer (kind=4) :: totft = 3, nft(3), id(5307), maxftnode = 1769, ift, ntag, ntotnd,  ndt, nnode, i, loc(1), nct,nuc
real (kind = 8) :: nv(3, 5307),rd(2,5307) 
real (kind=8) :: amu, ant0,ant, str, miu0, mius, tx,ty, dt,t,tcon,temp1,temp2,temp3, rn, rs
real (kind = 8) :: theta, j, minstrengthexcess, tmp1
real (kind = 8),allocatable,dimension(:) :: ss, ns, strengthexcess,nv1, ss0, ns0
real (kind=8),allocatable,dimension(:) :: sf2s,sf2n,sf2s0,sf2n0,sf2na
real (kind=8),allocatable,dimension(:) :: theta1,theta2,rn1,rt1,rn2,rt2
integer (kind = 4)::nuci(300), nucntag = 0
logical :: ny = .true.
logical :: quit 
amu = 3.20d10 ! Shear modulus
ant0 = 2.8d21*2.0d0 ! Viscosity 
str = 1.427d-14 ! Maximum shearing strain loading rate
mius = 0.4d0 
        open(1, file = './mesh/Mesh_general_info.txt', form = 'formatted', status = 'old') 
                read(1,*) temp1, temp2
                read(1,*) (nft(i), i = 1, totft)
        close(1)
        open(2, file = './mesh/nsmpnv.txt', form = 'formatted', status = 'old')
                read(2,*) (nv(1,i), nv(2,i),nv(3,i), i = 1, maxftnode*totft)
        close(2)
        open(3,file='Rate_direction.txt',form = 'formatted',status = 'old')
                read(3,*) (rd(1,i),rd(2,i),i=1,maxftnode*totft)
        close(3)
		
        ntotnd = sum(nft)
        write(*,*) 'Total fault nodes = ', ntotnd
        write(*,*) (nft(i), i = 1, totft )
        allocate(ss(ntotnd), ns(ntotnd), ss0(ntotnd), ns0(ntotnd), strengthexcess(ntotnd))        
		open(3, file = 'output4plot_dy.txt', form = 'formatted', status = 'old')
                read(3,*) (ss0(i), ns0(i), tmp1, tmp1, tmp1, i = 1, ntotnd)
        close(3)
        ntag = 0  
        do ift = 1, totft
                do i = 1, nft(ift)
                        ntag = ntag + 1
                        id(ntag) = (ift - 1)*maxftnode + i
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
do ndt = 1, 1000000000
        if (quit == .TRUE.) then
		exit
	endif
	t = t + dt
	if (mod(ndt,50) == 1) then
        	write (*,*) 't = ', t/tcon, ' years'
	endif
        do i = 1, ntotnd 
                theta = atan(nv(1, id(i))/nv(2,id(i))) ! in RAD
                theta = -theta + 11.71d0/180.0d0*3.1415926d0 - rd(2,id(i))/180.0d0*3.1415926d0
                if (theta >= 45.0d0/180.0d0*3.1415926d0) then 
                        theta = 45.0d0/180.0d0*3.1415926d0
                endif   
                ant = ant0*450.0d0/rd(1,id(i))     
                rs = str/450.0d0*rd(1,id(i)) * cos(2.0d0*theta) * ant
                rn = -str/450.0d0*rd(1,id(i)) * sin(2.0d0*theta) * ant
                ns(i) = (ns0(i) - -100.0d6 - rn) * exp(-t*amu/ant) + rn + -100.0d6 
                ss(i) = (ss0(i) - rs) * exp(-t*amu/ant) + rs 
                strengthexcess(i) = (abs(ns(i))*mius - ss(i))
		if (i==1.or.i==nft(1).or.i==1+nft(1).or.i==nft(1)+nft(2).or.i==(nft(1)+nft(2)+1).or.i==(nft(1)+nft(2)+nft(3))) then
			strengthexcess(i) = 100.0d6
		endif
        enddo
        minstrengthexcess = minval(strengthexcess)
	loc(1) = minloc(strengthexcess, dim = 1)
        if (minstrengthexcess<1.0d6) then 
!                dt = 1.0d0 * tcon
        endif
        do i = 1, ntotnd
	        temp3 = mius
                if ((i >1.and.i<nft(1)).or.(i>(1+nft(1)).and.i<(nft(1)+nft(2))).or.(i>(nft(1)+nft(2)+1).and.(i<(nft(1)+nft(2)+nft(3))))) then
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
  
open(4, file='finalstress.txt', form = 'formatted', status='unknown')
        write(4, '(2e18.7e4)') (ss(i), ns(i), i = 1, ntotnd)
close(4) 
  
open(5,file='interval.time', position='append')
        write(5,'(1e18.7e4)') t/tcon
close(5)

write(*,*) '========'
write(*,*) 'nuc loc =', loc(1) , 'min shearstrengthexcess',strengthexcess(loc(1))/1.0d6,'in MPa'
write(*,*) '========'
open(6, file = 'nucloc.txt', status = 'unknown')
        write(6, *) 1, loc(1)
close(6)
if (nucntag > 1) then 
	write(*,*) 'Many nucleation points', nucntag
	write(*,*) (nuci(i), i = 1, nucntag) 
endif

open(7, file ='numofnucpoints.txt', position ='append')
	write(7,*) nucntag
close(7)  
end program strbld
      
