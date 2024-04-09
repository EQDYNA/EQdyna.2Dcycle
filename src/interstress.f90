subroutine interstress(ic)

use globalvar
implicit none

integer (kind=4) :: ic, ida(maxftnode*ntotft), ift, ntag, i, j, nuc, ndt
real (kind = dp) :: shs(totftnode), ns(totftnode), ss0(totftnode), ns0(totftnode), &
	strengthexcess(totftnode) 
real (kind = dp) :: ant, dtinter, tinter, tcon, rn, rs
real (kind = dp) :: theta, minstrengthexcess
real (kind = dp),allocatable,dimension(:) :: theta1,theta2,rn1,rt1,rn2,rt2
integer (kind = 4)::nuci(300), nucntag = 0
character (len = 30) :: m1
logical :: quit 

write(m1,'(i6)') icstart
m1=trim(adjustl(m1))

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
!initialize
shs = 0.0d0
ns = 0.0d0
quit = .FALSE.
dtinter = 1.0d0 ! in years
tcon = 365.0d0 * 24.0d0 *3600.0d0
dtinter = dtinter * tcon
tinter = 0.0d0
do ndt = 1, 1000000000
	if (quit == .TRUE.) then
		exit
	endif
	tinter = tinter + dtinter
	! if (mod(ndt,100) == 1) then
        ! write (*,*) 'Interseismic t = ', tinter/tcon, ' years'
	! endif
	do i = 1, totftnode
		!theta = atan(nsmpnv(1, ida(i))/nsmpnv(2,ida(i))) ! in RAD
		theta = rd(2,ida(i))/180.0d0*pi
		if (theta >= 45.0d0/180.0d0*pi) then 
			theta = 45.0d0/180.0d0*pi
		endif   
		ant = ant0*str/rd(1,ida(i))     
		rs = rd(1,ida(i)) * cos(2.0d0*theta) * ant
		rn = -rd(1,ida(i)) * sin(2.0d0*theta) * ant
		ns(i) = (ns0(i) - ambientnorm - rn) * exp(-tinter*amu/ant) + rn + ambientnorm
		shs(i) = (ss0(i) - rs) * exp(-tinter*amu/ant) + rs 
		strengthexcess(i) = (abs(ns(i))*fric_fs - shs(i))
		if (i==1.or.i==nfnode(1).or.i==1+nfnode(1).or.i==nfnode(1)+nfnode(2).or.i==(nfnode(1)+nfnode(2)+1).or.i==(nfnode(1)+nfnode(2)+nfnode(3))) then
			strengthexcess(i) = 100.0d6
		endif
	enddo
	minstrengthexcess = minval(strengthexcess)
	loc = minloc(strengthexcess, dim = 1)
	nucntag = 0
	do i = 1, totftnode
		if ((i >1.and.i<nfnode(1)).or.(i>(1+nfnode(1)).and.i<(nfnode(1)+nfnode(2))).or.(i>(nfnode(1)+nfnode(2)+1).and.(i<(nfnode(1)+nfnode(2)+nfnode(3))))) then
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
write(*,*) '=                                                                   ='
write(*,*) '=     Interseismic takes                                            ='
write(*,'(X,A,40X,f7.2,3X,5A)')  '=', tinter/tcon, 'years'
write(*,*) '=     Nucleation occurs at on-fault node                            ='
write(*,'(X,A,40X,i5)') '=',loc

open(2001,file='interval.txt'//m1, position='append')
	write(2001,'(1e21.14)') tinter/tcon
close(2001)


! if (nucntag > 1) then 
	! write(*,*) 'Many nucleation points', nucntag
	! write(*,*) (nuci(i), i = 1, nucntag) 
! endif
end subroutine interstress
      
