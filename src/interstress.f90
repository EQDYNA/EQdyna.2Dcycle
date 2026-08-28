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
integer (kind = 4)::iFtLoop, iNodeAcc
logical :: isFtEnd, isFtInterior
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
	if (quit .eqv. .TRUE.) then
		exit
	endif
	tinter = tinter + dtinter
	if (mod(ndt,1000) == 1) then
        write (*,*) 'Interseismic t = ', tinter/tcon, ' years'
	endif
	do i = 1, totftnode
		if (C_mesh == 2) then
			! C_mesh=2: read loading rate and angle from Rate_direction.txt.
			! rd(1,i) = local maximum shear strain rate γ(x), in rad/s.
			! rd(2,i) = local "angle of compression" φ(x) (= local fault
			!   strike − max shear-strain-rate direction), degrees, per-node.
			! Implements paper Liu et al. 2022 eqs (1)-(3):
			!   γ_τ = γ cos(2φ),  γ_n = γ sin(2φ),
			!   η = ant0·str/γ  (η inversely ∝ γ, "minimum viscosity" form)
			!   σ_τ → η γ_τ,  σ_n → η γ_n + σ^a  asymptotically.
			! Fallback uniform η = ant0 at zero-rate nodes (fault endpoints).
			theta = rd(2,ida(i))/180.0d0*pi
			if (theta >= 45.0d0/180.0d0*pi) theta = 45.0d0/180.0d0*pi
			if (rd(1,ida(i)) > 0.0d0) then
				ant = ant0 * str / rd(1,ida(i))
			else
				ant = ant0
			endif
			rs = rd(1,ida(i)) * cos(2.0d0*theta) * ant
			rn = -rd(1,ida(i)) * sin(2.0d0*theta) * ant
		else
			! C_mesh=3: per-node loading from nsmpGeoPhys.txt
			! cols (1-indexed): 1=tx,2=ty,3=len,4=ftType,5=ftDip,
			!   6=ftLoadMaxShear (γ, s^-1, per node)
			!   7=ftLoadAngle   (φ, deg, per node; -999 = auto from tangent)
			!   8=ftLoadWt      (loading weight; 450 = paper-faithful baseline)
			!   9=ftVis         (Pa·s, per node; = ant0·str/γ for paper-faithful)
			!
			! Paper-faithful encoding (saf.gmsh.lite, derived from
			! paper.saf.A/Rate_direction.txt):
			!   ftVis = ant0·str/γ → ant = ftVis (since ftLoadWt=450)
			!   rs   = γ·cos(2θ)·ant = str·cos(2θ)·ant0  (γ cancels;
			!          identical to C_mesh=2 paper interstress branch above).
			! Generic encoding (test.subei, test.gulang): ftLoadWt may
			! deviate from 450 to encode fractional loading; formula still
			! holds.  Unified formula:
			!   ant = ftVis · 450 / ftLoadWt
			!   rs  =  γ · (ftLoadWt/450) · cos(2θ) · ant
			!   rn  = -γ · (ftLoadWt/450) · sin(2θ) · ant
			if (nsmpgp(7,ida(i)) > -998.0d0) then
				theta = nsmpgp(7,ida(i))/180.0d0*pi
			else
				theta = atan(nsmpgp(2,ida(i))/nsmpgp(1,ida(i)))
			endif
			if (theta >= 45.0d0/180.0d0*pi) theta = 45.0d0/180.0d0*pi
			ant = nsmpgp(9,ida(i)) * 450.0d0 / nsmpgp(8,ida(i))
			rs =  nsmpgp(6,ida(i)) / 450.0d0 * nsmpgp(8,ida(i)) * cos(2.0d0*theta) * ant
			rn = -nsmpgp(6,ida(i)) / 450.0d0 * nsmpgp(8,ida(i)) * sin(2.0d0*theta) * ant
		endif
		! --- debug: dump per-node loading/stress on first interseismic step ---
		if (ndt == 1 .and. ic == 1) then
			if (i == 1) then
				open(2099, file='debug_interstress_init.txt', status='replace')
				write(2099,'(A)') '# node  ida  ss0(Pa)  ns0(Pa)  theta(deg)  ant(Pa.s)  rs(Pa)  rn(Pa)'
			endif
			write(2099,'(2i6, 6e18.8)') i, ida(i), ss0(i), ns0(i), &
				theta*180.0d0/pi, ant, rs, rn
			if (i == totftnode) close(2099)
		endif
		ns(i) = (ns0(i) - ambientnorm - rn) * exp(-tinter*amu/ant) + rn + ambientnorm
		shs(i) = (ss0(i) - rs) * exp(-tinter*amu/ant) + rs
		! Cap the effective normal stress the same way the dynamic solver does
		! (faulting.f90: if (tnrm > minnorm) tnrm = minnorm). A fault cannot
		! carry tension: it opens, and its strength is zero.
		!
		! This branch previously took abs(ns), so a node driven tensile by the
		! interseismic loading (rn > |ambientnorm|, which happens wherever the
		! angle of compression is sufficiently negative) was handed strength
		! PROPORTIONAL TO HOW TENSILE IT WAS -- the more it opened, the stronger
		! it got. Nodes near the sign change instead got strength ~0 and sat
		! permanently above failure, which no rupture could relieve.
		if (ns(i) > minnorm) ns(i) = minnorm
		strengthexcess(i) = (-ns(i)*fric_fs - shs(i))
		! Fault endpoints are excluded from nucleation: slip tapers to zero
		! at a tip, so the stress state there is not meaningful. Generic over
		! ntotft -- this was hardcoded to the first THREE faults, which
		! silently excluded nothing on faults 4+ and, worse, left them out of
		! the nucleation search below entirely.
		isFtEnd = .false.
		iNodeAcc = 0
		do iFtLoop = 1, ntotft
			if (i == iNodeAcc + 1 .or. i == iNodeAcc + nfnode(iFtLoop)) isFtEnd = .true.
			iNodeAcc = iNodeAcc + nfnode(iFtLoop)
		enddo
		if (isFtEnd) strengthexcess(i) = 100.0d6
	enddo
	minstrengthexcess = minval(strengthexcess)
	loc = minloc(strengthexcess, dim = 1)
	nucntag = 0
	do i = 1, totftnode
		! Interior nodes of every fault, generic over ntotft. Previously this
		! tested only faults 1-3, so with ntotft > 3 the remaining faults could
		! never nucleate -- gulang.gmsh.lite (5 faults) and any 7-fault model
		! were silently running with part of the system unable to host an
		! earthquake.
		isFtInterior = .false.
		iNodeAcc = 0
		do iFtLoop = 1, ntotft
			if (i > iNodeAcc + 1 .and. i < iNodeAcc + nfnode(iFtLoop)) isFtInterior = .true.
			iNodeAcc = iNodeAcc + nfnode(iFtLoop)
		enddo
		if (isFtInterior) then
			if (shs(i) > fric_fs * (-ns(i))) then 
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
flush(6)

open(2001,file='interval.txt'//m1, position='append')
	write(2001,'(1e21.14)') tinter/tcon
close(2001)


! if (nucntag > 1) then 
	! write(*,*) 'Many nucleation points', nucntag
	! write(*,*) (nuci(i), i = 1, nucntag) 
! endif
end subroutine interstress
      
