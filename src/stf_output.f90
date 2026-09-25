! stf_output.f90 -- per-event source-time-function (STF) output.
!
! One file per run, stf.bin<icstart>, Fortran stream access (no record
! markers), native little-endian. It holds a file header with the static fault
! node table, followed by one record per earthquake cycle. Read it with
! scripts/stf_read.py; the layout below is the contract that reader implements.
!
! FILE HEADER (written once, when the file is new or empty)
!   char(8)   magic            'EQDYNSTF'
!   int32     version          1
!   int32     icstart          first cycle of this run segment
!   int32     ntotft           number of faults
!   int32     totftnode        number of fault nodes (the node table length)
!   int32     nvar             number of time-series variables (7)
!   real64    dt_solver        solver time step, s
!   int32     stfEvery         output decimation (dt_out = stfEvery*dt_solver)
!   real64    stfVmin          node-selection threshold on peak slip rate, m/s
!   real64    tSeqStartYr      sequence time at the start of this segment, yr
!   int32     originKnown      1 if tSeqStartYr is the true time since model
!                              start; 0 if a restart could not recover it
!   char(16)  varname(nvar)
!   char(16)  varunit(nvar)
!   node table, totftnode entries in the model's sequential node order n:
!     int32 fault, int32 local_index, real64 x, y (m),
!     real64 tx, ty (unit fault tangent), real64 length (m, tributary length)
!
! EVENT RECORD (one per cycle)
!   char(4)   'EVNT'
!   int64     nbytes           bytes that follow, up to and including the trailer
!   int32     eqid             cycle number, matches catalog.csv eqId
!   real64    time_yr          time of this event in the sequence, yr since
!                              model start (sum of intervals incl. this one)
!   real64    interval_yr      interseismic interval that preceded this event
!   int32     nuc_node         nucleation node, sequential index into the table
!   real64    nuc_x, nuc_y     nucleation node coordinates, m
!   real64    t0               time of the first sample, s (= dt_solver: the
!                              solver advances time before it first samples)
!   real64    dt_out           sample interval, s; sample k (0-based) is at
!                              t0 + k*dt_out
!   real64    t_end            dynamic time at which the event was stopped, s
!   int32     nt               samples per series
!   int32     nout             number of nodes written
!   int32     node(nout)       sequential node indices (1-based) into the table
!   real64    rupture_time(nout)   s, 1000 = never ruptured
!   real64    peak_slip_rate(nout) m/s
!   real32    data(nt, nout, nvar) -- time fastest, then node, then variable
!   int32     eqid             trailer, repeated for integrity
!
! Variables: 1 slip_rate |m/s|, 2 slip_rate_t signed along the tangent,
! 3 slip |m|, 4 slip_t signed along the tangent, 5 shear_stress Pa,
! 6 normal_stress Pa (negative = compression), 7 friction coefficient.
!
! Event length is not fixed: an event stops when the peak slip rate on the
! fault drops below 1 mm/s after 5 s (faulting.f90), capped at `term`. Each
! record therefore carries its own nt and t_end.

subroutine write_stf(ic)
	use globalvar
	implicit none
	integer (kind=4), intent(in) :: ic
	character (len=30) :: mm
	character (len=64) :: fname
	character (len=16) :: vname(stfNvar), vunit(stfNvar)
	logical :: exists
	integer (kind=8) :: fsize, nbytes
	integer (kind=4) :: i, j, k, n, nout, iv, u
	integer (kind=4), allocatable :: idx(:)
	real (kind=dp) :: xnuc, ynuc

	vname = [character(len=16) :: 'slip_rate', 'slip_rate_t', 'slip', 'slip_t', &
	         'shear_stress', 'normal_stress', 'friction']
	vunit = [character(len=16) :: 'm/s', 'm/s', 'm', 'm', 'Pa', 'Pa', '1']

	write(mm,'(i6)') icstart
	fname = 'stf.bin'//trim(adjustl(mm))
	inquire(file=fname, exist=exists, size=fsize)

	open(newunit=u, file=fname, access='stream', form='unformatted', &
	     status='unknown', position='append')

	if ((.not. exists) .or. fsize <= 0) then
		write(u) 'EQDYNSTF'
		write(u) int(1,4), int(icstart,4), int(ntotft,4), int(totftnode,4), int(stfNvar,4)
		write(u) dt1, int(stfEvery,4), stfVmin, tSeqStartYr, int(stfOriginKnown,4)
		write(u) vname, vunit
		n = 0
		do j = 1, ntotft
			do i = 1, nfnode(j)
				k = (j-1)*maxftnode + i
				n = n + 1
				! tangent as used in faulting.f90: tx = -ny = nsmpnv(2), ty = nx = -nsmpnv(1)
				write(u) int(j,4), int(i,4), x(1,nsmp0(1,k)), x(2,nsmp0(1,k)), &
				         nsmpnv(2,k), -nsmpnv(1,k), nsmpnv(3,k)
			enddo
		enddo
	endif

	! Nodes worth writing: those that actually slipped during this event.
	nout = count(stfVpeak(1:totftnode) > stfVmin)
	allocate(idx(nout))
	n = 0
	do k = 1, int(totftnode,4)
		if (stfVpeak(k) > stfVmin) then
			n = n + 1
			idx(n) = k
		endif
	enddo

	! nucleation node coordinates, same sequential ordering as the table
	xnuc = 0.0d0; ynuc = 0.0d0
	n = 0
	do j = 1, ntotft
		do i = 1, nfnode(j)
			k = (j-1)*maxftnode + i
			n = n + 1
			if (n == loc) then
				xnuc = x(1,nsmp0(1,k)); ynuc = x(2,nsmp0(1,k))
			endif
		enddo
	enddo

	nbytes = 4_8 + 8_8 + 8_8 + 4_8 + 8_8 + 8_8 + 8_8 + 8_8 + 8_8 + 4_8 + 4_8 &
	       + 20_8*int(nout,8) + 4_8*int(stfNvar,8)*int(nout,8)*int(stfNt,8) + 4_8

	write(u) 'EVNT', nbytes
	write(u) int(ic,4), tSeqYr, tInterYr, int(loc,4), xnuc, ynuc, &
	         dt1, dt1*dble(stfEvery), timedyna, int(stfNt,4), int(nout,4)
	if (nout > 0) then
		write(u) idx
		write(u) (output4plot(5,idx(n)), n = 1, nout)
		write(u) (stfVpeak(idx(n)), n = 1, nout)
		write(u) (((stfBuf(k, idx(n), iv), k = 1, stfNt), n = 1, nout), iv = 1, stfNvar)
	endif
	write(u) int(ic,4)
	close(u)
	deallocate(idx)


end subroutine write_stf
