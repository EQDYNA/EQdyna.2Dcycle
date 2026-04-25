subroutine loadGmshMesh
! Reads externally generated gmsh mesh (C_mesh=3) from:
!   meshGeneralInfo.txt, vert.txt, fac.txt, nsmp.txt
! Coordinates in those files are in km; converted to m here.
! Node indices in fac.txt and nsmp.txt are 0-based; converted to 1-based here.
! Padding rows in nsmp.txt are (0,0) and remain (0,0) in nsmp0.

use globalvar
implicit none

integer(kind=4) :: i, j, itmp1, itmp2, mxfn
real(kind=dp) :: km2m
logical :: fexist
character(len=60) :: fname

km2m = 1.0d3

! ---- meshGeneralInfo.txt ----
fname = 'meshGeneralInfo.txt'
inquire(file=trim(fname), exist=fexist)
if (.not. fexist) then
    write(*,*) '=                                                                   ='
    write(*,*) '= ERROR: meshGeneralInfo.txt not found.                             ='
    write(*,*) '=   Run "python3 meshgen.py" first, then re-run case.setup.         ='
    write(*,*) '=                                                                   ='
    stop
endif
open(unit=20, file=trim(fname), form='formatted', status='old')
    read(20,*) numnp, numel
    read(20,*) (nfnode(i), i=1, ntotft)
    read(20,*) xmin, xmax, ymin, ymax
close(20)
xmin = xmin * km2m
xmax = xmax * km2m
ymin = ymin * km2m
ymax = ymax * km2m

mxfn = maxval(nfnode)
write(*,*) '=     meshGeneralInfo.txt loaded                                    ='
write(*,'(X,A,36X,i8,A,i8)') '=', numnp, ' nodes ', numel
write(*,'(X,A,36X,10i6)') '=', (nfnode(i), i=1, ntotft)

! ---- vert.txt ----
fname = 'vert.txt'
inquire(file=trim(fname), exist=fexist)
if (.not. fexist) then
    write(*,*) '= ERROR: vert.txt not found. Run "python3 meshgen.py" first.'
    stop
endif
open(unit=21, file=trim(fname), form='formatted', status='old')
    do i = 1, numnp
        read(21,*) xnode0(1,i), xnode0(2,i)
        xnode0(1,i) = xnode0(1,i) * km2m
        xnode0(2,i) = xnode0(2,i) * km2m
    enddo
close(21)

! ---- fac.txt (0-indexed → 1-indexed) ----
fname = 'fac.txt'
inquire(file=trim(fname), exist=fexist)
if (.not. fexist) then
    write(*,*) '= ERROR: fac.txt not found. Run "python3 meshgen.py" first.'
    stop
endif
open(unit=22, file=trim(fname), form='formatted', status='old')
    do j = 1, numel
        read(22,*) ien0(1,j), ien0(2,j), ien0(3,j), ien0(4,j)
        ien0(1,j) = ien0(1,j) + 1
        ien0(2,j) = ien0(2,j) + 1
        ien0(3,j) = ien0(3,j) + 1
        ien0(4,j) = ien0(4,j) + 1
    enddo
close(22)

! ---- nsmp.txt (0-indexed → 1-indexed; (0,0) padding rows kept as 0) ----
fname = 'nsmp.txt'
inquire(file=trim(fname), exist=fexist)
if (.not. fexist) then
    write(*,*) '= ERROR: nsmp.txt not found. Run "python3 meshgen.py" first.'
    stop
endif
open(unit=23, file=trim(fname), form='formatted', status='old')
    do i = 1, mxfn * ntotft
        read(23,*) itmp1, itmp2
        if (itmp1 == 0 .and. itmp2 == 0) then
            nsmp0(1,i) = 0
            nsmp0(2,i) = 0
        else
            nsmp0(1,i) = itmp1 + 1
            nsmp0(2,i) = itmp2 + 1
        endif
    enddo
close(23)

write(*,*) '=     Gmsh mesh loaded successfully                                 ='

end subroutine loadGmshMesh
