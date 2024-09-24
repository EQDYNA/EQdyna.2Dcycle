subroutine loadGmshMesh
    use globalvar
    implicit none
    integer (kind = 4) :: i, j
    real (kind = dp) :: tmp1, tmp2, tmp3
    
    open(14, file = 'meshGeneralInfo.txt', form = 'formatted', status='old')
        read(14,*) numnp, numel
        read(14,*) (nfnode(i), i=1,ntotft)
        read(14,*) xmin, xmax, ymin, ymax
    close(14)
    totftnode = sum(nfnode)
    maxftnode = maxval(nfnode)
    xmin = xmin*1.d3
    xmax = xmax*1.d3
    ymin = ymin*1.d3
    ymax = ymax*1.d3
    
    if (debug==1) then 
        write(*,*) 'Num of nodes, cells are', numnp, numel
        write(*,*) 'Num of fault nodes are', (nfnode(i),i=1,ntotft)
        write(*,*) 'Total num of ft nodes is', totftnode
        write(*,*) 'Max num per ft node is', maxftnode
        write(*,*) 'Model ranges are', xmin, xmax, ymin, ymax
    endif  
    
    open(11,file='vert.txt', form = 'formatted', status='old')
        do i = 1,numnp
            read(11,*) xnode0(1,i), xnode0(2,i)
            if (i>5034 .and. i<5050) write(*,*) xnode0(1,i), xnode0(2,i), i
        enddo
    close(11)
    xnode0 = xnode0*1.0d3
    
    open(12,file='fac.txt', form = 'formatted', status='old')
        do i = 1, numel
            read(12,*) ien0(1,i), ien0(2,i), ien0(3,i), ien0(4,i)
            ien0(1:4,i) = ien0(1:4,i)+1
            !if (i==3663 .or. i==3662) write(*,*) ien0(1:4,i), i
        enddo
    close(12)
    
    open(14, file = 'nsmp.txt', form = 'formatted', status ='old')
        do i = 1,maxftnode*3
            read(14,*) nsmp0(1,i), nsmp0(2,i)
            nsmp0(1,i) = nsmp0(1,i)+1 ! adjust from Python index 0 rule
            nsmp0(2,i) = nsmp0(2,i)+1
        enddo 
    close(14)
    
    open(14, file = 'nsmpGeoPhys.txt', form = 'formatted', status ='old')
        do i = 1,maxftnode*3
            read(14,*) tmp1, tmp2, tmp3
            nsmpnv(1,i) = -tmp2
            nsmpnv(2,i) = tmp1 
            nsmpnv(3,i) = tmp3*1.0d3
        enddo
    close(14) 
    write(*,*) 'nsmpnv is ', nsmpnv(3,200)
    ! open(14, file = 'nsmpnv.txt', form = 'formatted', status = 'unknown')
    ! do j = 1, ntotft
        ! write(14, '(3e21.14)') (-dirvec(2,i,j), dirvec(1,i,j), dirvec(3,i,j), i = 1, ntemp)
    ! enddo
    ! close(14)
end subroutine loadGmshMesh

subroutine meshgen1
!... program to generate mesh for 2D models.
! This is a bend-stepover fault model with 2 faults, a simpified model for Altyn Tagh in China.
! This model is most complex 2D model so far in my study, as it combines bend and stepover.
! The scheme here is to use regular mesh as a base and then use degeneration technique
! to conform fault geometry. Triangular elements are introduced by degeneration
! from quadrilateral elements near faults. This sounds relatively simple to deal with
! complex fault geometry.
! Author: Benchun Duan (B.D.) March 2009.
!
!... revised to create smooth faults using natural cubic spline functions. B.D. 11/23/11
!
  use globalvar
  implicit none
!
!*********************************************************************************************************** 
! !!! To users: you may only need to modify within this block confined by 2 * lines!!!
!
!...MODEL GEOMETRY PARAMETER that most likely need to adjust by users for a specific model.
!    For example, element size, model range, fault corner coordinates, etc. B.D. 6/15/08
  !......faults 1 & 2 control points number given here; coordinates read in from files.
  real(kind = dp),allocatable,dimension(:,:,:) :: fxyr   !min, max x, y coordinates for faults
!**********************************************************************************************************
  real(kind = dp), allocatable,dimension(:)::xsp,ysp,hsp,asp,bsp,csp,dsp,esp  !for spline calculations
  real(kind = dp),allocatable,dimension(:,:)::f,xx,yy,hh
  character(len=155555)::filein
  logical,allocatable,dimension(:)::ny4flt,ny4dege,ny4abv,ny4abv2,ny4abv3,ny4abv4
  integer (kind=4) :: i,j,k,ix,iy,nx,ny,nxuni,nyuni,edgex1,edgex2,edgey1,edgey2,&
      ntemp,ynplot, np=100000, nft, nnode, nelement, ntag
  integer (kind=4),allocatable,dimension(:) :: line1,line2
  real (kind = dp) :: tol,xstep,ystep,xcoor,ycoor,temp1,temp2,temp3,temp4,temp5,fxmin,fxmax,fymin,fymax,pp(2,4)
  real (kind = dp),allocatable,dimension(:) :: xline,yline,ylinetmp
  real (kind = dp) :: xc(3)
  integer (kind = 4),allocatable,dimension(:,:,:)::nsmp
  real (kind = dp), allocatable, dimension(:,:,:)::dirvec  
!**************************************************************************
!!!! set number of faults.
  nft = ntotft !number of faults in model
!**************************************************************************
  allocate(fxyr(2,2,nft),ny4flt(nft),ny4dege(nft),ny4abv(nft),ny4abv2(nft),ny4abv3(nft),ny4abv4(nft))
  allocate(f(ftn,nft),xx(ftn,nft),yy(ftn,nft),hh(ftn,nft))
  f=0.0d0
  xx=0.0d0
  yy=0.0d0
  hh=0.0d0
!...calculate spline coefficients for the faults.

  do i=1,nft
    allocate(xsp(ftcn(i)),ysp(ftcn(i)),hsp(ftcn(i)-1),asp(ftcn(i)-2),bsp(ftcn(i)-2),csp(ftcn(i)-2),&
             dsp(ftcn(i)-2),esp(ftcn(i)-2))
    !...read control points (x,y) coordinate from predefined files. x-coor should be
    !   from min to max in these files: in sequence for spline! B.D. 11/23/11
      ! if (ntotft>3) then 
            ! write(*,*) 'Too many faults, in meshgen1 line 55'
            !  
      ! endif     
    if(i==1) then
      filein='x1_1.txt'
    elseif(i==2) then
      filein='x2_1.txt'
    elseif (i == 3) then
    filein = 'x3_1.txt'
    elseif (i == 4) then
    filein = 'x4.txt'
    endif
    open(15,file=filein,status='old')
    read(15,*) (xsp(j),ysp(j),j=1,ftcn(i))
    close(15)
    ! xsp = xsp * 1.0d3
    ! ysp = ysp * 1.0d3
    do j=1,ftcn(i)  !save for later interpolation use.
      xx(j,i) = xsp(j)
      yy(j,i) = ysp(j)
    enddo
    fxyr(1,1,i)=minval(xsp)
    fxyr(2,1,i)=maxval(xsp)
    fxyr(1,2,i)=minval(ysp)
    fxyr(2,2,i)=maxval(ysp)
    !...form a,b,c,d for natural cubic spline.
    do j=1,ftcn(i)-1
      hsp(j) = xsp(j+1) - xsp(j)
      hh(j,i) = hsp(j)  !save for later interpolation use.
    enddo
    do j=1,ftcn(i)-2
      asp(j) = hsp(j)
      if(j==1) asp(j)=0.0d0
      bsp(j) = 2.0d0 * (hsp(j) + hsp(j+1))
      csp(j) = hsp(j+1)
      if(j==ftcn(i)-2) csp(j)=0.0d0
      dsp(j) = 6.0d0 * ((ysp(j+2)-ysp(j+1))/hsp(j+1) - (ysp(j+1)-ysp(j))/hsp(j))
    enddo
    !...forward modify coefficients.
    csp(1) = csp(1)/bsp(1)
    dsp(1) = dsp(1)/bsp(1)
    do j=2,ftcn(i)-2
      temp1=bsp(j)-csp(j-1)*asp(j)
      csp(j) = csp(j)/temp1
      dsp(j) = (dsp(j) - dsp(j-1)*asp(j))/temp1
    enddo
    !...backward substitute for solution
    esp(ftcn(i)-2) = dsp(ftcn(i)-2)
    do j = ftcn(i)-3,1,-1
      esp(j) = dsp(j) - csp(j) * esp(j+1)
    enddo
    !...spline coefficients for ftcn(i) data points
    f(1,i)= 0.0d0
    f(ftcn(i),i) = 0.0d0
    do j=2,ftcn(i)-1
      f(j,i) = esp(j-1)
    enddo
    deallocate(xsp,ysp,hsp,asp,bsp,csp,dsp,esp)
  enddo ! enddo for i=1,nft
    
    if (debug == 1) write(*,*) 'Finish ft loop and calculate derived parameters'
    
    cellmin = dxy/2.0d0 !min grid size for dt for filing.f90
    tol = dxy/100.d0    !tolerance for comparison of equal
    
    if (ntotft>4) then 
        write(*,*) 'Too many faults, in meshgen1 line 114'
         
    endif 
    fxmin=min(fxyr(1,1,1), fxyr(1,1,2), fxyr(1,1,3), fxyr(1,1,4))
    fymin=min(fxyr(1,2,1), fxyr(1,2,2), fxyr(1,2,3), fxyr(1,2,4))
    fxmax=max(fxyr(2,1,1), fxyr(2,1,2), fxyr(2,1,3), fxyr(2,1,4))
    fymax=max(fxyr(2,2,1), fxyr(2,2,2), fxyr(2,2,3), fxyr(2,2,4))

    xmax = fxmax + term/2.0d0*vp
    xmin = fxmin - term/2.0d0*vp
    ymax = fymax + term/2.0d0*vp
    ymin = fymin - term/2.0d0*vp
  
    if (debug==1) write(*,*) 'Domain ranges from ', xmax, xmin, ymax, ymin
    
    nxuni = (fxmax-fxmin)/dxy
    xstep = dxy
    xcoor = fxmin
    do ix = 1,np
        xstep = xstep * rat
        xcoor = xcoor - xstep
        if(xcoor < xmin) exit
    enddo
    edgex1 = ix+1
    xstep = dxy
    xcoor = fxmax
    do ix = 1,np
        xstep = xstep * rat
        xcoor = xcoor + xstep
        if(xcoor > xmax) exit
    enddo
    edgex2 = ix+1 
    nx = edgex1 + nxuni + edgex2 +1
    !...determine num of nodes along y
    nyuni = (fymax-fymin+2.0d0*yext)/dxy
    ystep = dxy
    ycoor = fymin-yext
    do iy = 1,np
        ystep = ystep * rat
        ycoor = ycoor - ystep
        if(ycoor < ymin) exit
    enddo
    edgey1 = iy+1
    ystep = dxy
    ycoor = fymax+yext
    do iy = 1,np
        ystep = ystep * rat
        ycoor = ycoor + ystep
        if(ycoor > ymax) exit
    enddo
    edgey2 = iy+1 
    ny = edgey1 + nyuni + edgey2 +1
    if (debug==1) write(*,*) 'number of node along x, and y-coor', nx, ny
    !  
    ntemp = max(nxuni,nyuni) + 1
    write(*,*) 'maxftnode = ', ntemp
    !
    allocate(line1(ny+nft),line2(ny+nft),xline(nx),yline(ny),ylinetmp(ny),&
    nsmp(2,ntemp,nft),dirvec(3,ntemp,nft))    
    !...predetermine x-coor
    xline(edgex1+1) = fxmin
    xstep = dxy
    do ix = edgex1,1,-1
    xline(ix) = xline(ix+1) - xstep
    xstep = xstep * rat
    enddo
    xmin = xline(1)
    do ix = edgex1+2,edgex1+nxuni+1
    xline(ix) = xline(ix-1) + dxy
    enddo
    xstep = dxy
    do ix = edgex1+nxuni+2,nx
    xline(ix) = xline(ix-1) + xstep
    xstep = xstep * rat
    enddo
    xmax = xline(nx)
    !...predetermine y-coor
    yline(edgey1+1) = fymin-yext
    ystep = dxy
    do iy = edgey1,1,-1
    yline(iy) = yline(iy+1) - ystep
    ystep = ystep * rat
    enddo
    ymin = yline(1)
    do iy = edgey1+2,edgey1+nyuni+1
    yline(iy) = yline(iy-1) + dxy
    enddo
    ystep = dxy
    do iy = edgey1+nyuni+2,ny
    yline(iy) = yline(iy-1) + ystep
    ystep = ystep * rat
    enddo
    ymax = yline(ny)
           
    !...prepare for digitizing
    nnode = 0
    nelement = 0
    nfnode = 0
    line1 = 0
    line2 = 0
    nsmp = 0
    dirvec = 0.0d0
    ny4flt = .false.
    ny4dege = .false.
    ny4abv = .false.
    ny4abv2 = .false.
    ny4abv3 = .false.
    ny4abv4 = .false.  
    ylinetmp = yline 
  
    if (debug==1) write(*,*) 'digitize along constant x line (perpendicular to main fault)'
    do ix = 1, nx
        ylinetmp = yline
        xcoor = xline(ix)
        !write(*,*) 'ix,xcoor= ',ix,xcoor
        do iy = 1, ny
            ycoor = ylinetmp(iy)
                !....create nodes & elements together for general cases.
            nnode = nnode + 1
            line2(iy) = nnode
            xnode0(1,nnode) = xcoor
            xnode0(2,nnode) = ycoor
            
            if(ix >= 2 .and. iy >= 2) then
                nelement = nelement + 1
                ien0(1,nelement) = line2(iy)
                ien0(2,nelement) = line1(iy)
                ien0(3,nelement) = line1(iy-1)
                ien0(4,nelement) = line2(iy-1)          
                
                ! !....special treatments for elements above faults.
                ! do k=1,nft
                    ! if(ny4flt(k)) then
                        ! if(nfnode(k)==1) then
                            ! ien0(4,nelement) = line2(ny+k) !Assign master nodes.
                        ! else
                            ! ien0(3,nelement) = line1(ny+k) !Assign master nodes.
                            ! ien0(4,nelement) = line2(ny+k)
                            ! if(ny4abv(k)) then    !special quadrilateral above the fault
                                ! ien0(2,nelement) = line1(iy+1)
                                ! ien0(3,nelement) = line1(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                ! nelement = nelement + 1 
                                ! ien0(1,nelement) = line1(ny+k)
                                ! ien0(2,nelement) = line2(ny+k) !Assign master nodes.
                                ! ien0(3,nelement) = line2(iy)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                ! ny4dege(k) = .true.
                                ! ny4abv(k) = .false.
                            ! elseif(ny4abv2(k)) then !special quadrilateral and degeneration above fault for case 2
                                ! ien0(2,nelement) = line1(iy-1)
                                ! ien0(3,nelement) = line2(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                ! nelement = nelement + 1 
                                ! ien0(1,nelement) = line2(ny+k) !Assign master nodes.
                                ! ien0(2,nelement) = line1(iy-1) 
                                ! ien0(3,nelement) = line1(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)    
                                ! nelement = nelement + 1
                                ! ien0(1,nelement) = line2(iy)
                                ! ien0(2,nelement) = line1(iy)
                                ! ien0(3,nelement) = line1(iy-1)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                ! ny4abv2(k) = .false.
                            ! elseif(ny4abv3(k)) then    
                            ! ! ien0(1,nelement) = line2(iy)
                            ! ! ien0(2,nelement) = line1(iy)
                            ! ! ien0(3,nelement) = line1(ny+k)
                            ! ! ien0(4,nelement) = line2(ny+k)     
                                ! ien0(1,nelement) = line2(iy)
                                ! ien0(2,nelement) = line1(iy)
                                ! ien0(3,nelement) = line1(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                
                                ! nelement = nelement + 1
                                ! ien0(1,nelement) = line2(iy)
                                ! ien0(2,nelement) = line1(ny+k)
                                ! ien0(3,nelement) = line2(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)                
                                ! ny4abv3(k) = .false.
                            ! elseif(ny4abv4(k)) then    
                                ! ien0(1,nelement) = line2(iy)
                                ! ien0(2,nelement) = line1(iy)
                                ! ien0(3,nelement) = line2(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)
                                
                                ! nelement = nelement + 1
                                ! ien0(1,nelement) = line2(ny+k)
                                ! ien0(2,nelement) = line1(iy)
                                ! ien0(3,nelement) = line1(ny+k)
                                ! ien0(4,nelement) = ien0(3,nelement)    
                                ! ny4abv4(k) = .false.
                            ! endif
                        ! endif
                    ! endif
                ! enddo
            endif
    
            ! !if (debug==1) write(*,*) 'Loop over faults to see deciding if to degenerate an element'
            ! do k=1,nft
                ! if(.not.ny4flt(k) .and. ny4dege(k)) then
                    ! !the current element must be degenerated to a triangular
                    ! ien0(3,nelement) = line2(iy-1)
                    ! ien0(4,nelement) = ien0(3,nelement)
                    ! ny4dege(k) = .false.
                ! endif
            ! enddo
    !if (debug==1) write(*,*) 'conform to fault geometry, create split nodes, and fault-related elements.'
    
    do k=1,nft  ! loop over fault #k
        ny4flt(k)=.false.
        if(xcoor>=fxyr(1,1,k).and.xcoor<=fxyr(2,1,k) .and. &
            ycoor>=(fxyr(1,2,k)-yext).and.ycoor<=(fxyr(2,2,k)+yext)) then
            write(*,*) 'ft id is ',k 
            write(*,*) fxyr(1,1,k), fxyr(2,1,k),fxyr(1,2,k), fxyr(2,2,k)
            write(*,*) xcoor/1e3, ycoor/1e3
            do i=1,ftcn(k)-1
                write(*,*) ftcn(k), i
                if(xcoor>=(xx(i,k)-tol) .and. xcoor<=(xx(i+1,k)+tol)) then
                    if (i==1) then 
                        pp(1,2) = xx(i,k) 
                        pp(1,3) = xx(i+1,k)
                        pp(1,4) = xx(i+2,k)
                        pp(1,1) = xx(i,k) - (xx(i+1,k) - xx(i,k))
                        pp(2,2) = yy(i,k) 
                        pp(2,3) = yy(i+1,k)
                        pp(2,4) = yy(i+2,k)
                        pp(2,1) = yy(i,k) - (yy(i+1,k) - yy(i,k))
                    elseif (i==ftcn(k)-1) then 
                        pp(1,1) = xx(i-1,k)
                        pp(1,2) = xx(i,k) 
                        pp(1,3) = xx(i+1,k)
                        pp(1,4) = xx(i+1,k) + (xx(i+1,k) - xx(i,k))
                        pp(2,1) = yy(i-1,k)
                        pp(2,2) = yy(i,k) 
                        pp(2,3) = yy(i+1,k)
                        pp(2,4) = yy(i+1,k) + (xx(i+1,k) - xx(i,k))            
                    else 
                        pp(1,1) = xx(i-1,k)
                        pp(1,2) = xx(i,k) 
                        pp(1,3) = xx(i+1,k)
                        pp(1,4) = xx(i+2,k) 
                        pp(2,1) = yy(i-1,k)
                        pp(2,2) = yy(i,k) 
                        pp(2,3) = yy(i+1,k)
                        pp(2,4) = yy(i+2,k)                 
                    endif 
                    
                    call catmullrom(pp, xcoor, temp1, 0.5d0)
                    write(*,*) pp, xcoor, temp1
              ! temp1 = (f(i+1,k)*(xcoor-xx(i,k))**3+f(i,k)*(xx(i+1,k)-xcoor)**3)/(6.0d0*hh(i,k)) &
                     ! +(yy(i+1,k)/hh(i,k)-hh(i,k)*f(i+1,k)/6.0d0)*(xcoor-xx(i,k)) &
                     ! +(yy(i,k)/hh(i,k)-hh(i,k)*f(i,k)/6.0d0)*(xx(i+1,k)-xcoor)
                    if(temp1>yline(iy-1) .and. temp1<=yline(iy)) then
                    temp2 = temp1 - yline(iy-1)
                    temp3 = yline(iy) - temp1
                    if(temp3 <= temp2) then
                      ycoor = temp1
                      temp4 = yline(iy+1)
                      ylinetmp(iy+1) = temp4 -2.0d0*temp3/3.0d0
                      temp4 = yline(iy+2)
                      ylinetmp(iy+2) = temp4 -1.0d0*temp3/3.0d0
                      xnode0(2,nnode-1) = xnode0(2,nnode-1) -2.0d0*temp3/3.0d0
                      xnode0(2,nnode-2) = xnode0(2,nnode-2) -1.0d0*temp3/3.0d0
                      ny4flt(k) = .true.
                    endif
                    elseif(temp1>yline(iy) .and. temp1<=yline(iy+1)) then
                    temp2 = temp1 - yline(iy)
                    temp3 = yline(iy+1) - temp1
                    if(temp2 <= temp3) then
                      ycoor = temp1
                      temp4 = yline(iy+1)
                      ylinetmp(iy+1) = temp4 +2.0d0*temp2/3.0d0
                      temp4 = yline(iy+2)
                      ylinetmp(iy+2) = temp4 +1.0d0*temp2/3.0d0
                      xnode0(2,nnode-1) = xnode0(2,nnode-1) +2.0d0*temp2/3.0d0
                      xnode0(2,nnode-2) = xnode0(2,nnode-2) +1.0d0*temp2/3.0d0
                          ny4flt(k) = .true.
                    endif
                    endif
                    
                  if(ny4flt(k)) then
                    xnode0(2,nnode) = ycoor
                    nfnode(k) = nfnode(k) + 1
                    nsmp(1,nfnode(k),k) = nnode    !slave node
                    nnode = nnode + 1
                    nsmp(2,nfnode(k),k) = nnode    !master node
                    line2(ny+k) = nnode
                    xnode0(1,nnode) = xcoor
                    xnode0(2,nnode) = ycoor
                    if(nfnode(k)>1) then    !1st fault node, no need for this
                          ny4abv(k) = .false.
                          ny4abv2(k) = .false.
                            ny4abv3(k) = .false.
                          ny4abv4(k) = .false.    
                      if((xnode0(2,nsmp(1,nfnode(k)-1,k))-xnode0(2,line1(iy)))>tol) then !for down slope
                          !use ycoor to judge fault trace
                        !first, the current element is still quadrilateral, but with special nodes.
                     ien0(1,nelement) = line2(iy)    !this is same as usual
                      ien0(2,nelement) = line1(iy+1)    !the following two are special
                      ien0(3,nelement) = line1(iy)
                    ien0(4,nelement) = ien0(3,nelement)
                    nelement = nelement + 1
                    ien0(1,nelement) = line2(iy)    
                      ien0(2,nelement) = line1(iy)    
                      ien0(3,nelement) = line2(iy-1)
                      ien0(4,nelement) = ien0(3,nelement)    
                      !then,one more element is created by degeration, which is a triangular element.
                        nelement = nelement + 1
                        ien0(1,nelement) = line1(iy)
                        ien0(2,nelement) = line1(iy-1)
                        ien0(3,nelement) = line2(iy-1)
                        ien0(4,nelement) = ien0(3,nelement)
                            ny4abv(k) = .true.    !for special treatment above the fault
                          elseif((xnode0(2,line1(iy))-xnode0(2,nsmp(1,nfnode(k)-1,k)))>tol) then ! for up slope
                            !first, previous element degenerated to triangular
                            ien0(1,nelement-1) = line2(iy-1)
                            ien0(2,nelement-1) = line1(iy-2)
                            ien0(3,nelement-1) = line2(iy-2)
                            ien0(4,nelement-1) = ien0(3,nelement-1)
                            !then, current special quadrilateral
                            ien0(1,nelement) = line2(iy)
                            ien0(2,nelement) = line1(iy-1)
                            ien0(3,nelement) = line2(iy-1)
                            ien0(4,nelement) = ien0(3,nelement)
                            nelement = nelement + 1
                            ien0(1,nelement) = line2(iy-1)
                            ien0(2,nelement) = line1(iy-1)
                            ien0(3,nelement) = line1(iy-2)
                            ien0(4,nelement) = ien0(3,nelement)
                            ny4abv2(k) = .true. !for special treatment above fault: case 2 up slope
                          else
                            if (xnode0(2,line1(iy))>=ycoor) then 
                            ien0(1,nelement) = line2(iy)
                            ien0(2,nelement) = line1(iy)
                            ien0(3,nelement) = line1(iy-1)
                            ien0(4,nelement) = ien0(3,nelement)
                            
                            nelement = nelement + 1 
                            ien0(1,nelement) = line1(iy-1)
                            ien0(2,nelement) = line2(iy-1)
                            ien0(3,nelement) = line2(iy)
                            ien0(4,nelement) = ien0(3,nelement)    
                            ny4abv3(k) = .true. !for special treatment above fault: general fault related elements.                    
                            else
                            ien0(1,nelement) = line2(iy)
                            ien0(2,nelement) = line1(iy)
                            ien0(3,nelement) = line2(iy-1)
                            ien0(4,nelement) = ien0(3,nelement)                    

                            nelement = nelement + 1 
                            ien0(1,nelement) = line2(iy-1)
                            ien0(2,nelement) = line1(iy)
                            ien0(3,nelement) = line1(iy-1)
                            ien0(4,nelement) = ien0(3,nelement)    
                            ny4abv4(k) = .true. !for special treatment above fault: general fault related elements.
                            endif
                          endif                
                    endif
                    !...cal unit normal and tangent at each split-nodes pair: using spline
                          temp2=(f(i,k)*(xx(i+1,k)-xcoor)**2-f(i+1,k)*(xcoor-xx(i,k))**2)/(2*hh(i,k)) &
                               +(6.0d0*(yy(i,k)-yy(i+1,k))+hh(i,k)*hh(i,k)*(f(i+1,k)-f(i,k)))/(6.0d0*hh(i,k))
                          temp3=sqrt(temp2*temp2+1)
                          dirvec(1,nfnode(k),k)=1/temp3       !tx: x-comp unit tangent
                          dirvec(2,nfnode(k),k)=-temp2/temp3  !ty: y-comp unit tangent
                        if(nfnode(k)>2) then ! legnth associated with the previous fault node: linear segment
                      temp2 = xnode0(1,nsmp(1,nfnode(k)-1,k)) - xnode0(1,nsmp(1,nfnode(k)-2,k))
                      temp3 = xnode0(2,nsmp(1,nfnode(k)-1,k)) - xnode0(2,nsmp(1,nfnode(k)-2,k))
                      temp4 = sqrt(temp2*temp2+temp3*temp3)
                      temp2 = xnode0(1,nsmp(1,nfnode(k),k)) - xnode0(1,nsmp(1,nfnode(k)-1,k))
                      temp3 = xnode0(2,nsmp(1,nfnode(k),k)) - xnode0(2,nsmp(1,nfnode(k)-1,k))
                      temp5 = sqrt(temp2*temp2+temp3*temp3)
                      dirvec(3,nfnode(k)-1,k) = 0.5d0 * (temp4 + temp5)    !length associated with node
                    endif
                    exit
                  endif           
                endif
              enddo
            endif
        enddo  !enddo for k=1,nft
    !
    enddo    !iy
    line1 = line2
  enddo        !ix
    write(*,*) nnode, nelement
    if (debug==1) write(*,*) 'Finished ix,iy loop. Special calculation for fault split node length at two ends of each fault'

        do i=1,nft
        write(*,*) nfnode
    do j=1,2
      if(j==1) then
        temp2 = xnode0(1,nsmp(1,2,i)) - xnode0(1,nsmp(1,1,i))
        temp3 = xnode0(2,nsmp(1,2,i)) - xnode0(2,nsmp(1,1,i))
        temp4 = sqrt(temp2*temp2+temp3*temp3)
        k = 1
        write(*,*) j,i,temp2, temp3, temp4
      else
        temp2 = xnode0(1,nsmp(1,nfnode(i),i)) - xnode0(1,nsmp(1,nfnode(i)-1,i))
        temp3 = xnode0(2,nsmp(1,nfnode(i),i)) - xnode0(2,nsmp(1,nfnode(i)-1,i))
        temp4 = sqrt(temp2*temp2+temp3*temp3)
        k = nfnode(i)
        write(*,*) j,i,temp2,temp3,temp4,nfnode(i)
      endif      
      dirvec(3,k,i) = 0.5d0 * temp4    !txy: 'length' of the node
    enddo
  enddo

    write(*,*) 'Entire model region has:'
    write(*,*) nnode,' nodes; ',nelement,' elements.'
    if (plotmesh == 1) then 
        write(*,*) '=     Wrting out mesh data for matlab to visualize ...              ='
        open(11,file='vert.txt',status='unknown')
        write(11,'(1x,2e18.7e4)') ((xnode0(j,i),j=1,2),i=1,nnode)
        close(11)
        
        open(12,file='fac.txt',status='unknown')
        write(12,'(1x,4i15)') ((ien0(j,i),j=1,4),i=1,nelement)
        close(12)

        open(13,file='Mesh_general_info.txt',status='unknown')
        write(13,'( 2i15)') nnode,nelement
        write(13,'( 4i15)') (nfnode(k),k=1,nft)
        ! write(13,'( a60)') 'fault nodes'' x, y- xoor'
        ! write(13,'( a20)') 'fault #1'
        ! write(13,'( i10,2f12.1)') (nsmp(1,i,1),(xnode0(j,nsmp(1,i,1)),j=1,2),i=1,nfnode(1))
        ! write(13,'( a20)') 'fault #2'
        ! write(13,'( i10,2f12.1)') (nsmp(1,i,2),(xnode0(j,nsmp(1,i,2)),j=1,2),i=1,nfnode(2))
        close(13)
    endif

  ntemp = maxval(nfnode)
  ntag = 0
  do i = 1, nft
    do j = 1, ntemp
        ntag = ntag + 1
        nsmp0(1,ntag) = nsmp(1,j,i)
        nsmp0(2,ntag) = nsmp(2,j,i)
    enddo 
  enddo 
  ntag = 0
  do i = 1, nft
    do j = 1, ntemp
        ntag = ntag + 1
        nsmpnv(1,ntag) = -dirvec(2,j,i)
        nsmpnv(2,ntag) = dirvec(1,j,i)
        nsmpnv(3,ntag) = dirvec(3,j,i)
    enddo 
  enddo  
    if (plotmesh == 1) then  
        open(14, file = 'nsmp.txt', form = 'formatted', status ='unknown')
        do j = 1, nft 
             write(14,'(2i15)') (nsmp(1,i,j), nsmp(2,i,j), i = 1, ntemp)
        enddo
        close(14)

        open(14, file = 'nsmpnv.txt', form = 'formatted', status = 'unknown')
        do j = 1, nft
            write(14, '(3e21.14)') (-dirvec(2,i,j), dirvec(1,i,j), dirvec(3,i,j), i = 1, ntemp)
        enddo
        close(14)
    endif
  numnp = nnode
  numel = nelement

end subroutine meshgen1



              
      
      
