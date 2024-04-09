subroutine meshgen
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
  real(kind=8) :: dxy = 200.0d0	!element size along x, y in meter: basic parameter for resolution.
  real(kind=8) :: yext = 80.0d3  !distance in m to ensure large enough range in y beyond faults
  real(kind=8) :: xmin = -1169.0d3,xmax = 1077.0d3,ymin = -906.0d3,ymax = 979.0d3 !model boundary required
  real(kind=8) :: rat = 1.025d0	!ratio for enlarging element away from the fault
  !......faults 1 & 2 control points number given here; coordinates read in from files.
  integer(kind=4),allocatable,dimension(:) :: ftcn  !# of control points, should be given below
  real(kind=8),allocatable,dimension(:,:,:) :: fxyr   !min, max x, y coordinates for faults
!**********************************************************************************************************
  real(kind=8),allocatable,dimension(:)::x,y,h,a,b,c,d,e  !for spline calculations
  real(kind=8),allocatable,dimension(:,:)::f,xx,yy,hh
  character(len=155555)::filein
  logical,allocatable,dimension(:)::ny4flt,ny4dege,ny4abv,ny4abv2
  integer (kind=4) :: i,j,k,ix,iy,nx,ny,nxuni,nyuni,edgex1,edgex2,edgey1,edgey2,&
  	ntemp,ynplot,maxtotal,np=100000,ftn
  integer (kind=4),allocatable,dimension(:) :: line1,line2
  real (kind=8) :: pi=3.1415926535897931d0
  real (kind=8) :: tol,xstep,ystep,xcoor,ycoor,temp1,temp2,temp3,temp4,temp5,fxmin,fxmax,fymin,fymax
  real (kind=8),allocatable,dimension(:) :: xline,yline
  !
!**************************************************************************
!!!! set number of faults.
  nft = 3 !number of faults in model
!**************************************************************************
  allocate(nfn(nft),ftcn(nft),fxyr(2,2,nft),ny4flt(nft),ny4dege(nft),ny4abv(nft),ny4abv2(nft))
  ftcn(1)=15
  ftcn(2)=10
  ftcn(3)=80
  !ftcn(4)=19
  ftn=maxval(ftcn)
  allocate(f(ftn,nft),xx(ftn,nft),yy(ftn,nft),hh(ftn,nft))
  f=0.0d0
  xx=0.0d0
  yy=0.0d0
  hh=0.0d0
!...calculate spline coefficients for the faults.
  do i=1,nft
    allocate(x(ftcn(i)),y(ftcn(i)),h(ftcn(i)-1),a(ftcn(i)-2),b(ftcn(i)-2),c(ftcn(i)-2),&
             d(ftcn(i)-2),e(ftcn(i)-2))
    !...read control points (x,y) coordinate from predefined files. x-coor should be
    !   from min to max in these files: in sequence for spline! B.D. 11/23/11
    if(i==1) then
      filein='../Controlling_Points/x1_1.txt'
    elseif(i==2) then
      filein='../Controlling_Points/x2_1.txt'
    elseif (i == 3) then
	filein = '../Controlling_Points/x3_1.txt'
    elseif (i == 4) then
	filein = 'x4.txt'
    endif
    open(15,file=filein,status='old')
    read(15,*) (x(j),y(j),j=1,ftcn(i))
	close(15)
	x = x * 1.0d3
	y = y * 1.0d3
    do j=1,ftcn(i)  !save for later interpolation use.
      xx(j,i) = x(j)
      yy(j,i) = y(j)
    enddo
    fxyr(1,1,i)=minval(x)
    fxyr(2,1,i)=maxval(x)
    fxyr(1,2,i)=minval(y)
    fxyr(2,2,i)=maxval(y)
    !...form a,b,c,d for natural cubic spline.
    do j=1,ftcn(i)-1
      h(j) = x(j+1) - x(j)
      hh(j,i) = h(j)  !save for later interpolation use.
    enddo
    do j=1,ftcn(i)-2
      a(j) = h(j)
      if(j==1) a(j)=0.0d0
      b(j) = 2.0d0 * (h(j) + h(j+1))
      c(j) = h(j+1)
      if(j==ftcn(i)-2) c(j)=0.0d0
      d(j) = 6.0d0 * ((y(j+2)-y(j+1))/h(j+1) - (y(j+1)-y(j))/h(j))
    enddo
    !...forward modify coefficients.
    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)
    do j=2,ftcn(i)-2
      temp1=b(j)-c(j-1)*a(j)
      c(j) = c(j)/temp1
      d(j) = (d(j) - d(j-1)*a(j))/temp1
    enddo
    !...backward substitute for solution
    e(ftcn(i)-2) = d(ftcn(i)-2)
    do j = ftcn(i)-3,1,-1
      e(j) = d(j) - c(j) * e(j+1)
    enddo
    !...spline coefficients for ftcn(i) data points
    f(1,i)= 0
    f(ftcn(i),i) = 0
    do j=2,ftcn(i)-1
      f(j,i) = e(j-1)
    enddo
    deallocate(x,y,h,a,b,c,d,e)
  enddo ! enddo for i=1,nft
!
  !give max node #
  write(*,*) 'Max possible node #? please input:'
  read(*,*) maxtotal
  allocate(ien(5,maxtotal),xnode(3,maxtotal))
  xnode = 0.0d0
  ien = 1	!default for material properties: homogeneous material
  !...calculate derived parameters
  cellmin = dxy/2.0d0 !min grid size for dt for filing.f90
  tol = dxy/100.d0	!tolerance for comparison of equal
  fxmin=min(fxyr(1,1,1), fxyr(1,1,2), fxyr(1,1,3), fxyr(1,1,4))
  fymin=min(fxyr(1,2,1), fxyr(1,2,2), fxyr(1,2,3), fxyr(1,2,4))
  fxmax=max(fxyr(2,1,1), fxyr(2,1,2), fxyr(2,1,3), fxyr(2,1,4))
  fymax=max(fxyr(2,2,1), fxyr(2,2,2), fxyr(2,2,3), fxyr(2,2,4))
  !...determine num of nodes along x
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
  write(*,*) 'number of node along x, and y-coor',nx,ny
  pause 
  ntemp = max(nxuni,nyuni) + 1
	write(*,*) 'maxftnode = ', ntemp
  !
  allocate(line1(ny+nft),line2(ny+nft),xline(nx),yline(ny),&
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
  nfn = 0
  line1 = 0
  line2 = 0
  nsmp = 0
  dirvec = 0.0d0
  ny4flt = .false.
  ny4dege = .false.
  ny4abv = .false.
  ny4abv2 = .false.
  !...digitize along constant x line (perpendicular to main fault)
  do ix = 1, nx
    xcoor = xline(ix)
    write(*,*) 'ix,xcoor= ',ix,xcoor
    do iy = 1, ny
	ycoor = yline(iy)
        !....create nodes & elements together for general cases.
	nnode = nnode + 1
	line2(iy) = nnode
	xnode(1,nnode) = xcoor
	xnode(2,nnode) = ycoor
	if(ix >= 2 .and. iy >= 2) then
          nelement = nelement + 1
	  ien(1,nelement) = line2(iy)
	  ien(2,nelement) = line1(iy)
	  ien(3,nelement) = line1(iy-1)
	  ien(4,nelement) = line2(iy-1)	      
	  !....special treatments for elements above faults.
          do k=1,nft
  	  if(ny4flt(k)) then
  	    if(nfn(k)==1) then
  	      ien(4,nelement) = line2(ny+k)
  	    else
	      ien(3,nelement) = line1(ny+k)
	      ien(4,nelement) = line2(ny+k)
	      if(ny4abv(k)) then	!special quadrilateral above the fault
  	        ien(2,nelement) = line1(iy+1)
  	        ny4dege(k) = .true.
	        ny4abv(k) = .false.
              elseif(ny4abv2(k)) then !special quadrilateral and degeneration above fault for case 2
                ien(2,nelement) = line1(iy-1)
                nelement = nelement + 1
                ien(1,nelement) = line2(iy)
                ien(2,nelement) = line1(iy)
                ien(3,nelement) = line1(iy-1)
                ien(4,nelement) = ien(3,nelement)
                ny4abv2(k) = .false.
	      endif
	    endif
	  endif
          enddo
	endif
        do k=1,nft
	if(.not.ny4flt(k) .and. ny4dege(k)) then
	  !the current element must be degenerated to a triangular
	  ien(3,nelement) = line2(iy-1)
	  ien(4,nelement) = ien(3,nelement)
	  ny4dege(k) = .false.
	endif
        enddo
	!....conform to fault geometry, create split nodes, and fault-related elements.
        do k=1,nft  ! loop over fault #k
        ny4flt(k)=.false.
	if(xcoor>=fxyr(1,1,k).and.xcoor<=fxyr(2,1,k) .and. &
	  ycoor>=(fxyr(1,2,k)-yext).and.ycoor<=(fxyr(2,2,k)+yext)) then
	  do i=1,ftcn(k)-1
	    if(xcoor>=(xx(i,k)-tol) .and. xcoor<=(xx(i+1,k)+tol)) then
              temp1 = (f(i+1,k)*(xcoor-xx(i,k))**3+f(i,k)*(xx(i+1,k)-xcoor)**3)/(6.0d0*hh(i,k)) &
                     +(yy(i+1,k)/hh(i,k)-hh(i,k)*f(i+1,k)/6.0d0)*(xcoor-xx(i,k)) &
                     +(yy(i,k)/hh(i,k)-hh(i,k)*f(i,k)/6.0d0)*(xx(i+1,k)-xcoor)
	      if(temp1>yline(iy-1) .and. temp1<=yline(iy)) then
	        temp2 = temp1 - yline(iy-1)
	        temp3 = yline(iy) - temp1
	        if(temp3 <= temp2) then
	          ycoor = temp1
                  ny4flt(k) = .true.
	        endif
	      elseif(temp1>yline(iy) .and. temp1<=yline(iy+1)) then
	        temp2 = temp1 - yline(iy)
	        temp3 = yline(iy+1) - temp1
	        if(temp2 <= temp3) then
	          ycoor = temp1
                  ny4flt(k) = .true.
	        endif
	      endif
	      if(ny4flt(k)) then
	        xnode(2,nnode) = ycoor
	        nfn(k) = nfn(k) + 1
	        nsmp(1,nfn(k),k) = nnode	!slave node
	        nnode = nnode + 1
	        nsmp(2,nfn(k),k) = nnode	!master node
	        line2(ny+k) = nnode
	        xnode(1,nnode) = xcoor
	        xnode(2,nnode) = ycoor
	        if(nfn(k)>1) then	!1st fault node, no need for this
                  ny4abv(k) = .false.
                  ny4abv2(k) = .false.
	          if((xnode(2,nsmp(1,nfn(k)-1,k))-xnode(2,line1(iy)))>tol) then !for down slope
	          	!use ycoor to judge fault trace
	            !first, the current element is still quadrilateral, but with special nodes.
	 	    ien(1,nelement) = line2(iy)	!this is same as usual
	  	    ien(2,nelement) = line1(iy+1)	!the following two are special
	  	    ien(3,nelement) = line1(iy)
	  	    ien(4,nelement) = line2(iy-1)	!this is same as usual, too
	  	    !then,one more element is created by degeration, which is a triangular element.
	            nelement = nelement + 1
	            ien(1,nelement) = line1(iy)
	            ien(2,nelement) = line1(iy-1)
	            ien(3,nelement) = line2(iy-1)
	            ien(4,nelement) = ien(3,nelement)
                    ny4abv(k) = .true.	!for special treatment above the fault
                  elseif((xnode(2,line1(iy))-xnode(2,nsmp(1,nfn(k)-1,k)))>tol) then ! for up slope
                    !first, previous element degenerated to triangular
                    ien(1,nelement-1) = line2(iy-1)
                    ien(2,nelement-1) = line1(iy-2)
                    ien(3,nelement-1) = line2(iy-2)
                    ien(4,nelement-1) = ien(3,nelement-1)
                    !then, current special quadrilateral
                    ien(1,nelement) = line2(iy)
                    ien(2,nelement) = line1(iy-1)
                    ien(3,nelement) = line1(iy-2)
                    ien(4,nelement) = line2(iy-1)
                    ny4abv2(k) = .true. !for special treatment above fault: case 2 up slope
                  endif	            
	        endif
	        !...cal unit normal and tangent at each split-nodes pair: using spline
                  temp2=(f(i,k)*(xx(i+1,k)-xcoor)**2-f(i+1,k)*(xcoor-xx(i,k))**2)/(2*hh(i,k)) &
                       +(6.0d0*(yy(i,k)-yy(i+1,k))+hh(i,k)*hh(i,k)*(f(i+1,k)-f(i,k)))/(6.0d0*hh(i,k))
                  temp3=sqrt(temp2*temp2+1)
                  dirvec(1,nfn(k),k)=1/temp3       !tx: x-comp unit tangent
                  dirvec(2,nfn(k),k)=-temp2/temp3  !ty: y-comp unit tangent
                if(nfn(k)>2) then ! legnth associated with the previous fault node: linear segment
	          temp2 = xnode(1,nsmp(1,nfn(k)-1,k)) - xnode(1,nsmp(1,nfn(k)-2,k))
	          temp3 = xnode(2,nsmp(1,nfn(k)-1,k)) - xnode(2,nsmp(1,nfn(k)-2,k))
	          temp4 = sqrt(temp2*temp2+temp3*temp3)
	          temp2 = xnode(1,nsmp(1,nfn(k),k)) - xnode(1,nsmp(1,nfn(k)-1,k))
	          temp3 = xnode(2,nsmp(1,nfn(k),k)) - xnode(2,nsmp(1,nfn(k)-1,k))
	          temp5 = sqrt(temp2*temp2+temp3*temp3)
	          dirvec(3,nfn(k)-1,k) = 0.5d0 * (temp4 + temp5)	!length associated with node
	        endif
	        exit
	      endif	       
	    endif
	  enddo
	endif
        enddo  !enddo for k=1,nft
	!
    enddo	!iy
    line1 = line2
  enddo		!ix
!
  !...special calculation for fault split node's length at two ends of each fault
  do i=1,nft
    do j=1,2
      if(j==1) then
        temp2 = xnode(1,nsmp(1,2,i)) - xnode(1,nsmp(1,1,i))
        temp3 = xnode(2,nsmp(1,2,i)) - xnode(2,nsmp(1,1,i))
        temp4 = sqrt(temp2*temp2+temp3*temp3)
        k = 1
      else
        temp2 = xnode(1,nsmp(1,nfn(i),i)) - xnode(1,nsmp(1,nfn(i)-1,i))
        temp3 = xnode(2,nsmp(1,nfn(i),i)) - xnode(2,nsmp(1,nfn(i)-1,i))
        temp4 = sqrt(temp2*temp2+temp3*temp3)
        k = nfn(i)
      endif      
      dirvec(3,k,i) = 0.5d0 * temp4	!txy: 'length' of the node
    enddo
  enddo

  write(*,*) 'Entire model region has:'
  write(*,*) nnode,' nodes; ',nelement,' elements.'
  
  write(*,*) 'write out data for plotting?'
  write(*,*) '1 - yes,otherwise - no'
  read(*,*) ynplot
 ! ynplot = 0	!no output by default
  !***
  !...write out mesh data for matlab to plot
  !***
  if(ynplot==1) then
    open(11,file='vert.txt',status='unknown')
    write(11,'(1x,2e18.7e4)') ((xnode(j,i),j=1,2),i=1,nnode)
    close(11)
    open(12,file='fac.txt',status='unknown')
    write(12,'(1x,4i15)') ((ien(j,i),j=1,4),i=1,nelement)
    close(12)
  endif

!  write(*,*) 'plotting done, stop!'
!  pause
  
  open(13,file='Mesh_general_info.txt',status='unknown')
  write(13,'( 2i15)') nnode,nelement
  write(13,'( 4i15)') (nfn(k),k=1,nft)
 ! write(13,'( a60)') 'fault nodes'' x, y- xoor'
 ! write(13,'( a20)') 'fault #1'
 ! write(13,'( i10,2f12.1)') (nsmp(1,i,1),(xnode(j,nsmp(1,i,1)),j=1,2),i=1,nfn(1))
 ! write(13,'( a20)') 'fault #2'
 ! write(13,'( i10,2f12.1)') (nsmp(1,i,2),(xnode(j,nsmp(1,i,2)),j=1,2),i=1,nfn(2))
  close(13)

  ntemp = maxval(nfn)
  open(14, file = 'nsmp.txt', form = 'formatted', status ='unknown')
  do j = 1, nft 
         write(14,'(2i15)') (nsmp(1,i,j), nsmp(2,i,j), i = 1, ntemp)
  enddo
  close(14)

  open(14, file = 'nsmpnv.txt', form = 'formatted', status = 'unknown')
  do j = 1, nft
        write(14, '(3e18.7e4)') (-dirvec(2,i,j), dirvec(1,i,j), dirvec(3,i,j), i = 1, ntemp)
  enddo
  close(14)

!  
end subroutine meshgen

	          
	  
	  
