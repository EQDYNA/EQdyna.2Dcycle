MODULE globalvar
  implicit none
  save
  integer, parameter :: dp = selected_real_kind(15,307)
  !...inout/output units
  integer (kind=4) :: iingl=11,iinft=12,ioutgl=13, &
  	   ioutrt=15,ioutst=16,ioutsl=17,iinpa=18,iinst=19
  !...measure portion CPU time
  real (kind = dp), dimension (6) :: timeused=0.0d0
  integer (kind=4) :: gethrtime, time1, time2
  real (kind = dp) :: timedyna=0.0d0, cellmin
  external gethrtime	!this function in the SUN Fortran library
  !...execution control parameters	   	   
  integer (kind=4) :: iexec,irank,numseq,ndout,nsd = 2,numnp,ndof = 2, &
           nlvect,nltftn,nptslf,numeg
  !...scalar variables (eps. for quadrilateral elements)
  integer (kind=4) :: neq,nen=4,ned=2,nee=8,nesd=2,nrowsh=3, &
           nrowb=3,nrowc=4,nstr=3,nint,numel,imass,numat = 1	   
  !...element arrays
  integer (kind=4),allocatable,dimension(:) :: mat	    
  integer (kind=4),allocatable,dimension(:,:) :: ien,id
  integer (kind=4),allocatable,dimension(:,:,:) :: lm
  real (kind = dp),allocatable,dimension(:) :: eledet
  real (kind = dp),allocatable,dimension(:,:) :: elemass
  real (kind = dp),allocatable,dimension(:,:,:) :: eleb
  !... equations,node mass,weight,solutions,coordinate,shape
  real (kind = dp),allocatable,dimension(:) :: alhs,brhs,fnms,w
  real (kind = dp),allocatable,dimension(:,:) :: d,v
  real (kind = dp),allocatable,dimension(:,:) :: x
  real (kind = dp),allocatable,dimension(:,:,:) :: shl
  !...definitions for faults
  integer (kind=4) ::  nfnodemx, iopt = 1
  real (kind = dp),allocatable,dimension(:) :: ftimestr
  integer (kind=4) :: nplpts = 1, nuciloc(20), numnuc, nstep1, sumnfnode
  integer (kind=4),allocatable,dimension(:) :: nstep,ndprt,nsprt,nhplt,niter
  real (kind=4),allocatable,dimension(:) :: alpha,beta,gamma,dt
   !...material properties
  real (kind = dp),allocatable,dimension(:) :: e,rho,rdampm,rdampk,th, pois
  real (kind = dp),allocatable,dimension(:,:,:) :: c
  !...time histories	   	    
  integer (kind=4) :: locplt
  integer (kind=4),allocatable,dimension(:,:) :: idhist
  real (kind = dp),allocatable,dimension(:,:) :: dout
  !...hourglass control arrays
  real (kind = dp),allocatable,dimension(:,:) :: ss,phi
  logical :: quitdriver = .FALSE.
  integer (kind = 4):: maxtotal
  real (kind = dp) :: pi=3.1415926535897931d0
  real (kind = dp) :: dxy, yext, xmin, xmax, ymin, ymax, rat
  real (kind = dp) :: fric_fs, fric_fd, fric_fv, fric_fini
  real (kind = dp) :: critd0, critv0, critt0, vrupt0, radius, term, dt1
  real (kind = dp) :: amu, ant0, str, rou, youngs, poisr
  integer (kind = 4) :: C_mesh, ntotft, friclaw, ien0(5,10000000),nsmp0(2,10000)
  integer (kind = 8) :: totftnode, maxftnode, loc, icstart, icend, debug, ftn, plotmesh
  real (kind = dp) :: xnode0(2,10000000), nsmpnv(3, 10000), fnft(10000), fric(50,10000)
  real (kind = dp) :: ambientnorm, minnorm = -10.0d6, vp, vs, lambda
  real (kind = dp),allocatable, dimension(:,:) :: fistr, output4plot, rd! Shear, normal, slip, sliprate, miu
  integer (kind=8),allocatable, dimension(:) :: nfnode, ftcn
  
	     
end MODULE globalvar 	      
