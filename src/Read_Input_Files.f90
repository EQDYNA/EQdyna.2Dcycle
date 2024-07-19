! #1 readgloabl -----------------------------------------------------------
subroutine readglobal
! This subroutine is read information from FE_global.txt
	use globalvar
	implicit none
	integer (kind = 4) :: i
	logical::file_exists
	
	INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
	write(*,*) '=                                                                   ='
	write(*,*) '=     Checking in FE_global.txt ...                                 ='
	if (file_exists .eqv. .false.) then
		write(*,*) '=                                                                   ='
		write(*,*) '=     FE_Global.txt is required but missing ...                     ='	
		 
	endif 
	!INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
	! if (file_exists .eqv. .false.) then
		! write(*,*) 'FE_Global.txt is still missing, so exiting EQdyna'
		!  
	! endif 
	
	open(unit = 1001, file = 'FE_Global.txt', form = 'formatted', status = 'old')
		read(1001,*) C_mesh
		read(1001,*) ntotft
		read(1001,*) friclaw
		read(1001,*)
		read(1001,*) dt1
		read(1001,*) term 		
		read(1001,*) icstart, icend
		read(1001,*)
		read(1001,*) fric_fs
		read(1001,*) fric_fd
		read(1001,*) fric_fv
		read(1001,*) fric_fini
		read(1001,*) critd0
		read(1001,*) critv0
		read(1001,*) critt0
		read(1001,*) vrupt0
		read(1001,*) radius 
		read(1001,*)
		read(1001,*) vp
		read(1001,*) vs
		read(1001,*) rou
		read(1001,*) eta0 !2.8d21*2.0d0 ! Viscosity 
		read(1001,*) maxShearStrainLoadRate  !1.427d-14 ! Maximum shearing strain loading rate	
		read(1001,*) 
		read(1001,*) ambientnorm
		read(1001,*) debug
		read(1001,*) plotmesh
	close(1001)
end subroutine readglobal 
! #2 readmodelgeometry -------------------------------------------------
subroutine readmodelgeometry
! This subroutine is read information from FE_global.txt
	use globalvar
	implicit none

	logical::file_exists
	
	INQUIRE(FILE="FE_Model_Geometry.txt", EXIST=file_exists)
	write(*,*) '=                                                                   ='
	write(*,*) '=     Checking in FE_Model_Geometry.txt ...                         ='	
	if (file_exists .eqv. .false.) then
		write(*,*) '=                                                                   ='
		write(*,*) '=     FE_Model_Geometry.txt is required but missing ...             ='		
		 
	endif 
	!INQUIRE(FILE="FE_Model_Geometry.txt", EXIST=file_exists)
	! if (file_exists .eqv. .false.) then
		! write(*,*) 'FE_Model_Geometry.txt is still missing, so exiting EQdyna'
		!  
	! endif 
	
	open(unit = 1002, file = 'FE_Model_Geometry.txt', form = 'formatted', status = 'old')
		! read(1002,*) xmin, xmax
		! read(1002,*) ymin, ymax
		! read(1002,*) 
		read(1002,*) yext
		read(1002,*) rat
		read(1002,*) dxy 
	close(1002)
end subroutine readmodelgeometry

! #3 readfaultgeometry
subroutine readfaultgeometry
! This subroutine is read information from FE_Fault_Geometry.txt
	use globalvar
	implicit none

	logical::file_exists
	integer(kind=4)::me,i
	
	INQUIRE(FILE="FE_Fault_Geometry.txt", EXIST=file_exists)
	write(*,*) '=                                                                   ='
	write(*,*) '=     Checking in FE_Fault_Geometry.txt ...                         ='
	if (file_exists .eqv. .false.) then
		write(*,*) '=                                                                   ='
		write(*,*) '=     FE_Fault_Geometry.txt is required but missing ...             ='	
		 
	endif 
	! INQUIRE(FILE="FE_Fault_Geometry.txt", EXIST=file_exists)
	! if (file_exists .eqv. .false.) then
		! write(*,*) 'FE_Fault_Geometry.txt is still missing, so exiting EQdyna'
		!  
	! endif 
	
	open(unit = 1003, file = 'FE_Fault_Geometry.txt', form = 'formatted', status = 'old')
		read(1003,*) (ftcn(i), i = 1,ntotft)
	close(1003)
end subroutine readfaultgeometry

! #4 readmaterial --------------------------------------------------------
subroutine readmaterial
! This subroutine is read information from FE_Material.txt
	use globalvar
	implicit none

	logical::file_exists
	integer(kind=4):: me, i, j 
	
	
	INQUIRE(FILE="FE_Material.txt", EXIST=file_exists)
	write(*,*) 'Checking FE_Material.txt by the master procs'
	if (file_exists .eqv. .false.) then
		write(*,*) 'FE_Material.txt is required but missing ...'
		 
	endif 

	INQUIRE(FILE="FE_Material.txt", EXIST=file_exists)
	if (file_exists .eqv. .false.) then
		write(*,*) 'FE_Material.txt is still missing, so exiting EQdyna'
		 
	endif 

	
	open(unit = 1004, file = 'FE_Material.txt', form = 'formatted', status = 'old')
			read(1004,*) 
	close(1002)
end subroutine readmaterial

! #5 readfric --------------------------------------------------------
subroutine readfric
! This subroutine is read information from FE_Material.txt
	use globalvar
	implicit none


	logical::file_exists
	integer(kind=4):: me, i, j 
	
	
	INQUIRE(FILE="FE_Fric.txt", EXIST=file_exists)
	write(*,*) 'Checking FE_Fric.txt by the master procs'
	if (file_exists .eqv. .false.) then
		write(*,*) 'FE_Fric.txt is required but missing ...'
		 
	endif 

	INQUIRE(FILE="FE_Fric.txt", EXIST=file_exists)
	if (file_exists .eqv. .false.) then
		write(*,*) 'FE_Fric.txt is still missing, so exiting EQdyna'
		 
	endif 

	
	open(unit = 1005, file = 'FE_Fric.txt', form = 'formatted', status = 'old')
	if (friclaw == 1) then ! slip-weakening
		
	endif 
	if (friclaw == 3) then ! aging law
		do i = 1, 4
			read(1005,*)
		enddo
		
	endif 	
	if (friclaw == 4) then ! strong-rate weakening
		do i = 1, 4
			read(1005,*)
		enddo
		
	endif 		
	if (friclaw == 5) then ! strong-rate weakening + termop
		do i = 1, 4
			read(1005,*)
		enddo
	
	endif 	
	close(1005)
end subroutine readfric