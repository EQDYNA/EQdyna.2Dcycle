SUBROUTINE qdcshg(xl,det,shg,nel,lquad,xs)
	use globalvar
	implicit none
	!====================================================================
	!### program to calculate global derivatives of shape functions and
	!        jacobian determinants for a four-node quadrilateral element
	!
	!        xl(j,i)    = global coordinates
	!        det(l)     = jacobian determinant
	!        shl(1,i,l) = local ("xi") derivative of shape function
	!        shl(2,i,l) = local ("eta") derivative of shape function
	!        shl(3,i,l) = local  shape function
	!        shg(1,i,l) = x-derivative of shape function
	!        shg(2,i,l) = y-derivative of shape function
	!        shg(3,i,l) = shl(3,i,l)
	!        xs(i,j)    = jacobian matrix
	!                 i = local node number or global coordinate number
	!                 j = global coordinate number
	!                 l = integration-point number
	!              nint = number of integration points, eq. 1 or 4
	!=====================================================================
	logical :: lquad
	integer (kind=4) :: i,j,k,l,nel
	real (kind = dp) :: temp
	real (kind = dp),dimension(nint) :: det
	real (kind = dp),dimension(2,2) :: xs
	real (kind = dp),dimension(nesd,nen) :: xl
	real (kind = dp),dimension(nrowsh,nen,nint) :: shg  
	!
	!
	!*** equal local to global,first. ***
	!
	do l=1,nint
		do i=1,nen
			do k=1,3
				shg(k,i,l) = shl(k,i,l)
			enddo
		enddo
	enddo
	!
	!*** loop over integration points ***
	!
	do l=1,nint
		!......degenration to triangular
		if (.not.lquad) then	!degeneration
			do i=1,3
				shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
				shg(i,4,l) = 0.0d0
			enddo
		endif
		!......pre-Jacobi matrix
		do j=1,2
			do i=1,2
				temp = 0.0d0
				do k=1,nen
					temp = temp + shg(i,k,l) * xl(j,k)
				enddo
				xs(i,j) = temp
			enddo
		enddo
		!......determinant
		det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
		if (det(l) <= 0.0d0) then
			write(*,1000) nel
			write(*,*) xl(1,1), x(2,1)
			write(*,*) xl(1,2), x(2,2)
			write(*,*) xl(1,3), x(2,3)
			write(*,*) xl(1,4), x(2,4)
			stop	!non-positive det results in termination
		endif
		!......Jacobi matrix
		do j=1,2
			do i=1,2
				xs(i,j) = xs(i,j)/det(l)
			enddo
		enddo
		!......global derivatives of shape function
		do i=1,nen
			temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
			shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
			shg(1,i,l) = temp
		enddo
	enddo
	!
	1000 format('1','non-positive determinant in element number  ',i10)
	!
end SUBROUTINE qdcshg
