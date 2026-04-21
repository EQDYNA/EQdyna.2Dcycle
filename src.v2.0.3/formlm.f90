SUBROUTINE formlm
	use globalvar
	implicit none
	!
	!### program to form lm array
	!
	integer (kind=4) :: i,j,k,node
	!
	do k=1,numel
		do j=1,nen
			node=ien(j,k)
			do i=1,ned
				lm(i,j,k) = id(i,node)
			enddo
		enddo
	enddo
end SUBROUTINE formlm
