subroutine catmullrom(p, xxx, CCC, alpha0)
	
	use globalvar
	implicit none
	
	real (kind = 8) :: p(2,4), xxx, y, alpha0, t(4), tt, A1, A2, A3, B1, B2, CCC
	integer (kind = 4) :: i
	
	t = 0.0d0
	do i = 2,4
		t(i) = ((p(1,i)-p(1,i-1))**2 + (p(2,i)-p(2,i-1))**2)**(alpha0/2.0d0) + t(i-1) 
	enddo 
	
	tt = (t(3)-t(2))* (xxx-p(1,2))/(p(1,3)-p(1,2)) + t(2)
	
	A1 = (t(2) - tt)/(t(2)-t(1))*p(2,1) + (tt - t(1))/(t(2) - t(1))*p(2,2)
	A2 = (t(3) - tt)/(t(3)-t(2))*p(2,2) + (tt - t(2))/(t(3) - t(2))*p(2,3)
	A3 = (t(4) - tt)/(t(4)-t(3))*p(2,3) + (tt - t(3))/(t(4) - t(3))*p(2,4)
	
	B1 = (t(3) - tt)/(t(3)-t(1))*A1 + (tt - t(1))/(t(3) - t(1))*A2
	B2 = (t(4) - tt)/(t(4)-t(2))*A2 + (tt - t(2))/(t(4) - t(2))*A3
	
	CCC = (t(3) - tt)/(t(3)-t(2))*B1 + (tt - t(2))/(t(3) - t(2))*B2
	
end subroutine catmullrom