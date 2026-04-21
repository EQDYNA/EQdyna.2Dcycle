SUBROUTINE slip_weak(slip,frictmp,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! friction law for fault dynamics. B.D. 8/19/06
  !
  real (kind = dp) :: xmu,slip
  real (kind = dp),dimension(2) :: frictmp
  !
  if(slip == 0.0) then
    xmu = frictmp(1)	!xmu is frictmptional coefficient, node by node on fault
  elseif(slip < critd0) then
    xmu = frictmp(1) - (frictmp(1) - frictmp(2))*slip/critd0
  else
    xmu = frictmp(2)
  endif
  !
end SUBROUTINE slip_weak
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE slip_rate_weak(slip,rate,frictmp,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! and rate-weakening frictmption law for fault dynamics.
  ! original references of the law are: Madariaga et al.
  ! (1998) and Aagaard et al. (2001): same linear form for
  ! slip-weakening law and rate-weakening, but the largest
  ! one to be used at a given time.
  ! B.D. 12/29/16
  !
  real (kind = dp) :: xmu,slip,rate,xmu1,xmu2
  real (kind = dp),dimension(3) :: frictmp
  !
  !...slip-weakening
  if(slip == 0.0) then
    xmu1 = frictmp(1) !xmu is frictmptional coefficient, node by node on fault
  elseif(slip < critd0) then
    xmu1 = frictmp(1) - (frictmp(1) - frictmp(2))*slip/critd0
  else
    xmu1 = frictmp(2)
  endif
  !...rate-weakening
  if(rate == 0.0) then
    xmu2 = frictmp(3) !restrengthing mu (mur), should be smaller than mus (frictmp(1)). 
  elseif(rate < critv0) then
    xmu2 = frictmp(3) - (frictmp(3) - frictmp(2))*rate/critv0
  else
    xmu2 = frictmp(2)
  endif
  !...actual mu is the bigger of the above two
  xmu = max(xmu1,xmu2)
  !
end SUBROUTINE slip_rate_weak

!

!================================================
SUBROUTINE slip_weak_healing(slip,frictmp,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! frictmption law for fault dynamics. B.D. 8/19/06
  !
  real (kind = dp) :: xmu,slip
  real (kind = dp),dimension(2) :: frictmp
  !
  if(slip == 0.0) then
    xmu = frictmp(1)	!xmu is frictmptional coefficient, node by node on fault
  elseif(slip < critd0) then
    xmu = frictmp(1) - (frictmp(1) - frictmp(2))*slip/critd0
  else
    xmu =0.56 
  endif
  !
end SUBROUTINE slip_weak_healing
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE slip_rate_caltech(slip,sliprate,frictmp,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! frictmption law for fault dynamics. B.D. 8/19/06
  !
  real (kind = dp) :: xmu,slip,sliprate

  real (kind = dp),dimension(2) :: frictmp
  !
  if(slip == 0.0) then
    xmu = frictmp(1)	!xmu is frictmptional coefficient, node by node on fault
  elseif(slip < critd0) then
    xmu = frictmp(1) - (frictmp(1) - frictmp(2))*slip/critd0
!	write(*,*)"I am here !"
  elseif(sliprate<critv0 .and. slip >=critd0)then
    xmu =0.6 
  else
	xmu=frictmp(2)  
  endif
  !
end SUBROUTINE slip_rate_caltech
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE slip_rate_zaifeng(slip,trupt,frictmp,xmu)
!  use globalvar
!  implicit none
!  !
!  !### subroutine to implement linear slip-weakening
!  ! frictmption law for fault dynamics. B.D. 8/19/06
!  !
!  real (kind = dp) :: xmu,slip,trupt
!  real (kind = dp):: critv0=0.0
!  real (kind = dp),dimension(2) :: frictmp
!  !
!  if(slip == 0.0) then
!    xmu = frictmp(1)	!xmu is frictmptional coefficient, node by node on fault
!  elseif(slip < critd0) then
!    xmu = frictmp(1) - (frictmp(1) - frictmp(2))*slip/critd0
!  elseif(trupt<(2*critt0) .and. slip>=critd0 )then
!   	xmu=frictmp(2)+(frictmp(1)-frictmp(2))*(trupt-critt0)/critt0
!  else
!	xmu=frictmp(1)  
!  endif
!  !
!end SUBROUTINE slip_rate_zaifeng
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE time_weak(trupt,frictmp,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear time-weakening
  ! frictmption law for fault dynamics. B.D. 8/19/06
  !
  real (kind = dp) :: xmu,trupt
  real (kind = dp),dimension(2) :: frictmp
  !
  if(trupt <= 0.0) then
    xmu = frictmp(1)
  elseif(trupt < critt0) then
    xmu = frictmp(1) - (frictmp(1) - frictmp(2))*trupt/critt0
  else
    xmu = frictmp(2)
  endif
  !
end SUBROUTINE time_weak

