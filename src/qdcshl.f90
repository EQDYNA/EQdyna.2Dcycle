SUBROUTINE qdcshl
  use globalvar
  implicit none
  !
  !### program to calculate integration-rule weights, shape functions
  !        and local derivatives for a four-node quadrilateral element
  !
  !               s,t = local element coordinates ("xi", "eta", resp.)
  !        shl(1,i,l) = local ("xi") derivative of shape function
  !        shl(2,i,l) = local ("eta") derivative of shape function
  !        shl(3,i,l) = local  shape function
  !              w(l) = integration-rule weight
  !                 i = local node number
  !                 l = integration point number
  !              nint = number of integration points, eq. 1 or 4
  !
  integer (kind=4) :: l,i
  real (kind = dp),dimension(4) :: ra=(/-0.5d0,0.5d0,0.5d0,-0.5d0/), &
  				sa=(/-0.5d0,-0.5d0,0.5d0,0.5d0/)
  real (kind = dp) :: g,r,s,tempr,temps
  !
  g = 0.0d0
  w(1) = 4.0d0
  if (nint == 4) then
    g = 2.0d0/sqrt(3.0d0)
    w(1) = 1.0d0
    w(2) = 1.0d0
    w(3) = 1.0d0
    w(4) = 1.0d0
  endif
  !
  do l=1,nint
    r = g*ra(l)
    s = g*sa(l)
    do i=1,4
      tempr = 0.5d0 + ra(i)*r
      temps = 0.5d0 + sa(i)*s
      shl(1,i,l) = ra(i)*temps
      shl(2,i,l) = tempr*sa(i)
      shl(3,i,l) = tempr*temps
    enddo
  enddo
  !
end SUBROUTINE qdcshl
