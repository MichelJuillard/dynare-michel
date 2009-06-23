subroutine TV(X,T,V,n)
!     COMPUTATIONAL SUBROUTINE: TV=T*V T upper triangular; V vector

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: V(n,1)

real(8), intent(out) :: X(n,1) 

integer(4) :: i,k

real(8) :: stemp



do i=1,n
stemp = 0.0
do  k = 0,n-i
stemp = stemp + T(i,k+i) * V(k+i,1)
end do
X(i,1)=stemp
end do

return
end subroutine TV
