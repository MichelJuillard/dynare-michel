subroutine TU(X,T,U,n)
!     COMPUTATIONAL SUBROUTINE: TU=T*U; T upper triangular matrix; U strictly upper triangular


implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: U(n,n)

real(8), intent(out) :: X(n,n)
integer(4) :: i,j,k

real(8) :: stemp


X=0.0*T


do i=1,n-1
do j=i+1,n
stemp = 0.0
do  k = i,j-1
stemp = stemp + T(i,k) * U(k,j)
end do
X(i,j)=stemp
end do
end do

return
end subroutine TU
