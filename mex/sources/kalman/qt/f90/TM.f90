subroutine TM(X,T,M,n)
!     COMPUTATIONAL SUBROUTINE: TM=T*M T upper triangular; M arbitrary


implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: M(n,n)

real(8), intent(out) :: X(n,n)

integer(4) :: i,j,k

real(8) :: stemp


do i=1,n
do j=1,n
stemp = 0.0
do  k = 0,n-i
stemp = stemp + T(i,k+i) * M(k+i,j)
end do
X(i,j)=stemp
end do
end do

return
end subroutine TM
