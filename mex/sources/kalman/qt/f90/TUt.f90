subroutine TUt(X,T,Ut,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*Ut T upper triangular; U striclty lower triangular

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: Ut(n,n)

real(8), intent(out)  :: X(n,n)
integer(4) :: i,j,k,h

real(8) :: stemp


X=0.0*T
do i=1,n
do j=1,n-1
h=max(i,j)
stemp = 0.0
do  k = 0,n-h
stemp = stemp + T(i,k+h) * Ut(k+h,j)
end do
X(i,j)=stemp
end do
end do

return
end subroutine TUt