subroutine TD(X,T,D,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*D T upper triangular; D diagonal

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: D(n,n)

real(8), intent(out) :: X(n,n)
integer(4) :: i,j


X=0.0*T
do i=1,n
do j=i,n
X(i,j)=T(i,j)*D(j,j)
end do
end do

return
end subroutine TD
