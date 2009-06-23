subroutine QT2T(X,QT,n)
!     COMPUTATIONAL SUBROUTINE: extracts upper triangular from Quasi triangular matrix

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: QT(n,n)
real(8), intent(out) :: X(n,n)
integer(4) :: i,j


X=0.0*QT
do i=1,n
do j=i,n
X(i,j)=QT(i,j)
end do
end do

return
end subroutine QT2T
