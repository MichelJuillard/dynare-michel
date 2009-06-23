subroutine QT2Ld(X,QT,n)
!     COMPUTATIONAL SUBROUTINE: extracts lower diagonal from Quasi triangular matrix

!     COMPUTATIONAL SUBROUTINE: extracts lower diagonal from Quasi triangular matrix
implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: QT(n,n)
real(8), intent(out) :: X(n,n)
integer(4) :: i





X=0.0*QT


do i=2,n
if (QT(i,i-1)/=0.0) then
X(i,i-1)=QT(i,i-1)
end if
end do

return
end subroutine QT2Ld
