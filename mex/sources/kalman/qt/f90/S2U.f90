subroutine S2U(X,S,n)
!     COMPUTATIONAL SUBROUTINE: S2U extracts striclty upper triangular from symmetric

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: S(n,n)
real(8), intent(out) :: X(n,n)
integer(4) :: i,j





X=0.0*S
do i=1,n
do j=i+1,n
X(i,j)=S(i,j)
end do
end do

return
end subroutine S2U
