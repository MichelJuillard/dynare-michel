subroutine LdM(X,Ld,M,n)
!     COMPUTATIONAL SUBROUTINE: LdM=Ld*M; Ld lower diagonal matrix; M arbitrary

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: Ld(n,n)
real(8), intent(in) :: M(n,n)
real(8), intent(out) :: X(n,n)
integer(4) :: i,j


X=0.0*M

do i=2,n
if ((Ld(i,i-1)/=0.0)) then
do j=1,n
X(i,j)=Ld(i,i-1)*M(i-1,j)
end do
end if
end do

return
end subroutine LdM