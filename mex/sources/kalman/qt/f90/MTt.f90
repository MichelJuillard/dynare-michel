subroutine MTt(X,M,Tt,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*Ut T upper triangular; U striclty lower triangular

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: M(n,n)
real(8), intent(in) :: Tt(n,n)
real(8), intent(out) :: X(n,n)

integer(4) :: i,j,k

real(8) :: stemp



X=0.0*M

do i=1,n
do j=1,n
stemp = 0.0
do  k = j,n
stemp = stemp + M(i,k+i) * Tt(k+i,j)
end do
X(i,j)=stemp
end do

end do

return
end subroutine MTt