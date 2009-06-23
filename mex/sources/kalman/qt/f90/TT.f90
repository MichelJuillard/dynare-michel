subroutine TT(X, T1,T2,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*Ut T upper triangular; U striclty lower triangular

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T1(n,n)
real(8), intent(in) :: T2(n,n)

real(8), intent(out) :: X(n,n)
integer(4) :: i,j,k,h

real(8) :: stemp

X=0.0*T1

do i=1,n
do j=1,n
stemp = 0.0
do  k = 0,n-i
stemp = stemp + T1(i,k+i) * T2(k+i,j)
end do
X(i,j)=stemp
end do

end do

return
end subroutine TT