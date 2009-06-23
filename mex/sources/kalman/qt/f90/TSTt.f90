subroutine TSTt(X,T,S,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*Ut T upper triangular; U striclty lower triangular

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: T(n,n)
real(8), intent(in) :: S(n,n)

real(8), intent(out) :: X(n,n)
integer(4) :: i,j,k,h

real(8) :: stemp


X=0.0*S

do i=1,n
do j=1,n
stemp = 0.0
do  h = i,n
do  k = j,n
stemp = stemp + T(i,h) * S(h,k) * T(j,k)
end do
end do
X(i,j)=stemp
X(j,i)=stemp

end do

end do

return
end subroutine TSTt
