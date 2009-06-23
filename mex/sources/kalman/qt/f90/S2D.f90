subroutine S2D(X,S,n)
!     COMPUTATIONAL SUBROUTINE TUt=T*Ut T upper triangular; U striclty upper triangular

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: S(n,n)
real(8), intent(out) :: X(n,n)

integer(4) :: i



X=0.0*S
do i=1,n
X(i,i)=sqrt(S(i,i))
end do
return
end subroutine S2D
