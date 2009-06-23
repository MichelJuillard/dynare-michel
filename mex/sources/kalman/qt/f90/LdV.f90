subroutine LdV(X,Ld,v,n)
implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: ld(n,n)
real(8), intent(in) :: v(n,1)
real(8), intent(out) :: X(n,1)
integer(4) :: i

X=0*v

do i=2,n
X(i,1)=Ld(i,i-1)*v(i-1,1)
end do

return
end subroutine LdV
