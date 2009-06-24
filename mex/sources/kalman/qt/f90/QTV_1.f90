subroutine QTV_1(X,QT,V,n)
!     COMPUTATIONAL SUBROUTINE: X=QT*V QT upper quasi-triangular; V vector

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: QT(n,n)
real(8), intent(in) :: V(n,1)

real(8), intent(out) :: X(n,1) 

integer(4) :: i,k

real(8) :: stemp

X=0*V

do i=1,n
stemp = 0.0
if (i > 1 .AND. (QT(i,i-1)/= 0.0)) then
stemp = QT(i,i-1)*v(i-1,1)
end if 

do  k = 0,n-i
stemp = stemp + QT(i,k+i) * V(k+i,1)
end do
X(i,1)=stemp
end do

return
end subroutine QTV_1


