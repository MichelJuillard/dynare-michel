subroutine QTSQTt(X,QT,S,n)
!     COMPUTATIONAL SUBROUTINE X=QT*SQT' QT upper quasi-triangular; S symmetric

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: QT(n,n)
real(8), intent(in) :: S(n,n)

real(8), intent(out) :: X(n,n)
integer(4) :: i,j,k,h,k0

real(8) :: stemp

X=0*S
do  i=1,n
do j=i, n
if (i > 1 .AND. (QT(i,i-1)/= 0.0)) then
h=i-1
else
h=i
end if 
if (j > 1 .AND. (QT(j,j-1)/= 0.0)) then
k=j-1
k0=k 
else
k=j
k0=k
end if 
stemp=0
do while (h <= n)
do while (k <= n)
stemp=stemp+QT(i,h)*S(h,k)*QT(j,k)
k=k+1
end do
k=k0
h=h+1
end do
X(i,j)=stemp
if (i /= j) then
X(j,i)=stemp
end if 
end do
end do




return
end subroutine QTSQTt
