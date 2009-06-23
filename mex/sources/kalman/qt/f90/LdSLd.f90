subroutine LdSLd(X,Ld,S,n)
!     COMPUTATIONAL SUBROUTINE: LdM=Ld*S*Ld; Ld lower diagonal matrix; S symmetric

implicit none

integer(4), intent(in) :: n
real(8), intent(in) :: Ld(n,n)
real(8), intent(in) :: S(n,n)
real(8), intent(out) :: X(n,n)
integer(4), dimension(1,n) :: vk


integer(4) :: i,j, jj, h


h=0.0
do  i=2,n
if ((Ld(i,i-1)/= 0.0)) then
h=h+1.0
vk(1,h)=i
end if
end do
if (h==0.0) then
X=0*S
return
end if
do  j=1,h
do  jj=1,h
X(vk(1,j),vk(1,jj))=Ld(vk(1,j),vk(1,j)-1)*S(vk(1,j)-1,vk(1,jj)-1)*Ld(vk(1,jj),vk(1,jj)-1)
end do
end do

return
end subroutine LdSLd
