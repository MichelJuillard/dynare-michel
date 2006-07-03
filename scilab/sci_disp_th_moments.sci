function [stk,txt,top]=sci_disp_th_moments()
RHS=[]
for k=1:rhs
  RHS=[stk(top)(1),RHS]
  top=top-1
end
top=top+1
stk=list('disp_th_moments'+rhsargs(RHS),'0','?','?','?')
