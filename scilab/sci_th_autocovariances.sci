function [stk,txt,top]=sci_th_autocovariances()
RHS=[]
for k=1:rhs
  RHS=[stk(top)(1),RHS]
  top=top-1
end
top=top+1
stk=list('th_autocovariances'+rhsargs(RHS),'0','?','?','?')
