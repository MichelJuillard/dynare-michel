function [stk,txt,top]=sci_jacob_a()
RHS=[]
for k=1:rhs
  RHS=[stk(top)(1),RHS]
  top=top-1
end
top=top+1
stk=list('jacob_a'+rhsargs(RHS),'0','?','?','1')
