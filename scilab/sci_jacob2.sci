function [stk,txt,top]=sci_jacob2()
RHS=[]
for k=1:rhs
  RHS=[stk(top)(1),RHS]
  top=top-1
end
top=top+1
stk=list('jacob2'+rhsargs(RHS),'0','?','?','?')
