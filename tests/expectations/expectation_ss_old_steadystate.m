function [ys_, check_] = expectation_ss_old_steadystate(ys_orig_, exo_)
    ys_=zeros(6,1);
    global M_
    ys_(4)=0;
    ys_(6)=0;
    ys_(5)=0.3333333333333333;
    ys_(3)=((1/M_.params(1)-(1-M_.params(4)))/(M_.params(3)*ys_(5)^(1-M_.params(3))))^(1/(M_.params(3)-1));
    ys_(1)=ys_(5)^(1-M_.params(3))*ys_(3)^M_.params(3);
    ys_(2)=ys_(1)-M_.params(4)*ys_(3);
    M_.params(5)=(1-M_.params(3))*ys_(1)/(ys_(2)*ys_(5)^(1+M_.params(6)));
    check_=0;
end
