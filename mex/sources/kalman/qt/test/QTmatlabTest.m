size=100;
AA=rand(size);
BB=rand(size);
[QT1 QT2 U1 U2]=qz(AA,BB);

a=rand(size,1);

Loops=10000

% Calling QT without use of Sylv Vector and General Matrix
t = clock;  
[AAAA]=qtamvm(QT1,a,Loops);
QTcpp_noSylv_TaInnerLoop_time=etime(clock, t)

% Calling QT using of Sylv Vector and General Matrix
t = clock;  
[AAAA]=qtmvm(QT1,a,Loops);
QTcppTaInnerLoop_time=etime(clock, t)

t = clock;  
[AAAA]=qtmvm_sub(QT1,a,Loops);
QTcppTaInnerLoop_time=etime(clock, t)

t = clock;  
for tt=1:Loops%0
[AAAA]=qtmvm_sub(QT1,a,1);
end
QTcppTaOuterLoop_time=etime(clock, t)

t = clock;  
[AAAA]=gmvm(QT1,a,Loops);
GMcppTaInnrLoop_time=etime(clock, t)

t = clock;  
for tt=1:Loops%0
[AAAA]=gmvm(QT1,a,1);
end
GMcppTaOuterLoop_time=etime(clock, t)

t = clock;  
for tt=1:Loops%0
MTA=QT1*a;
end
matlabTa_time=etime(clock, t)

