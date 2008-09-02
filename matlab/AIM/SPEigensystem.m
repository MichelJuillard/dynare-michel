%  [w,rts,lgroots] = SPEigensystem(a,uprbnd)
%
%  Compute the roots and the left eigenvectors of the companion
%  matrix, sort the roots from large-to-small, and sort the
%  eigenvectors conformably.  Map the eigenvectors into the real
%  domain. Count the roots bigger than uprbnd.

function [w,rts,lgroots,flag_trouble] = SPEigensystem(a,uprbnd,rowsLeft) 
opts.disp=0; 
try
    [w,d]   = eigs(a',rowsLeft,'LM',opts);
    rts     = diag(d);
    mag     = abs(rts);
    [mag,k] = sort(-mag);
    rts     = rts(k);
catch
    %disp('Catch in SPE');
    %pause(0.5);
    %aStr=datestr(clock);
    %eval(['save ' regexprep(aStr,' ','')  ' a']);
    try
        [w,d]=eig(a');
    catch
        lasterr
        w=[];rts=[];lgroots=[];
        flag_trouble=1;
        return
    end
    rts     = diag(d);
    mag     = abs(rts);
    [mag,k] = sort(-mag);
    rts     = rts(k);
end
flag_trouble=0; 

%ws=SPSparse(w);
ws=sparse(w);
ws       = ws(:,k);

%  Given a complex conjugate pair of vectors W = [w1,w2], there is a
%  nonsingular matrix D such that W*D = real(W) + imag(W).  That is to
%  say, W and real(W)+imag(W) span the same subspace, which is all
%  that aim cares about. 

ws = real(ws) + imag(ws);

lgroots = sum(abs(rts) > uprbnd);

w=full(ws);