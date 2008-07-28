function  [q,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
                       SPQZ(aa,bb)
%  [q,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
%                       SPQZ(aa,bb)
%
%  Solve a linear perfect foresight model using the matlab eig
%  function to find the invariant subspace associated with the big
%  roots.  This procedure will fail if the companion matrix is
%  defective and does not have a linearly independent set of
%  eigenvectors associated with the big roots.
% 
%  Input arguments:
% 
%    h         Structural coefficient matrix (neq,neq*(nlag+1+nlead)).
%    neq       Number of equations.
%    nlag      Number of lags.
%    nlead     Number of leads.
%    condn     Zero tolerance used as a condition number test
%              by numeric_shift and reduced_form.
%    uprbnd    Inclusive upper bound for the modulus of roots
%              allowed in the reduced form.
% 
%  Output arguments:
% 
%    b         Reduced form coefficient matrix (neq,neq*nlag).
%    rts       Roots returned by eig.
%    ia        Dimension of companion matrix (number of non-trivial
%              elements in rts).
%    nexact    Number of exact shiftrights.
%    nnumeric  Number of numeric shiftrights.
%    lgroots   Number of roots greater in modulus than uprbnd.
%    aimcode     Return code: see function aimerr.
[neq,aCols]=size(aa);
[bRows,bCols]=size(bb);
if(~and(neq ==aCols,and(neq== bRows,neq ==bCols)))
error(SPQZ:'aa, bb must be square and same dim')
end
nlag=0
nlead=1
b=[];
rts=[];
ia=[];
nexact=[];
nnumeric=[];
lgroots=[];
aimcode=[];

condn=1.0e-8
uprbnd=1.0e16

% Initialization.
nexact   = 0;
nnumeric = 0;
lgroots  = 0;
iq       = 0;
aimcode    = 0;
b=0;
qrows = neq*nlead;
qcols = neq*(nlag+nlead);
bcols = neq*nlag;
q        = zeros(qrows,qcols);
rts      = zeros(qcols,1);

% Compute the auxiliary initial conditions and store them in q.
h=[aa bb];


[h,q,iq,nexact] = SPExact_shift(h,q,iq,qrows,qcols,neq);
   if (iq>qrows) 
      aimcode = 61;
      return;
   end

[h,q,iq,nnumeric] = SPNumeric_shift(h,q,iq,qrows,qcols,neq,condn);
   if (iq>qrows) 
      aimcode = 62;
      return;
   end

%  Build the companion matrix.  Compute the stability conditions, and
%  combine them with the auxiliary initial conditions in q.  

[a,ia,js] = SPBuild_a(h,qcols,neq);

if (ia ~= 0)
   [w,rts,lgroots] = SPEigQZ(a,uprbnd,min(length(ia),qrows-iq+1));
%enough in w to fill q
   q = SPCopy_w(q,w,js,iq,qrows);
end


