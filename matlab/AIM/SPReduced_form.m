% [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,b,condn);
%
% Compute reduced-form coefficient matrix, b.

function [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,condn);
b=[];
%qs=SPSparse(q);
qs=sparse(q);
left = 1:qcols-qrows;
right = qcols-qrows+1:qcols;
nonsing = rcond(full(qs(:,right))) > condn;
if(nonsing)
    qs(:,left) = -qs(:,right)\qs(:,left);
    b = qs(1:neq,1:bcols);
    b = full(b);
else  %rescale by dividing row by maximal qr element
    %'inverse condition number small, rescaling '
    themax=max(abs(qs(:,right)),[],2);
    oneover = diag(1 ./ themax);
    nonsing = rcond(full(oneover *qs(:,right))) > condn;
    if(nonsing)
        qs(:,left) = -(oneover*qs(:,right))\(oneover*(qs(:,left)));
        b = qs(1:neq,1:bcols);
        b = full(b);
    end
end

