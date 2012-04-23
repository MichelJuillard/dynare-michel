function [mean,variance] = GetPosteriorMeanVariance(M,drop)
    
    DirectoryName = CheckPath('metropolis',M.dname);
    o=load([ DirectoryName '/'  M.fname '_mh_history']);
    NbrDraws = sum(o.record.MhDraws(:,1));
    NbrFiles = sum(o.record.MhDraws(:,2));
    NbrBlocks = o.record.Nblck;
    mean = 0;
    variance = 0;
    z = [];
    
    fname1 = [DirectoryName '/' M.fname];
    nkept = 0;
    for i=1:NbrBlocks
        n = 0;
        for j=1:NbrFiles
            o = load([fname1 '_mh' int2str(j) '_blck' int2str(i)]);
            m = size(o.x2,1);
            if n + m < drop*NbrDraws
                n = n + m;
                continue
            elseif n < drop*NbrDraws
                k = drop*NbrDraws - n + 1;
                x2 = o.x2(k:end,:);
            else
                x2 = o.x2;
            end
            z =[z; x2];        
            p = size(x2,1);
            mean = (nkept*mean + sum(x2)')/(nkept+p);
            x = bsxfun(@minus,x2,mean');
            variance = (nkept*variance + x'*x)/(nkept+p);
            n = n + m;
            nkept = nkept + p;
        end
    end
    
    