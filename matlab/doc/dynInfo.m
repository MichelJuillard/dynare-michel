function dynInfo(fun)
[pathstr, name, ext] = fileparts(which(fun));
if strcmp(ext(2:end),'m')
    block = get_internal_doc_block(name,pathstr);
    if ~isempty(block)
        fid = fopen([fun '.texi'],'wt');
        for i=1:size(block,1)
            fprintf(fid,'%s\n',deblank(block(i,:)));
        end
        fclose(fid);
        disp(' ')
        disp(' ')
        system(['makeinfo --plaintext --no-split --no-validate ' fun '.texi']);
        delete([fun '.texi']);
    else
        disp('No documentation for this routine!')
    end
else
    disp('Not a known matlab/octave routine!')
end