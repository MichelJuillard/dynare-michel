function build_internal_documentation()
% The name of the function should be explicit...

datafiles  = [];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'initialize_dataset'}];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'descriptive_statistics'}];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'compute_stdv'}];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'compute_cova'}];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'compute_corr'}];
datafiles = [ datafiles ; {'../../matlab/utilities/dataset'}, {'compute_acov'}];

estimationfiles = [];

simulationfiles = [];

miscfiles = [];

% Data.
if exist('data.texi')
    delete('data.texi')
end
fid = fopen('data.texi','w');
if rows(datafiles)
    for i=1:rows(datafiles)
        block = get_internal_doc_block(datafiles{i,2},datafiles{i,1});
        if isempty(block), continue, end
        for j=1:rows(block)
            fprintf(fid,[block(j,:) '\n']);
        end
        fprintf(fid,'\n\n\n');
    end
end
fclose(fid);

% Estimation.
if exist('estimation.texi')
    delete('estimation.texi')
end
fid = fopen('estimation.texi','w');
if rows(estimationfiles)
    for i=1:rows(datafiles)
        block = get_internal_doc_block(datafiles{i,2},datafiles{i,1});
        if isempty(block), continue, end
        for j=1:rows(block)
            fprintf(fid,[block(j,:) '\n']);
        end
        fprintf(fid,'\n\n\n');
    end
end
fclose(fid);

% Simulation.
if exist('simulation.texi')
    delete('simulation.texi')
end
fid = fopen('simulation.texi','w');
if rows(simulationfiles)
    for i=1:rows(datafiles)
        block = get_internal_doc_block(datafiles{i,2},datafiles{i,1});
        if isempty(block), continue, end
        for j=1:rows(block)
            fprintf(fid,[block(j,:) '\n']);
        end
        fprintf(fid,'\n\n\n');
    end
end
fclose(fid);

% Miscellaneous.
if exist('misc.texi')
    delete('misc.texi')
end
fid = fopen('misc.texi','w');
if rows(miscfiles)
    for i=1:rows(datafiles)
        block = get_internal_doc_block(datafiles{i,2},datafiles{i,1});
        if isempty(block), continue, end
        for j=1:rows(block)
            fprintf(fid,[block(j,:) '\n']);
        end
        fprintf(fid,'\n\n\n');
    end
end
fclose(fid);