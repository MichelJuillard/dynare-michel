function [d1,d2] = get_date_of_a_file(filename)
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
    info = dir(filename);
    d1 = info.datenum;
    if nargout>1
        d2 = info.date;
    end