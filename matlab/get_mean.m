function y0 = get_mean(varargin);
global oo_ M_

ys_ = oo_.steady_state; 
lgy_ = M_.endo_names;


lgobs_= [];
mfys=[];
for j=1:length(varargin),
    dum = strmatch(varargin{j},lgy_,'exact');
    mfys = [mfys dum];
end

y0 = ys_(mfys);
