function plot_ms_variance_decomposition(M_, options_, vd, title_, graph_save_formats, ...
                                        TeX, varargin)
% plot the variance decomposition of shocks
%
% Inputs
%       M_
%		shocks: matrix of the individual shocks Tx(KxK)with J=number of shocks
%
% Optional Inputs
%		'data': the actual data, TxK with K=number of data series
%		'steady': the steady state value, TxK
%		'shock_names': to specify the names of the shocks
%		'series_names': to specify the names of the different series
%		'dates': pass a date vector to use, otherwise will just index on 1:T
%		'colors': Jx3 list of the rgb colors to use for each shock
%
% Example:
% plot_historic_decomposition(shocks,'VD','shock_names',shock_names,'series_names',series_names)

% Copyright (C) 2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

    nvars = M_.endo_nbr;
    endo_names = M_.endo_names;

    names = {};
    tex_names = {};
    m = 1;
    for i=1:M_.orig_endo_nbr
        tex_name = deblank(M_.endo_names_tex(i,:));
        if ~isempty(tex_name)
            names{m} = deblank(endo_names(i,:));
            tex_names{m} = tex_name;
            m = m + 1;
        end
    end
    
    dims = size(vd);
    if length(dims) == 3
  	[T,K,J] = dims;
  	shocks = vd;
    else
        T = dims(1);
        K = nvars;
        J = nvars;
        temp_vd = zeros(T,K,J);
        for i=1:nvars
            for j=1:nvars
                temp_vd(:,i,j) = vd(:,((j-1) + ((i-1)*nvars)+1));
            end
        end
        shocks = temp_vd;
    end
    
    for i=1:nvars
        shock_names{i} = endo_names(i,:);
        series_names{i} = endo_names(i,:);
    end

    if nargin < 2
        title_ = '';
    end

    x = [1:T]; plot_dates = 0;
    
    data = 0;
    steady = 0;
    
    colors = [ .1 .1 .75
               .8 0 0
               1 .7 .25
               1 1 0	
               .5 1 .5 
               .7 .7 .1
               .5 .6 .2
               .1 .5 .1];
    
    % overide the defaults with optional inputs
    for i=1:length(varargin)
        if strcmpi(varargin{i},'data')
            data = varargin{i+1};
        elseif strcmpi(varargin{i},'steady')
            steady = varargin{i+1};
        elseif strcmpi(varargin{i},'shock_names')
            shock_names = varargin{i+1};
        elseif strcmpi(varargin{i},'series_names')
            series_names = varargin{i+1};
        elseif strcmpi(varargin{i}, 'dates')
            x = varargin{i+1}; plot_dates = 1;
        elseif strcmpi(varargin{i},'colors')
            colors = varargin{i+1};
        end
    end
    
    % add an extra period to the time series
    x(T+1) = x(T) + (x(T) - x(T-1));

    
    figure('Name',title_)
    for k=1:K
        % Go through each series
        subplot(K,1,k);
        sshocks = shocks(:,k,:);
        hold on
        % plot the stacked shocks
        for t=1:T
            % step through each time period
            pos_position = 0; neg_position = 0;
            xt = [x(t) x(t) x(t+1) x(t+1)];
            for j=1:J
                % stack each shock
                st = sshocks(t,1,j);
                if st < 0
                    yi = st+neg_position;
                    y = [neg_position yi yi neg_position];
                    neg_position = yi;
                else
                    yi = st+pos_position;
                    y = [pos_position yi yi pos_position];
                    pos_position = yi;
                end
                fill(xt,y,colors(j,:));
                XY(t,j,:) = y;
            end
        end
        if data
            plot(x(2:end)',data(:,k),'k','LineWidth',2.5);
        end
        if steady
            plot(x(2:end)',steady(:,k), '--k','LineWidth',2.25);
        end
        if k==K
	    if exist('OCTAVE_VERSION')
                legend(shock_names,'Location','SouthOutside');
	    else 
                legend(shock_names,'Location','BestOutside','Orientation','horizontal'); 
            end
        end
        
        hold off
        if plot_dates
            datetick 'x';
        end
        xlim([min(x),max(x)])
        ylim([0 , 1])
        grid on
        title(series_names{k});
        %suptitle(title_);
    end
    dyn_save_graph([options_.ms.output_file_tag filesep 'Output' ...
        filesep 'Variance_Decomposition'], 'MS-Variance-Decomposition', ...
        graph_save_formats, TeX,names,tex_names,'Variance decomposition');
end
