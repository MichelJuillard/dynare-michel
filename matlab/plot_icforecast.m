function plot_icforecast(Variables)
% Build plots for the conditional forecasts.
%
% INPUTS 
%  o Variables     [char]        m*x array holding the names of the endogenous variables to be plotted. 
%
% OUTPUTS
%  None.
% 
% SPECIAL REQUIREMENTS
%  This routine has to be called after imcforecast.m.

% Copyright (C) 2006-2009 Dynare Team
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

load conditional_forecasts;

for i=1:size(Variables,1)
    eval(['ci1 = forecasts.cond.ci.' Variables(i,:) ';'])
    eval(['m1 = forecasts.cond.mean.' Variables(i,:) ';'])
    eval(['ci2 = forecasts.uncond.ci.' Variables(i,:) ';'])
    eval(['m2 = forecasts.uncond.mean.' Variables(i,:) ';'])
    build_figure(Variables(i,:),ci1,ci2,m1,m2);
end

function build_figure(name,ci1,ci2,m1,m2)
    figure('Name',['Conditional forecast: ' name '.']);
    H = length(m1);
    h1 = area(1:H,ci1(2,1:H));
    set(h1,'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))
    set(h1,'FaceColor',[.9 .9 .9])
    hold on
    h2 = area(1:H,ci1(1,1:H));
    set(h2,'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))
    set(h2,'FaceColor',[1 1 1])
    plot(1:H,m1,'-k','linewidth',3)
    plot(1:H,m2,'--k','linewidth',3)
    plot(1:H,ci2(1,:),'--k','linewidth',1)
    plot(1:H,ci2(2,:),'--k','linewidth',1)
    axis tight
    hold off