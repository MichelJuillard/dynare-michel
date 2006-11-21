function x0 = stab_map_cyc(Period, alpha2, alpha, pprior)
%
% function x0 = stab_map_(Period, alpha2, alpha)
%
% Mapping of target cyclicity region in the prior ranges applying
% Monte Carlo filtering techniques
%
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models
% I. Mapping stability, MIMEO, 2005.
%
% INPUTS
% Period = range of periodicity for the targete cycles [Tmin Tmax]
% alpha2 =  significance level for bivariate sensitivity analysis
% [abs(corrcoef) > alpha2]
% alpha =  magnitude if cyclical eigenvalue w.r.t. max eigenvalue
% (default=0.9)
%
% OUTPUT: 
% x0: one parameter vector for which the model provides targeted cycle.
%
% GRAPHS
% 1) Pdf's of marginal distributions under the ciclicity (dotted
%     lines) and non-cyclicity (solid lines) regions
% 2) Cumulative distributions of: 
%   - cyclic subset (dotted lines) 
%   - uncyclic subset (solid lines)
% 3) Bivariate plots of significant correlation patterns 
%  ( abs(corrcoef) > alpha2) under the cyclic/uncyclic subsets
%
%
% Copyright (C) 2005 Marco Ratto
% THIS PROGRAM WAS WRITTEN FOR MATLAB BY
% Marco Ratto,
% Unit of Econometrics and Statistics AF
% (http://www.jrc.cec.eu.int/uasa/),
% IPSC, Joint Research Centre
% The European Commission,
% TP 361, 21020 ISPRA(VA), ITALY
% marco.ratto@jrc.it 
%
% ALL COPIES MUST BE PROVIDED FREE OF CHARGE AND MUST INCLUDE THIS COPYRIGHT
% NOTICE.
%
% USES: stab_map_1
%       stab_map_2

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

dr_ = oo_.dr;
if isfield(dr_,'ghx'),
    ys_ = oo_.dr.ys;
    nspred = size(dr_.ghx,2);
    nboth = dr_.nboth;
    nfwrd = dr_.nfwrd;
end
fname_ = M_.fname;

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

if nargin==0,
    Nsam=2000; %2^13; %256;
end
if nargin<4,
    pprior=1;
end
if nargin<3,
    alpha=0.9;
end
if nargin<2,
    alpha2=0.3;
end

if pprior,
    load([fname_,'_prior'])
    delete([fname_,'_cyc*.*']);
    delete([fname_,'_cyclic*.*']);
else
    load([fname_,'_mc'])
    delete([fname_,'_cyc_mc*.*']);
    delete([fname_,'_cyclic_mc*.*']);
end

Nsam = size(lpmat,1);    


if length(istable)>0,
    thex=istable;
    for j=1:length(istable),
        PerX{j}=[];
        MaxEig = max(abs(egg(1:nspred,istable(j))));
        ic = find(imag(egg(1:nspred,istable(j))));
        %i=find( abs(egg( ic ,istable(j)) )>0.9*MaxEig); %only consider complex dominant eigenvalues 
        i=find( abs(egg( ic ,istable(j)) )>alpha*MaxEig); %only consider complex dominant eigenvalues 
        if ~isempty(i),
            i=i(1:2:end);
            thedum=[];
            for ii=1:length(i),
                idum = ic( i(ii) );
                thedum(ii)=2*pi/abs(angle(egg(idum,istable(j))));
            end
            if thedum<Period(1) | thedum>Period(2),
                thex(j)=0;
            else
                PerX{j}=[PerX{j} thedum(find(thedum>Period(1) & thedum<Period(2)))];
            end
            %             [dum, icx]=max(thedum);
            %             icy(j) = ic( i(icx) );
            %             PerX(j)=max(thedum);
            %             if PerX(j)<Period(1) | PerX(j)>Period(2),
            %                 thex(j)=0;
            %             end
        else
            thex(j)=0;
        end
    end
    inocyc=istable(find(thex==0));   % non-cyclic params
    icyc=thex(find(thex));   % cyclic params
    
    if ~isempty(icyc)
        if pprior, 
            cycnam='cyc';
            cycnam2='cyclic';
        else
            cycnam='cyc_mc';
            cycnam2='cyclic_mc';
        end
        proba = stab_map_1(lpmat, icyc, inocyc, cycnam);
        
        stab_map_2(lpmat(icyc,:),alpha2, cycnam2);
        x0=0.5.*(bayestopt_.ub(1:nshock)-bayestopt_.lb(1:nshock))+bayestopt_.lb(1:nshock);
        x0 = [x0; lpmat(icyc(1),:)'];
        
    else
        disp('None of the parameter values gives the targeted cyclicity!')
        x0=[];
    end
    %stab_map_2(lpmat(iunstable,:),alpha2, 0);
    
else
    if length(iunstable)==0,
        disp('All parameter values in the specified ranges are stable!')
    else
        disp('All parameter values in the specified ranges are unstable!')        
    end
    
end


