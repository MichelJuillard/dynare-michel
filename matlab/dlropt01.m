function [NVR,ALPHA,opts,separ]= dlropt0(y,z,IRW,method,nvr,alpha,nvr0,alpha0,opts,ALG,output,P0,p);
% DLROPT  Hyper-parameter estimation for DLR
%
% [nvr,alpha,opts,parse]=dlropt(y,z,TVP,meth,nvrc,alphac,nvr0,alpha0,opts,ALG,tab,P0)
%                               1 2  3   4    5     6     7     8     9   10  11  12
%
% y: Time series (*)
% z: Regressors (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% meth: Estimation method ('ml')
%     'ml': Maximum Likelihood
%     'f#': Sum of squares of the #-step-ahead forecasting errors
% nvrc: Constraints for each NVR (-2)
%       -2: Free estimation
%       -1: Constrained estimation (all parameters with nvrc=-1 are equal)
%      >=0: NVR constrained to this value (it is not estimated)
% alphac: Constraints for each alpha (-2, -1, or >=0 as for nvrc) (1)
% nvr0: Initial NVR hyper-parameters (0.0001)
% alpha0: Initial alpha hyper-parameters (1)
% opts: Optimisation options. Type 'help foptions' for details
% ALG: Optimisation algorithm: 0=fmins, 1=fminu, 2=leastsq (not ML) (0)
% tab: Display: 0=none, 1=tabulate results, 2=update window (2)
% P0: Initial P matrix (1e5)
%
% nvr: Estimated NVR hyper-parameters
% alpha: Estimated alpha hyper-parameters
% opts: Returned optimisation options. Type 'help foptions' for details
% parse: Standard Error of hyper-parameters (omit to reduce computation time)
%
% Example: dlropt(y, [ones(size(u)) u], 0, [], [0 -2])
%   regression type model y = c1 + c2(t)*u, with an RW model for
%   both parameters and c1 assumed constant (NVR fixed at zero)
%
% See also DLR, FCAST, STAND, FOPTIONS, FMINS

% Copyright (c) 2004 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

y=y(:);  % ensure y is a column, 28/6/99, JT

Z=[y, z];  %23/9/98, JT
[n, cZ]= size(Z); 

p=[];

steps=0;

P0= diag(P0);
zlabel='Z1'; 

nnvr=1;

Miss= isnan(Z(:,1)); 

F= zeros(sum(IRW+1)); j=1;
for i= 1:cZ-1
    if IRW(i), f= [alpha(i) 1; 0 1]; else, f= alpha(i); end 
    F(j:j+IRW(i), j:j+IRW(i))= f; j= j+IRW(i)+1;
end

% Initial conditions ... strange initial condition for no alpha
% optimisation ....   *******************   WT 11/04
%b0= log10([alpha0./(1-alpha0) nvr0]);
b0=log10(nvr0);

t1= clock; 

% Optimisation update windows, JT, revised 3/8/99

% ALG is 0
% Matlab is newer than 5
% output is 2
% steps is 0
if output==2
    update_panel=figure('Units', 'Normalized','Position', [0.425 0.5 0.15 0.15], ...
        'Number', 'Off', 'Name', 'Working...', 'MenuBar', 'None');
    stop_alg=uicontrol(update_panel, 'Style', 'Text', ...
        'Units', 'Normalized', 'Position', [0.1 12/15 0.8 2/15]);
    stop_meth=uicontrol(update_panel, 'Style', 'Text', ...
        'Units', 'Normalized', 'Position', [0.1 9/15 0.8 2/15]);
    stop_info=uicontrol(update_panel, 'Style', 'Text', 'Units', ...
        'Normalized', 'Position', [0.1 6/15 0.8 2/15], 'Tag', 'info');
    stop_but=uicontrol(update_panel, 'Style', 'Push', 'String', 'S T O P', ...
        'Call', 'set(gcbf,''UserData'', 1);', 'Units', 'Normalized', ...
        'Position', [0.1 1/15 0.8 4/15]);
    set(update_panel, 'UserData', 0);  % before stop button is pressed
    set(stop_alg, 'String', 'ALG = fmins')
    set(stop_meth, 'String', 'Log-likelihood')
end
if exist('stop_info')==0, stop_info=[]; end
[lb, opts]= fstop('dlrfun01', b0, opts, [], ...
    Z, nvr, IRW, Miss, F, P0, steps, ALG, p, stop_info);  % JT, 30/01/04, fstop

if exist('update_panel'), close(update_panel), end

separ= nan*ones(size(lb));

separ=seml('dlrfun01', lb, Z, nvr, IRW, Miss, F, P0, steps, ALG, p, 0);

tiempo= etime(clock, t1); 

NVR=lb;

SEPARn=separ;
SEPARa=0;

ALPHA= alpha';

NVR=10.^NVR;

% Display results
if output
    disp(' ')
    if tiempo>60
        t1= [int2str(fix(tiempo/60)) ' minutes,  ' ...
                num2str(fix(rem(tiempo, 60))) ' seconds.']; 
    else, t1= [num2str(tiempo) ' seconds.']; 
    end
    hora= fix(clock); if hora(5)<10, pega= '0'; else, pega= []; end
    disp('TIME DOMAIN ESTIMATION')
    if steps==0, disp('METHOD: MAXIMUM LIKELIHOOD'),    % DP, 10/02/99
    else, disp(['METHOD: SUM OF SQUARES ' int2str(steps) '-STEPS-AHEAD FORECAST ERRORS'])
    end 
    if ALG==0, disp('OPTIMISER: FMINS');  % DP, 10/02/99
    elseif ALG==1, disp('OPTIMISER: FMINU'); 
    elseif ALG==2, disp('OPTIMISER: LEASTSQ'); 
    end 
    disp(['Date: ' date ' / Time: ' num2str(hora(4)) ':' pega num2str(hora(5))])
    disp(t1)
    disp('RW   NVR       -S.E.       +S.E.')
    disp([num2str(IRW,1) '  ' num2str(NVR,4) '  ' ,...
            num2str(10.^(lb-SEPARn),3) '  ' num2str(10.^(lb+SEPARn),3)])
    
    
    if steps==0, disp(['Log-likelihood: ' num2str(-0.5*opts(8)-length(y)/2*(log(2*pi)+1))])
    else, disp(['Sum of Squares of ' int2str(steps) '-steps-ahead forecast errors: ' num2str(opts(8))])
    end
    disp(' ')
end

% end of m-file