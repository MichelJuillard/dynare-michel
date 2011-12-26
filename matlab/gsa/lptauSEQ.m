function [lpmat] = lptauSEQ(Nsam,Nvar)

% [lpmat] = lptauSEQ(Nsam,Nvar)
%
% Generates a sample from a sobol sequence of length Nsam for a
% number of parameters Nvar
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


clear lptau
lpmat = zeros(Nsam, Nvar);
seed = int64(0);
for j=1:Nsam,
    [v, seed] = qmc_sequence(Nvar, seed, 0);
    lpmat(j,:) = v;
end
return
