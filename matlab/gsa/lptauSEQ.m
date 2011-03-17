function [lpmat] = lptauSEQ(Nsam,Nvar)

% [lpmat] = lptauSEQ(Nsam,Nvar)
%
% function call LPTAU and generates a sample of dimension Nsam for a
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
for j=1:Nsam,
    lpmat(j,:)=LPTAU(j,Nvar);
end
return
