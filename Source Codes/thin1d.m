function [AT, AI, CT, CI] = thin1d(R, B, S, i1, i2, i3)
%%
% THIN1D: Create parameter matrices for interior cell subdivisions of LAS2G
%         for the case k1 or k2=1 and for LAS3G for k1=k2=1, k1=k3=1, or
%         k2=k3=1.
%
% INPUTS:
%
%   R:      Matrix R of size 9*9, generated in 'dcvmt2'
%
%   B:      Matrix B of size 4*4, generated in 'dcvmt2'
%
%   S:      Matrix S of size 9*4, generated in 'dcvmt2'
%
%   i1, i2, i3: indexes into the global covariance matrix R to be extracted
%               to form the local covariance matrix.
%
% OUTPUTS:
%
%   AT:     Matrix AT of 2x3, by solving RT*AT=ST
%
%   AI:     Matrix AI of 3x3, by solving RI*AI=SI
%
%   CT:     Matrix CT of 3x3, which is a lower triangular, satisfying
%           CT*CT'=BB-ST'*AT
%
%   CI:     Matrix CI of 3x3, which is a lower triangular, satisfying
%           CI*CI'=BB-SI'*AI
%
% -------------------------------------------------------------------------
% Author:       Dr Pengpeng He
% Organisation: University of Dundee
% Email:        phe001@dundee.ac.uk
% Website:
%   <a href="matlab: 
%   web('http://discovery.dundee.ac.uk/en/persons/pengpeng-he')">Author's Site</a>
% -------------------------------------------------------------------------
%
% Created:      Feb 16, 2024
% Version:      v1.0
% Revision history:
%  none
%


%%
% the 4th row/column of matrix involving subdivisions (i.e., B, and S)
% should be removed as the 4th subdivision will be calculated to preserve
% the upwards average, see P229, Gordon book, for more details.
n  = 3;                         % 4 subdivisions, so only 3 needed here

% extract R
RT(1, 1) = R(i2, i2);                   % RT (edges)
RT(1, 2) = R(i2, i3);
RT(2, 2) = R(i3, i3);
RT(2, 1) = RT(1, 2);

RI(1, 1) = R(i1, i1);                   % RI (interior)
RI(1, 2) = R(i1, i2);
RI(1, 3) = R(i1, i3);
RI(2, 1) = RI(1, 2);
RI(2, 2) = R(i2, i2);
RI(2, 3) = R(i2, i3);
RI(3, 1) = RI(1, 3);
RI(3, 2) = RI(2, 3);
RI(3, 3) = R(i3, i3);

% Solve for A
ST = S([i2,i3], 1:n);
AT  = RT \ ST;

SI = S([i1,i2,i3], 1:n);
AI  = RI \ SI;

% Update B (dimension reduced)
B2 = B(1:n, 1:n);
BBT = B2 - ST' * AT;
BBI = B2 - SI' * AI;

% Cholesky Decomposition
[CT, flagT] = chol(BBT);
CT = CT';
if flagT
    warning('thin1d: Cholesky decomposition of matrix BBT is NOT symmetric!');
end

[CI, flagI] = chol(BBI);
CI = CI';
if flagI
    warning('thin1d: Cholesky decomposition of matrix BBI is NOT symmetric!');
end


