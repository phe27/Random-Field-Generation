function [AI, CI] = intr2d(R, B, S)
%%
% INTR2D: Create parameter matrices for interior cell subdivisions used by LAS2G.
%
% INPUTS:
%
%   R:      Matrix R of size 9*9, generated in 'dcvmt2'
%
%   B:      Matrix B of size 4*4, generated in 'dcvmt2'
%
%   S:      Matrix S of size 9*4, generated in 'dcvmt2'
%
% OUTPUTS:
%
%   AI:     Matrix AI of size 9*3, by solving R*AI=S
%
%   CI:     Matrix CI of size 3*3, which is a lower triangular, satisfying
%           CI*CI'=BB-S'*AI
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

% Solve for A
S2 = S(:, 1:n);
AI  = R \ S2;

% Update B (dimension reduced)
B2 = B(1:n, 1:n);
BB = B2 - S2' * AI;

% Cholesky Decomposition
[Ltmp, flag] = chol(BB);
CI = Ltmp';
if flag
    warning('INTR2D: Cholesky decomposition of matrix BB is NOT symmetric!');
end


