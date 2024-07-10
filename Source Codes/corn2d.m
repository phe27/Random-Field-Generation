function [AC, CC] = corn2d(R, B, S)
%%
% CORN2D: Create parameter matrices for corner cell subdivisions used by LAS2G.
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
%   AC:     Matrix AC of size 4*3*4, by solving R*AC=S
%
%   CC:     Matrix CC of size 3*3*4, which is a lower triangular, satisfying
%           CC*CC'=BB-S'*AC
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


%% loop each corner
% the 4th row/column of matrix involving subdivisions (i.e., B, and S)
% should be removed as the 4th subdivision will be calculated to preserve
% the upwards average, see P229, Gordon book, for more details.
n  = 3;                         % 4 subdivisions, so only 3 needed here

% bottom left: 5, 6, 8, 9 no. parent cells extracted
% bottom right: 4, 5, 7, 8 no. parent cells extracted
% top left: 2, 3, 5, 6 no. parent cells extracted
% top right: 1, 2, 4, 5 no. parent cells extracted
ind = [5 6 8 9;
       4 5 7 8;
       2 3 5 6;
       1 2 4 5];                % index no. for each side

AC = zeros(4, n, 4);
CC = zeros(n, n, 4);

for ii = 1 : 4
    % Solve for A
    R2 = R(ind(ii,:), ind(ii,:));
    S2 = S(ind(ii,:), 1:n);
    Atmp = R2 \ S2;
    AC(:, :, ii) = Atmp;

    % Update B (dimension reduced)
    B2 = B(1:n, 1:n);
    BB = B2 - S2' * Atmp;

    % Cholesky Decomposition
    [Ltmp, flag] = chol(BB);
    CC(:, :, ii) = Ltmp';
    if flag
        warning(['CORN2D: Cholesky decomposition of matrix BB is ' ...
                 'NOT symmetric! [Index=%d]'], ii);
    end
end
