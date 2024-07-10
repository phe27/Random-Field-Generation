function [C0, CT1, CI1, CC, CS, CI, ...
              AT1, AI1, AC, AS, AI, ...
          m, k1, k2] = las2i(nx, ny, xl, yl, mxm, mxk, thx, thy, var)
%%
% LAS2I: Initialize parameters for LAS2G.
%
%   This routine sets up the matrices required by LAS2G to construct
%   realizations of the random field. The general recursive field
%   construction follows the relationship
%
%                   {Z^j} = [A^T]{Z^(j-1)} + [C]{U}
%
%   where {Z^j} is a vector of length 4 representing the values assigned to
%   the 2 x 2 cell subdivision and {Z^(j-1)} are the parent cell values in
%   some neighbourhood. The following figure illustrates the cell
%   subdivision, the neighbourhood, and the numbering scheme for an interior
%   cell (special subsets of the neighbourhood are used for the corners and
%   sides):
%
%                ----------------------------------------
%                |            |            |            |
%                |            |            |            |
%                |     7      |     8      |      9     |
%                |            |            |            |
%                |            |            |            |
%                |------------|------------|------------|
%                |            | 3   |    4 |            |
%                |            |     |      |            |
%                |     4      |-----5------|      6     |
%                |            |     |      |            |
%                |            | 1   |    2 |            |
%                |------------|------------|------------|
%                |            |            |            |
%                |            |            |            |
%                |     1      |     2      |     3      |
%                |            |            |            |
%                |            |            |            |
%                ----------------------------------------
%
%   We see that for the upper left corner, the parent cell neighbourhood
%   used consists of just cells {4,5,7,8} and similarly for the other
%   corners and sides.
%
%   The first stage of the simulation involves the direct generation of a
%   k1 x k2 cell array, where k1 and k2 are integers which satisfying the
%   decomposition nx = k1*2**m, ny = k2*2**m for a common factor 2**m. The
%   integers k1 and k2 are chosen to be as large as possible while
%   requiring the product k1*k2 to be less than or equal to mxk. Note that
%   the direct simulation involves the inversion of a mxk x mxk matrix (at
%   the upper limit) and so mxk should not be overly large. This
%   formulation is somewhat less restrictive than simply requiring nx and
%   ny to be powers of 2. Also nx and ny do not have to be equal. However
%   n1 and n2 cannot still be chosen arbitrarily, for example the set
%   (nx,ny) = (336,256) results in k1 = 21, k2 = 16, m = 4 which is not
%   acceptable here (since k1*k2 > mxk for mxk = 256), while the set
%   (nx,ny) = (160,256) is acceptable since k1 = 10, k2 = 16, m = 4. In
%   general it may be easier to choose k1, k2, and m before specifying nx
%   and ny. The maximum value of m is set by the calling routine in the 
%   argument mxm.
%
% INPUTS:
%
%   nx, ny: number of cells to discretize the field in the x and y
%           directions respectively (corresponding to the first and second
%           indices of Z respectively). Both nx and ny must have the form
%           nx = k1 * 2**m and ny = k2 * 2**m where m is common to both and
%           k1 and k2 are positive integers satisfying k1*k2 <= mxk.
%           Generally k1 and k2 are chosen to be as large as possible and
%           still satisfy the above requirements so the the first stage
%           involves directly simulating a k1 x k2 field by inversion of a
%           covariance matrix. A potential example is (nx,ny) = (160,256)
%           which gives k1 = 5, k2 = 8, and m = 5. Note that in general
%           nx and ny cannot be chosen arbitrarily - it is usually best to
%           choose m first then k1 and k2 so as to satisfy or exceed the
%           problem requirements. Note that because of the requirements on
%           k1*k2, nx cannot be more than mxk times as big (or smaller)
%           than ny.
%
%   xl, yl: physical dimensions of the process.
%
%   mxm:    integer giving the largest value that m can take. An error is
%           generated if the process size is such that m > mxm.
%
%   mxk:    maximum value of k1 x k2
%
%   thx:    correlation length in x-direction
%
%   thy:    correlation length in y-direction
%
%   var:    point variance
%
% OUTPUTS:
%
%   C0:     vector containing the lower triangular matrix of the Cholesky
%           decomposition of the covariance matrix for the initial stage(0)
%           of k1 x k2 cells.
%
%   CT1, CI1:   vectors containing the lower triangular values of the
%               Cholesky decomposition of the covariance matrix for the
%               k1 or k2 = 1 special case.
% 
%   CC:     vector containing the lower triangular values of the Cholesky
%           decomposition of the covariance matrix for the corner cell 2x2
%           subdivisions.
%
%   CS:     vector containing the lower triangular values of the Cholesky
%           decomposition of the covariance matrix for the side cell 2x2
%           subdivisions.
% 
%   CI:     vector containing the lower triangular values of the Cholesky
%           decomposition of the covariance matrix for the interior cell
%           2x2 subdivisions.
%
%   AT1, AI1:   arrays containing the best linear estimation coefficients
%               for the k1 or k2 = 1 special case.
%
%   AC:     array containing the best linear estimation coefficients for 
%           the corner cell subdivisions.
%
%   AS:     array containing the best linear estimation coefficients for 
%           the side cell subdivisions.
%
%   AI:     array containing the best linear estimation coefficients for 
%           the interior cell subdivisions.
%
%   m:      the number of 2 x 2 subdivisions to perform. It is an error for
%           m to be greater than mxm.
%
%   k1, k2: integers giving the size of the initial field (see C0). It is
%           an error for the product k1*k2 to exceed mxk.
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


n = 3;                                  % 4 subdivitions - 1
[k1, k2, m] = decnn(nx, ny, mxm, mxk);  % decompose n1 and n2
fprintf(1, 'LAS2I: k1=%d, k2=%d, m=%d\n', k1, k2, m);


%% form initial covariance matrix and compute its cholesky decomp
% kseed = iseed(kseed);                   % initialize internal generator

t1 = xl / k1;                           % cell size in x-direction
t2 = yl / k2;                           % cell size in y-direction
[Q, R] = dcvit2(thx, thy, var, t1, t2, k1, k2);
[Ltmp, flag] = chol(Q);
C0 = Ltmp';                             % lower triangular
if flag
    error(['LAS2I: Cholesky decomposition of stage 0 covariance matrix is ' ...
           'NOT symmetric!']);
end


%% setup for subsequent stages
AT1 = zeros(2, n);
AI1 = zeros(3, n);
CT1 = zeros(n, n);
CI1 = zeros(n, n);
nm = 1;

% if the following is not excuted, AT1, AI1, CT1, CI1 will all be zero
if (k1 == 1 || k2 == 1) && m > 0
    % special k1, k2 = 1 case
    t1 = 0.5 * t1;
    t2 = 0.5 * t2;

    % get basic cov matrices
    [R, B, S] = dcvmt2(thx, thy, var, t1, t2);
    i2 = 5;
    if k1 == 1
        i1 = 2;
        i3 = 8;
    else
        i1 = 4;
        i3 = 6;
    end
    [AT1, AI1, CT1, CI1] = thin1d(R, B, S, i1, i2, i3);
    nm = 2;
end

AC = zeros(4, n, 4, m);                         % if m=0, AC=[]
CC = zeros(n, n, 4, m);                         % if m=1, last dim omitted
AS = zeros(6, n, 4, m);
CS = zeros(n, n, 4, m);
AI = zeros(9, n, m);
CI = zeros(n, n, m);

% if m=0, the following will not be excuted, and AC, CC, AS, CS, AI, CI
% will all be empty
for kk = nm : m
    t1 = 0.5 * t1;
    t2 = 0.5 * t2;
    [R, B, S] = dcvmt2(thx, thy, var, t1, t2);  % get basic cov matrices

    % corner, side, interior parameters
    [ACtmp, CCtmp] = corn2d(R, B, S);
    [AStmp, CStmp] = side2d(R, B, S);
    [AItmp, CItmp] = intr2d(R, B, S);

    if m==1
        AC = ACtmp;
        CC = CCtmp;
        AS = AStmp;
        CS = CStmp;
        AI = AItmp;
        CI = CItmp;
    else                                        % m>1
        AC(:, :, :, kk) = ACtmp;
        CC(:, :, :, kk) = CCtmp;
        AS(:, :, :, kk) = AStmp;
        CS(:, :, :, kk) = CStmp;
        AI(:, :, kk)    = AItmp;
        CI(:, :, kk)    = CItmp;
    end
end


