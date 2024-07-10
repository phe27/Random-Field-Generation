function [Q, R] = dcvit2(thx, thy, var, dx, dy, k1, k2)
%%
% DCVIT2: Compute initial k1*k2 x k1*k2 covariance matrix and a 9 x 9
%         submatrix between local averages. Used by LAS2G.
%
%   This routine forms the covariance matrix Q between local averages
%   arranged on a k1 x k2 grid using a user provided variance function. The
%   grid and element numbering schemes appear as follows (for k1=k2=5)
%
%              | dx |
%         ---  --------------------------
%          dy  | 21 | 22 | 23 | 24 | 25 |
%         ---  --------------------------
%              | 16 | 17 | 18 | 19 | 20 |
%              --------------------------           ----------------
%              | 11 | 12 | 13 | 14 | 15 |           |  7 |  8 |  9 |
%              --------------------------           ----------------
%              |  6 |  7 |  8 |  9 | 10 |           |  4 |  5 |  6 |
%              --------------------------           ----------------
%              |  1 |  2 |  3 |  4 |  5 |           |  1 |  2 |  3 |
%              --------------------------           ----------------
%
%                      Q Array                          R Array
%
%   with the k1 x k2 array numbering on the left and a basic 3 x 3 subarray
%   shown on the right. If we call Z_i the local average of a random process
%   over the i'th cell, then for E[Z] = 0, the elements of Q are defined by
%
%                       Q(i,j) = E[Z_i*Z_j]
%
%   Note that when calculating covariance, assume that bottom left conner
%   is the origin (0,0), and positive x-dir is to the right, and positive
%   y-dir is upward.
%
% INPUTS:
%
%   thx:    correlation length in x-dir
%
%   thy:    correlation length in y-dir
%
%   var:    point variance
%
%   dx:     x-dimension of the local average cells
%
%   dy:     y-dimension of the local average cells
%
%   k1:     number of local average cells in x-dir
%
%   k2:     number of local average cells in y-dir
%
% OUTPUT:
%
%   Q:      array of size at least k1*k2 x k1*k2 which on output will
%           contain the covariance matrix between the local average
%           elements shown above on the left.
%
%   R:      array of size at least 9 x 9 which on output will contain the
%           covariance matrix between the local average elements shown
%           above on the right.
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


%% find Q
ntol = k1 * k2;
Q = zeros(ntol, ntol);

% only lower part here
for ii = 1 : ntol
    nox1 = mod(ii, k1); nox1(nox1==0) = k1;     % make sure nox1<=k1
    ax = [nox1-1, nox1] * dx;                   % x limits of the ii'th cell
    noy1 = ceil(ii/k1);
    ay = [noy1-1, noy1] * dy;                   % y limits of the ii'th cell
    for jj = 1 : ii
        nox2 = mod(jj, k1); nox2(nox2==0) = k1; % make sure nox2<=k1
        bx = [nox2-1, nox2] * dx;               % x limits of the jj'th cell
        noy2 = ceil(jj/k1);
        by = [noy2-1, noy2] * dy;               % y limits of the jj'th cell

        Q(ii, jj) = dcvab2(thx, thy, var, ax, ay, bx, by);
        Q(jj, ii) = Q(ii, jj);                  % symmetric
    end
end


%% find R
ntol = 3 * 3;
R = zeros(ntol, ntol);

% only lower part here
for ii = 1 : ntol
    nox1 = mod(ii, 3); nox1(nox1==0) = 3;       % make sure nox1<=3
    ax = [nox1-1, nox1] * dx;                   % x limits of the ii'th cell
    noy1 = ceil(ii/3);
    ay = [noy1-1, noy1] * dy;                   % y limits of the ii'th cell
    for jj = 1 : ii
        nox2 = mod(jj, 3); nox2(nox2==0) = 3;   % make sure nox2<=3
        bx = [nox2-1, nox2] * dx;               % x limits of the jj'th cell
        noy2 = ceil(jj/3);
        by = [noy2-1, noy2] * dy;               % y limits of the jj'th cell

        R(ii, jj) = dcvab2(thx, thy, var, ax, ay, bx, by);
        R(jj, ii) = R(ii, jj);                  % symmetric
    end
end
