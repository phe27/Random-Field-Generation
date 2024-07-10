function [R, B, S] = dcvmt2(thx, thy, var, dx, dy)
%%
% DCVMT2: Compute covariance matrix between local averages used by LAS2G.
%
%   This routine computes the three covariance matrices, R(9x9), B(4x4),
%   and S(9x4), required by LAS2G. The matrices R and B represent the
%   covariances between a 3 x 3 array of equal sized local average elements
%   and a 2 x 2 subarray respectively, numbered as follows;
%
%            |-- dx--|   R Array                         B Array
%      ---   -------------------------           -------------------------
%       |    |       |       |       |           |       |       |       |
%       dy   |   7   |   8   |   9   |           |       |       |       |
%       |    |       |       |       |           |       |       |       |
%      ---   -------------------------           -------------------------
%            |       |       |       |           |       |       |       |
%            |   4   |   5   |   6   |           |   3   |   4   |       |
%            |       |       |       |           |       |       |       |
%            -------------------------           -------------------------
%            |       |       |       |           |       |       |       |
%            |   1   |   2   |   3   |           |   1   |   2   |       |
%            |       |       |       |           |       |       |       |
%            -------------------------           -------------------------
%
%                                       Figure 1
%
%   Notice that the array B is just equal to selected elements of the array
%   R (ie. B(1,1) = R(1,1), B(1,4) = R(1,5), etc.). S represents the 9x4
%   covariance matrix between a doubly large array and an interior subdivided
%   2x2 array as shown, where the numbering of the larger array is the same
%   as on the left in Figure 1.
%
%              |--2dx--|
%        ---   -------------------------
%         |    |       |       |       |
%        2dy   |       |       |       |
%         |    |       |       |       |
%        ---   -------------------------
%              |       | 3 | 4 |       |
%              |       |-------|       |      Figure 2
%              |       | 1 | 2 |       |
%              -------------------------
%              |       |       |       |
%              |       |       |       |
%              |       |       |       |
%              -------------------------
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
%   dx:     x-dimension of the subdivided elements. The dimension of the
%           parent elements is assumed to be 2*dx.
%
%   dy:     y-dimension of the subdivided elements. The dimension of the
%           parent elements is assumed to be 2*dy.
%
% OUTPUT:
%
%   R:      array of size at least 9 x 9 which on output will contain the
%           covariance matrix between the local average elements shown
%           above in Figure 1, left.
%
%   B:      array of size at least 4 x 4 which on output will contain the
%           covariance matrix between the local average elements shown
%           above in Figure 1, right.
%
%   S:      array of size at least 9 x 4 which on output will contain the
%           covariances between the set of subdivided central elements and
%           the parent 3 x 3 set of doubly large elements (see Figure 2)
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


B = zeros(4, 4);
R = zeros(9, 9);
S = zeros(9, 4);

%% find B
B(1,1) = dcvab2(thx, thy, var, [0,1]*dx, [0,1]*dy, [0,1]*dx, [0,1]*dy);
B(1,2) = dcvab2(thx, thy, var, [0,1]*dx, [0,1]*dy, [1,2]*dx, [0,1]*dy);
B(1,3) = dcvab2(thx, thy, var, [0,1]*dx, [0,1]*dy, [0,1]*dx, [1,2]*dy);
B(1,4) = dcvab2(thx, thy, var, [0,1]*dx, [0,1]*dy, [1,2]*dx, [1,2]*dy);

B(2,1) = B(1,2);
B(2,2) = B(1,1);
B(2,3) = B(1,4);
B(2,4) = B(1,3);

B(3,1) = B(1,3);
B(3,2) = B(2,3);
B(3,3) = B(1,1);
B(3,4) = B(1,2);

B(4,1) = B(1,4);
B(4,2) = B(2,4);
B(4,3) = B(3,4);
B(4,4) = B(1,1);


%% find R
fac = 2;
R(1,1) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [0,1]*fac*dx, [0,1]*fac*dy);
R(1,2) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [1,2]*fac*dx, [0,1]*fac*dy);
R(1,3) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [2,3]*fac*dx, [0,1]*fac*dy);
R(1,4) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [0,1]*fac*dx, [1,2]*fac*dy);
R(1,5) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [1,2]*fac*dx, [1,2]*fac*dy);
R(1,6) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [2,3]*fac*dx, [1,2]*fac*dy);
R(1,7) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [0,1]*fac*dx, [2,3]*fac*dy);
R(1,8) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [1,2]*fac*dx, [2,3]*fac*dy);
R(1,9) = dcvab2(thx, thy, var, [0,1]*fac*dx, [0,1]*fac*dy, [2,3]*fac*dx, [2,3]*fac*dy);

R(2,1) = R(1,2);
R(2,2) = R(1,1);
R(2,3) = R(1,2);
R(2,4) = R(1,5);
R(2,5) = R(1,4);
R(2,6) = R(1,5);
R(2,7) = R(1,8);
R(2,8) = R(1,7);
R(2,9) = R(1,8);

R(3,1) = R(1,3);
R(3,2) = R(2,3);
R(3,3) = R(1,1);
R(3,4) = R(1,6);
R(3,5) = R(1,5);
R(3,6) = R(1,4);
R(3,7) = R(1,9);
R(3,8) = R(1,8);
R(3,9) = R(1,7);

R(4,1) = R(1,4);
R(4,2) = R(2,4);
R(4,3) = R(3,4);
R(4,4) = R(1,1);
R(4,5) = R(1,2);
R(4,6) = R(1,3);
R(4,7) = R(1,4);
R(4,8) = R(1,5);
R(4,9) = R(1,6);

R(5,1) = R(1,5);
R(5,2) = R(2,5);
R(5,3) = R(3,5);
R(5,4) = R(4,5);
R(5,5) = R(1,1);
R(5,6) = R(1,2);
R(5,7) = R(1,5);
R(5,8) = R(1,4);
R(5,9) = R(1,5);

R(6,1) = R(1,6);
R(6,2) = R(2,6);
R(6,3) = R(3,6);
R(6,4) = R(4,6);
R(6,5) = R(5,6);
R(6,6) = R(1,1);
R(6,7) = R(1,6);
R(6,8) = R(1,5);
R(6,9) = R(1,4);

R(7,1) = R(1,7);
R(7,2) = R(2,7);
R(7,3) = R(3,7);
R(7,4) = R(4,7);
R(7,5) = R(5,7);
R(7,6) = R(6,7);
R(7,7) = R(1,1);
R(7,8) = R(1,2);
R(7,9) = R(1,3);

R(8,1) = R(1,8);
R(8,2) = R(2,8);
R(8,3) = R(3,8);
R(8,4) = R(4,8);
R(8,5) = R(5,8);
R(8,6) = R(6,8);
R(8,7) = R(7,8);
R(8,8) = R(1,1);
R(8,9) = R(1,2);

R(9,1) = R(1,9);
R(9,2) = R(2,9);
R(9,3) = R(3,9);
R(9,4) = R(4,9);
R(9,5) = R(5,9);
R(9,6) = R(6,9);
R(9,7) = R(7,9);
R(9,8) = R(8,9);
R(9,9) = R(1,1);


%% find S
S(1,1) = dcvab2(thx, thy, var, [0,2]*dx, [0,2]*dy, [2,3]*dx, [2,3]*dy);
S(1,2) = dcvab2(thx, thy, var, [0,2]*dx, [0,2]*dy, [3,4]*dx, [2,3]*dy);
S(1,3) = dcvab2(thx, thy, var, [0,2]*dx, [0,2]*dy, [2,3]*dx, [3,4]*dy);
S(1,4) = dcvab2(thx, thy, var, [0,2]*dx, [0,2]*dy, [3,4]*dx, [3,4]*dy);

S(2,1) = dcvab2(thx, thy, var, [2,4]*dx, [0,2]*dy, [2,3]*dx, [2,3]*dy);
S(2,2) = S(2,1);
S(2,3) = dcvab2(thx, thy, var, [2,4]*dx, [0,2]*dy, [2,3]*dx, [3,4]*dy);
S(2,4) = S(2,3);

S(3,1) = S(1,2);
% S(3,1) = dcvab2(thx, thy, var, [4,6]*Dx, [0,2]*Dy, [2,3]*Dx, [2,3]*Dy);
S(3,2) = S(1,1);
% S(3,2) = dcvab2(thx, thy, var, [4,6]*Dx, [0,2]*Dy, [3,4]*Dx, [2,3]*Dy);
S(3,3) = S(1,4);
S(3,4) = S(1,3);

S(4,1) = dcvab2(thx, thy, var, [0,2]*dx, [2,4]*dy, [2,3]*dx, [2,3]*dy);
S(4,2) = dcvab2(thx, thy, var, [0,2]*dx, [2,4]*dy, [3,4]*dx, [2,3]*dy);
S(4,3) = S(4,1);
S(4,4) = S(4,2);

S(5,1) = dcvab2(thx, thy, var, [2,4]*dx, [2,4]*dy, [2,3]*dx, [2,3]*dy);
S(5,2) = S(5,1);
S(5,3) = S(5,1);
S(5,4) = S(5,1);

S(6,1) = S(4,2);
S(6,2) = S(4,1);
S(6,3) = S(6,1);
S(6,4) = S(4,3);

S(7,1) = S(1,3);
S(7,2) = S(1,4);
S(7,3) = S(1,1);
S(7,4) = S(1,2);

S(8,1) = S(2,3);
S(8,2) = S(2,4);
S(8,3) = S(2,1);
S(8,4) = S(2,2);

S(9,1) = S(1,4);
S(9,2) = S(1,3);
S(9,3) = S(1,2);
S(9,4) = S(1,1);

