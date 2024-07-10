function [nrfx, nrfy] = chknxe(nxe, nye, mxk)
%%
% CHKNXE: Check nxe and nye to ensure compatiblity with LAS2G.
%
%   This routine checks nxe and nye to ensure that they can be written in
%   the form:
%                       nxe = k1*(2**m)
%                       nye = k2*(2**m)
%
%   If not, nxe and nye are adjusted upwards until they do satisfy the
%   above equations for some non-negative integer m and positive
%   integers (k1, k2) with k1*k2 <= mxk.
%
% INPUTS:
%
%   nxe:    integer giving the number of cells in the x-direction in the
%           random fields. nrfx and nrfy should have the following form
%           nx = k1*(2**m), ny = k2*(2**m), where m is a non-negative
%           integer, and k1 is a positive integer with k1*k2 < mxk.
%
%           Generally, k1 and k2 are chosen to be as large as possible and
%           still satisfy the above requirements so the the first stage
%           involves directly simulating a k1 x k2 field by inversion of a
%           covariance matrix. An example is (nx,ny) = (160,256), which
%           gives k1 = 5, k2 = 8, and m = 5. Note that in general nrfx and
%           nrfy cannot be chosen arbitrarily - it is usually best to
%           choose m first then k1 and k2 so as to satisfy or exceed the
%           problem requirements. Note that because of the requirements on
%           k1*k2, nx cannot be more than mxk times as big (or smaller)
%           than ny.
%
%   nye:    see nx above.
%
%   mxk:    maximum value of k1 x k2
%
% OUTPUTS:
%
%   nrfx:   integer giving the number of cells in the x-direction in the
%           random field. This is normally equal to nxe, but if nxe is not
%           an integer of the form nxe = k1*(2**m), where m is a 
%           non-negative integer, and k1 is a positive integer with
%           k1*k2 < mxk then nrfx is the next larger possible integer that
%           does satisfy this requirement of LAS2G.
%
%   nrfy:   integer giving the number of cells in the y-direction in the
%           random field. This is normally equal to nxe, but if nxe is not
%           an integer of the form nye = k2*(2**m), where m is a 
%           non-negative integer, and k1 is a positive integer with
%           k1*k2 < mxk then nrfy is the next larger possible integer that
%           does satisfy this requirement of LAS2G.
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
% compute minimum m
aa = nxe*nye / mxk;
am = log(aa) / (2*log(2));
if am < 0
    am = 0;
end
m  = fix(am + 0.9);                     % round m up
tm = 2^m;

% now compute minimum k's
ak1 = nxe / tm;
k1  = fix(ak1 + 0.9);
ak2 = nye / tm;
k2  = fix(ak2 + 0.8);

% new random field dimensions
nrfx = k1 * fix(tm + 0.5);
nrfy = k2 * fix(tm + 0.5);

% note if we're changing anything
if (nrfx~=nxe) || (nrfy~=nye)
         
    warning(['Incompatible provided number of elements (%d, %d)\n' ...
             '         Adjusting random field size to (%d, %d)\n'], ...
             nxe, nye, nrfx, nrfy);
end

