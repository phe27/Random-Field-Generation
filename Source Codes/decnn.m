function [k1, k2, m] = decnn(n1, n2, mxm, mxk)
%%
% DECNN: Decompose n1 and n2 to the form of n1 = k1*2**m, n2 = k2*2**m.
%       
%   The decomposition n1 = k1*2**m, n2 = k2*2**m for a common factor 2**m.
%   The integers k1 and k2 are chosen to be as large as possible while
%   requiring the product k1*k2 to be less than or equal to mxk. Note that
%   the direct simulation involves the inversion of a mxk x mxk matrix
%   (at the upper limit) and so mxk should not be overly large. This
%   formulation is somewhat less restrictive than simply requiring n1 and
%   n2 to be powers of 2. Also n1 and n2 do not have to be equal. However
%   n1 and n2 cannot still be chosen arbitrarily, for example the set
%   (n1,n2) = (336,256) results in k1 = 21, k2 = 16, m = 4 which is not
%   acceptable here (since k1*k2 > mxk for mxk = 256), while the set
%   (n1,n2) = (160,256) is acceptable since k1 = 10, k2 = 16, m = 4. In
%   general it may be easier to choose k1, k2, and m before specifying n1
%   and n2. In the event that an unacceptable (k1,k2,m) combination is
%   selected, IERR is set to -1 and control returned to the calling routine.
%   The maximum value of m is set by the calling routine in the argument 
%   mxm.
%
% INPUTS:
%
%   n1, n2: number of cells to discretize the field in the x and y
%           directions respectively (corresponding to the first and second
%           indices of Z respectively). Both n1 and n2 must have the form
%           n1 = k1 * 2**m and n2 = k2 * 2**m where m is common to both and
%           k1 and k2 are positive integers satisfying k1*k2 <= mxk.
%           Generally k1 and k2 are chosen to be as large as possible and
%           still satisfy the above requirements so the the first stage
%           involves directly simulating a k1 x k2 field by inversion of a
%           covariance matrix. A potential example is (n1,n2) = (160,256)
%           which gives k1 = 5, k2 = 8, and m = 5. Note that in general
%           N1 and N2 cannot be chosen arbitrarily - it is usually best to
%           choose m first then k1 and k2 so as to satisfy or exceed the
%           problem requirements. Note that because of the requirements on
%           k1*k2, n1 cannot be more than mxk times as big (or smaller)
%           than n2.
%
%   mxm:    integer giving the largest value that M can take. An error is
%           generated if the process size is such that m > mxm.
%
%   mxk:    max number of k1 x k2 in the initial field
%
% OUTPUTS:
%
%   k1, k2: integers giving the size of the initial field (see C0). It is
%           an error for the product k1*k2 to exceed MXK.
%
%   m:      the number of 2 x 2 subdivisions to perform. It is an error for
%           M to be greater than mxm.
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


k1 = n1;
k2 = n2;

for m = 0 : mxm
    kk = k1 * k2;
    if kk <= mxk
        break;
    end

    j1 = floor(k1/2);
    j2 = floor(k2/2);
    if (2*j1 ~= k1) || (2*j2 ~= k2)
        error(['Cannot determine an acceptable combination of k1, k2, and m! ' ...
               'Try changing N1 and N2.']);
    end

    k1 = j1;
    k2 = j2;
end