function z = las2g(nx, ny, C0, CT1, CI1, CC, CS, CI, ...
                               AT1, AI1, AC, AS, AI, m, k1, k2)
%%
% LAS2G: Produce a 2-D quadrant symmetric stationary Gaussian random field
%        using the local subdivision algorithm.
%
%   This routine creates a zero-mean realization of a 2-D random process
%   given its variance function. Each discrete value generated herein
%   represents the local average of a realization of the process over the
%   area dx x dy, where (dx,dy) is the grid spacing of the desired field.
%   The construction of the realization proceeds recursively as follows:
%
%      1) generate a low resolution field of size k1 x k2. If (n1 x n2) is
%         the desired field resolution, then k1 and k2 are determined so
%         that n1 = k1*2**m, n2 = k2*2**m, where 2**m is a common factor
%         and k1*k2 <= mxk. Generally k1 and k2 are maximized where possible.
%         This is refered to subsequently as the Stage 0 generation.
%      2) subdivide the domain m times by dividing each cell into 4 equal
%         parts (2 x 2). In each subdivision, new random values are
%         generated for each new cell. The parent cells of the previous
%         stage are used to obtain best linear estimates of the mean of each
%         new cell so that the spatial correlation is approximated. Also
%         upwards averaging is preserved (i.e., the average of the 4 new
%         values is the same as their parent cell value). Only parent cells
%         in a neighbourhood of 3 x 3 are considered (the approximation to
%         the spatial correlation).
%
%   The output z here contains values for all subdivisions. Each column is
%   the realization of one subdivision. For each realization (column), the
%   numbering of this 2-D domain is as follows (for example, this is the
%   the 3rd subdivision with nx=5 and ny=3):
%
%                  |-------|-------|-------|-------|-------|
%                  |  11   |  12   |  13   |  14   |  15   |
%                  |z(11,3)|z(12,3)|z(13,3)|z(14,3)|z(15,3)|
%                  |-------|-------|-------|-------|-------|
%                  |   6   |   7   |   8   |   9   |  10   |
%                  |z(6,3) |z(7,3) |z(8,3) |z(9,3) |z(10,3)|
%                  |-------|-------|-------|-------|-------|
%                  |   1   |   2   |   3   |   4   |   5   |
%                  |z(1,3) |z(2,3) |z(3,3) |z(4,3) |z(5,3) |
%                  |-------|-------|-------|-------|-------|
%
%   The linear estimation of the mean is accomplished by using the
%   covariance between local averages over each cell, consistent with the
%   goal of producing a local average field. Note that this conditioning
%   process implies that the construction of cells near the edge of the
%   boundary will require the use of values which are, strictly speaking,
%   outside the boundary of the field in question. This is handled by using
%   special reduced neighbourhoods along the boundaries (equivalent to
%   saying that what goes on beyond the boundary has no effect on the
%   process within the boundary).
%
%   Note that this routine sets up a number of parameters on the first call
%   and thus the time required to produce the first realization is
%   substantially greater than on subsequent calls.
%
% INPUTS:
%
%   nx:     integer giving the number of cells in the x-direction in the
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
%   ny:     see nx above.
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
% OUTPUTS:
%
%   z:      matrix of size at least (nx*ny) x (m+1) which on output will
%           contain the random) field. each column corresponds to one
%           subdiviions.
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


%% set parameters
n  = 3;                                 % 4 subdivisions - 1
nn = nx * ny;
kk = k1 * k2;
z  = zeros(nn, m+1);                    % each column is one stage
zz = zeros(4, 1);                       % store 4 subdivisions temporarily


%% generate realizations
% stage 0 field
U = randn(kk, 1);                       % kk*1 standard normal random numbers
z(1:kk, 1) = C0 * U;                    % values of inital k1*k2 cells

% generate stage 1, 2, ..., m sub-fields
nm  = 1;
k1m = k1;
k2m = k2;

for ii = nm : m
    % extract A and C for this stage `ii`
    if m==1
        ACm = AC;
        CCm = CC;
        ASm = AS;
        CSm = CS;
        AIm = AI;
        CIm = CI;
    else
        ACm = AC(:, :, :, ii);
        CCm = CC(:, :, :, ii);
        ASm = AS(:, :, :, ii);
        CSm = CS(:, :, :, ii);
        AIm = AI(:, :, ii);
        CIm = CI(:, :, ii);
    end

    k1m = 2 * k1m;                      % new field sizes
    k2m = 2 * k2m;


    %% 4 corners
    % corner #1: bottom left
    U = randn(n, 1);
    ic   = 1;                           % corner 1 index
    pind = [1, 2, k1m/2+1, k1m/2+2];    % indexes of 4 parent cells
    zp      = z(pind, ii);              % 4 parent cell values
    zz(1:3) = (ACm(:,:,ic))' * zp + CCm(:,:,ic) * U;
    zz(4)   = 4*zp(ic) - sum(zz(1:3));

    z(1, ii+1)     = zz(1);             % map local (zz) to global (z)
    z(2, ii+1)     = zz(2);
    z(k1m+1, ii+1) = zz(3);
    z(k1m+2, ii+1) = zz(4);

    % corner #2: bottom right
    U = randn(n, 1);
    ic   = 2;                           % corner 2 index
    pind = [k1m/2-1, k1m/2, ...
            k1m-1,   k1m];              % indexes of 4 parent cells
    zp      = z(pind, ii);              % 4 parent cell values
    zz(1:3) = (ACm(:,:,ic))' * zp + CCm(:,:,ic) * U;
    zz(4)   = 4*zp(ic) - sum(zz(1:3));

    z(k1m-1, ii+1)   = zz(1);           % map local (zz) to global (z)
    z(k1m, ii+1)     = zz(2);
    z(2*k1m-1, ii+1) = zz(3);
    z(2*k1m, ii+1)   = zz(4);

    % corner #3: top left
    U = randn(n, 1);
    ic   = 3;                           % corner 3 index
    pind = [1+(k2m/2-2)*k1m/2, ...
            1+(k2m/2-2)*k1m/2+1, ...
            1+(k2m/2-1)*k1m/2, ...
            1+(k2m/2-1)*k1m/2+1];       % indexes of 4 parent cells
    zp      = z(pind, ii);              % 4 parent cell values
    zz(1:3) = (ACm(:,:,ic))' * zp + CCm(:,:,ic) * U;
    zz(4)   = 4*zp(ic) - sum(zz(1:3));

    z(1+(k2m-2)*k1m, ii+1)   = zz(1);   % map local (zz) to global (z)
    z(1+(k2m-2)*k1m+1, ii+1) = zz(2);
    z(1+(k2m-1)*k1m, ii+1)   = zz(3);
    z(1+(k2m-1)*k1m+1, ii+1) = zz(4);

    % corner #4: top right
    U = randn(n, 1);
    ic   = 4;                           % corner 4 index
    pind = [(k2m/2-1)*k1m/2-1, ...
            (k2m/2-1)*k1m/2, ...
            (k2m/2-0)*k1m/2-1, ...
            (k2m/2-0)*k1m/2];           % indexes of 4 parent cells
    zp      = z(pind, ii);              % 4 parent cell values
    zz(1:3) = (ACm(:,:,ic))' * zp + CCm(:,:,ic) * U;
    zz(4)   = 4*zp(ic) - sum(zz(1:3));

    z((k2m-1)*k1m-1, ii+1) = zz(1);     % map local (zz) to global (z)
    z((k2m-1)*k1m, ii+1)   = zz(2);
    z((k2m-0)*k1m-1, ii+1) = zz(3);
    z((k2m-0)*k1m, ii+1)   = zz(4);

    %% 4 sides
    % side #1: bottom
    is = 1;                                 % side 1 index
    for jj = 2 : k1m/2-1
        U = randn(n, 1);
        pind = [jj-1,     jj,         jj+1, ...
                jj+k1m/2-1, jj+k1m/2, jj+k1m/2+1];  % 6 parent cell indexes
        zp      = z(pind, ii);              % 6 parent cell values
        zz(1:3) = (ASm(:,:,is))' * zp + CSm(:,:,is) * U;
        zz(4)   = 4*zp(is+1) - sum(zz(1:3));

        z(2*(jj-1)+1, ii+1)     = zz(1);    % map local (zz) to global (z)
        z(2*(jj-1)+2, ii+1)     = zz(2);
        z(2*(jj-1)+1+k1m, ii+1) = zz(3);
        z(2*(jj-1)+2+k1m, ii+1) = zz(4);
    end

    % side #2: left
    is = 2;                                 % side 2 index
    for jj = k1m/2+1 : k1m/2 : (k2m/2-2)*k1m/2+1
        U = randn(n, 1);
        pind = [jj-k1m/2, jj-k1m/2+1, ...
                jj,       jj+1, ...
                jj+k1m/2, jj+k1m/2+1];      % 6 parent cell indexes
        zp      = z(pind, ii);              % 6 parent cell values
        zz(1:3) = (ASm(:,:,is))' * zp + CSm(:,:,is) * U;
        zz(4)   = 4*zp(is+1) - sum(zz(1:3));

        z(4*(jj-1)+1, ii+1)     = zz(1);    % map local (zz) to global (z)
        z(4*(jj-1)+2, ii+1)     = zz(2);
        z(4*(jj-1)+1+k1m, ii+1) = zz(3);
        z(4*(jj-1)+2+k1m, ii+1) = zz(4);
    end

    % side #3: right
    is = 3;                                 % side 3 index
    for jj = 2*k1m/2 : k1m/2 : (k2m/2-1)*k1m/2
        U = randn(n, 1);
        pind = [jj-k1m/2-1, jj-k1m/2, ...
                jj-1,       jj, ...
                jj+k1m/2-1, jj+k1m/2];      % 6 parent cell indexes
        zp      = z(pind, ii);              % 6 parent cell values
        zz(1:3) = (ASm(:,:,is))' * zp + CSm(:,:,is) * U;
        zz(4)   = 4*zp(is+1) - sum(zz(1:3));

        z(4*(jj-k1m/2)+2*(k1m/2-1)+1, ii+1)     = zz(1);    % local to global
        z(4*(jj-k1m/2)+2*(k1m/2-1)+2, ii+1)     = zz(2);
        z(4*(jj-k1m/2)+2*(k1m/2-1)+1+k1m, ii+1) = zz(3);
        z(4*(jj-k1m/2)+2*(k1m/2-1)+2+k1m, ii+1) = zz(4);
    end

    % side #4: top
    is = 4;                                 % side 4 index
    for jj = 2+(k2m/2-1)*k1m/2 : k1m/2*k2m/2-1
        U = randn(n, 1);
        pind = [jj-k1m/2-1, jj-k1m/2, jj-k1m/2+1, ...
                jj-1,       jj,       jj+1];    % 6 parent cell indexes
        zp   = z(pind, ii);                 % 6 parent cell values
        zz(1:3) = (ASm(:,:,is))' * zp + CSm(:,:,is) * U;
        zz(4)   = 4*zp(is+1) - sum(zz(1:3));

        ctol = k1m/2*(k2m/2-1);             % total parent cell no except last row
        z(4*ctol+2*(jj-1-ctol)+1, ii+1)     = zz(1);    % local to global
        z(4*ctol+2*(jj-1-ctol)+2, ii+1)     = zz(2);
        z(4*ctol+2*(jj-1-ctol)+1+k1m, ii+1) = zz(3);
        z(4*ctol+2*(jj-1-ctol)+2+k1m, ii+1) = zz(4);
    end

    %% interior cells
    for jj = 2 : k1m/2-1
        for kk = 2 : k2m/2-1
            U = randn(n, 1);
            pcind = (kk-1)*k1m/2 + jj;      % central parent cell index
            pind  = [pcind-k1m/2-1, pcind-k1m/2, pcind-k1m/2+1, ...
                     pcind-1,       pcind,       pcind+1, ...
                     pcind+k1m/2-1, pcind+k1m/2, pcind+k1m/2+1];    % 9 parent cell ind
            zp      = z(pind, ii);          % 9 parent cell values
            zz(1:3) = AIm' * zp + CIm * U;
            zz(4)   = 4*zp(5) - sum(zz(1:3));

            ctol = k1m/2*(kk-1);            % total parent cell no before central cell row
            z(4*ctol+2*(jj-1)+1, ii+1)     = zz(1);     % local to global
            z(4*ctol+2*(jj-1)+2, ii+1)     = zz(2);
            z(4*ctol+2*(jj-1)+1+k1m, ii+1) = zz(3);
            z(4*ctol+2*(jj-1)+2+k1m, ii+1) = zz(4);
        end
    end

end



