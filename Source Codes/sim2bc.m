function [simrf1, simrf2] = sim2bc(rf1, rf2, rc, thx, thy, dx, dy, nrfx, nrfy, ...
                                   nxe, nye, mxm, mxk, lfst, lpltrf, irf, nrf)
%%
% SIM2BC: Simulate two random fields, which are possibly intercorrelated.
% 
%   Individual fields are generated as standard Gaussian random fields
%   using the 2-D Local Average Subdivision (LAS) algorithm. These fields
%   are then combined to produce correlated standard Gaussian fields using
%   the Covariance Matrix Decomposition approach. Finally, the individual
%   fields are transformed so that they have the desired marginal
%   distributions. These transformations are as follows:
%
%        P(x,y) = mean + sd*G(x,y)                  if normally dist
%        P(x,y) = exp{ log-mean + log-sd*G(x,y) }   if lognormally dist
%        P(x,y) = a + 0.5*(b-a)*[ 1 + tanh((m + s*G(x,y))/2*pi) ]
%                                                   if bounded dist
%   where P(x,y) is the desired random field; G(x,y) is one of the standard
%   correlated Gaussian fields.
%
%   If the field is deterministic, the entire field is simply set to its mean.
%
%   Only the last subdivision (the end column of `simrf`) will be taken
%   here, which is only a single column, say `z` here. The numbering of 
%   this original 2-D random field (`z`) is below (here nxfx=5 and nyfy=3, 
%   thus `z` is of size 15*1):
%
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         | z(11) | z(12) | z(13) | z(14) | z(15) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |         Original rf
%         | z(6)  | z(7)  | z(8)  | z(9)  | z(10) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         | z(1)  | z(2)  | z(3)  | z(4)  | z(5)  |
%         |-------|-------|-------|-------|-------|
%
%   To make the random field consistent with matrix layout (`z(1,1)` is the
%   top left corner), the original rf is reshaped to size of `(nrfy,nrfx)`
%   (example below is `z(1:3,1:5)`) as below.
%
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         |z(1,1) |z(1,2) |z(1,3) |z(1,4) |z(1,5) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |         New rf
%         |z(2,1) |z(2,2) |z(2,3) |z(2,4) |z(2,5) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         |z(3,1) |z(3,2) |z(3,3) |z(3,4) |z(3,5) |
%         |-------|-------|-------|-------|-------|
%
%   This can be done using: ` z = reshape(z, [nrfx,nrfy]); z = flipud(z') `
%
% INPUTS:
%
%   rf1:    real vector of length at least 7 containing the parameters of
%           random field 1. Notably,
%              rf1(1) = mean,
%              rf1(2) = standard deviation,
%              rf1(3) = distribution type;
%                       = 0.0 if rf1 is deterministic (at mean value)
%                       = 1.0 if rf1 is normally distributed
%                       = 2.0 if rf1 is lognormally distributed
%                       = 3.0 if rf1 is bounded
%              rf1(4) = lower bound (if bounded)
%              rf1(5) = upper bound (if bounded)
%              rf1(6) = m parameter (if bounded)
%              rf1(7) = s parameter (if bounded)
%           If rf1 is bounded, rf1(1) and rf1(2) are ignored and the values
%           rf1(4) through rf1(7) completely describe the distribution.
%
%   rf2:    real vector of length at least 7 containing the parameters of
%           random field 2. See `rf1` for what the various elements of rf2
%           contain.
%
%   rc:     upper triangular matrix, which is the Cholesky decomposition of
%           the correlation coefficient matrix between the two random
%           fields. Unit matrix if not correlated.
%
%   varfnc: character string containing the name of the variance function
%           controlling the random fields. Possible variance functions are:
%           `dlavx2` - 2-D exponentially decaying (Markov) model requires
%                      X- and Y-direction scales of fluctuation.
%           `dlafr2` - 2-D isotropic fractional Gaussian noise model
%                      requires (H,delta) as parameters. In this case, thx
%                      is H, and delta is the minimum element dimension.
%           `dlsep2` - 2-D separable (1D x 1D) Markov model requires X- and
%                      Y-direction scales of fluctuation.
%           `dlsfr2` - 2-D separable fractional Gaussian noise model
%                      requires (H_x,H_y,delta) as parameters. In this case,
%                      thx is H_x, thy is H_y, and delta is the minimum
%                      element dimension.
%           `dlspx2` - 2-D separable Gaussian decaying model requires X-
%                      and Y-direction scales of fluctuation.
%
%   thx:    real value giving the x-direction scale of fluctuation (or,
%           at least, this is a parameter of the variance function).
%
%   thy:    real value giving the y-direction scale of fluctuation (or,
%           at least, this is a parameter of the variance function).
%
%   nxe:    integer giving the number of elements describing the random
%           fields needed for analysis in the x-direction (horizontally).
%
%   nye:    integer giving the number of elements describing the random
%           fields needed for analysis in the y-direction (horizontally).
%
%   nrfx:   integer giving the number of cells in the x-direction in the
%           random fields rf1 and rf2. This is normally equal to nxe, but
%           if nxe is not an integer of the form nxe = k1*(2**m), where m
%           is a non-negative integer, and k1 is a positive integer with
%           k1*k2 < mxk (see chknxe) then nrfx is the next larger possible
%           integer that does satisfy this requirement of LAS2G.
%
%   nrfy:   integer giving the number of cells in the y-direction in the
%           random fields rf1 and rf2. This is normally equal to nxe, but
%           if nxe is not an integer of the form nye = k2*(2**m), where m
%           is a non-negative integer, and k2 is a positive integer with
%           k1*k2 < mxk (see chknxe) then nrfy is the next larger possible
%           integer that does satisfy this requirement of LAS2G.
%
%   dx:     real value giving the physical size of an element in the
%           x-direction.
%
%   dy:     real value giving the physical size of an element in the
%           y-direction.
%
%   mxm:    integer giving the largest value that M can take. An error is
%           generated if the process size is such that m > mxm.
%
%   mxk:    maximum value of k1 x k2
%
%   lfst:   logical flag that will be true only when this is the first time
%           to enter sim2bc.
%
%   lpltrf: flag if a random field will be displayed
%
%   irf:    if a random field will be displayed, which one? 1 -> rf1;
%           2 -> rf2
%
%   nrf:    how many random fields? 1 or 2
%
% OUTPUTS:
%
%   simrf1: real matrix of size nrfx x nrfy which on output will contain
%           the (optionally random) field 1.
%
%   simrf2: real matrix of size nrfx x nrfy which on output will contain
%           the (optionally random) field 2. if nrf==1, simrf2=[].
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
% Version:      v1.3
% Revision history:
%  1.1  reshape `simrf` here to match matrix layout in matlab. See above
%       descriptions. Replaced some double for loops. (Mar 22/2024)
%  1.2  now plotted random fields are actual ones, rather than standard
%       normal ones previously. (Mar 26/2024)
%  1.3	added title for random field plot. (Jun 20/2024)
%


%%
global C0 CT1 CI1 CC CS CI AT1 AI1 AC AS AI m k1 k2;

xl = nrfx * dx;                                 % entire physical field size
yl = nrfy * dy;


%% initialize las2i
if lfst                                         % run only once
    var = 1;                                    % std norm first
    [C0, CT1, CI1, CC, CS, CI, ...
         AT1, AI1, AC, AS, AI, ...
     m, k1, k2] = las2i(nrfx, nrfy, xl, yl, mxm, mxk, thx, thy, var);
end


%% generate standard normal
simrf1 = las2g(nrfx, nrfy, C0, CT1, CI1, CC, CS, CI, ...
                               AT1, AI1, AC, AS, AI, m, k1, k2);
simrf1 = reshape(simrf1(:,end), [nrfx,nrfy]);
simrf1 = flipud(simrf1');                       % now match matrix layout

if nrf==2
    simrf2 = las2g(nrfx, nrfy, C0, CT1, CI1, CC, CS, CI, ...
                                   AT1, AI1, AC, AS, AI, m, k1, k2);
    simrf2 = reshape(simrf2(:,end), [nrfx,nrfy]);
    simrf2 = flipud(simrf2');
elseif nrf==1
    simrf2 = [];
end


%% combine for correlated fields if 2 random fields
if ~isempty(rc) && nrf==2                       % cross-correlated
    simrf2 = rc(1,2)*simrf1 + rc(2,2)*simrf2;
end

% % plot random field? ONLY standard normal here
% if lpltrf
%     if irf==1
%         simrfa = simrf1(1:nye, 1:nxe);          % actual rf size for analysis
%     elseif irf==2 && nrf==2
%         simrfa = simrf2(1:nye, 1:nxe);
%     else
%         error('SIM2BC: Random field number incorrect! irf must be 1 or 2.\n');
%     end
% 
%     im = pltfld(simrfa, nxe, nye, dx, dy);      % plot random field
% end


%% convert to final fields
% random field #1
if rf1(3)==0                                    % deterministic
    simrf1 = ones(size(simrf1)) * rf1(1);
elseif rf1(3)==1			                    % normal
    simrf1 = rf1(1) + rf1(2)*simrf1;
elseif rf1(3)==2                                % lognormal
    rf1(5) = log(1 + (rf1(2))^2/(rf1(1))^2);    % var of log-rf1
    rf1(4) = log(rf1(1)) - 0.5*rf1(5);          % mean of log-rf1
    rf1(5) = sqrt(rf1(5));                      % sd of log-rf1
    simrf1 = exp( rf1(4) + rf1(5)*simrf1 );
elseif rf1(3)==3                                % bounded
    twopi  = 2*pi;
    simrf1 = rf1(4) + 0.5*(rf1(5)-rf1(4)) * ...
             (1 + tanh((rf1(6)+rf1(7)*simrf1)/twopi));
else
    error('SIM2BC: Unknown distribution type of Random field 1! Check!');
end

% random field #2
if nrf==2
    if rf2(3)==0                                    % deterministic
        simrf2 = ones(size(simrf2)) * rf2(1);
    elseif rf2(3)==1			                    % normal
        simrf2 = rf2(1) + rf2(2)*simrf2;
    elseif rf2(3)==2                                % lognormal
        rf2(5) = log(1 + (rf2(2))^2/(rf2(1))^2);    % var of log-rf1
        rf2(4) = log(rf2(1)) - 0.5*rf2(5);          % mean of log-rf1
        rf2(5) = sqrt(rf2(5));                      % sd of log-rf1
        simrf2 = exp(rf2(4) + rf2(5)*simrf2);
    elseif rf2(3)==3                                % bounded
        twopi  = 2*pi;
        simrf2 = rf2(4) + 0.5*(rf2(5)-rf2(4)) * ...
                 (1 + tanh((rf2(6)+rf2(7)*simrf2)/twopi));
    else
        error('SIM2BC: Unknown distribution type of Random field 2! Check!');
    end
end

% plot random field? Final fields
if lpltrf
    if irf==1
        simrfa = simrf1(1:nye, 1:nxe);          % actual rf size for analysis
		im = pltfld(simrfa, nxe, nye, dx, dy, ...
					'Random Field, $rf1$');		% plot random field
    elseif irf==2 && nrf==2
        simrfa = simrf2(1:nye, 1:nxe);
		im = pltfld(simrfa, nxe, nye, dx, dy, ...
					'Random Field, $rf2$');		% plot random field
    else
        error('SIM2BC: Random field number incorrect! irf must be 1 or 2.\n');
    end

    
end
