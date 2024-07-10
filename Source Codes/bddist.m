function [fnpdf, mn, sd] = bddist(a, b, m, s)
%%
% BDDIST: Function handle for bounded distribution. See Fenton & Griffiths
%         (2008) P60 for more details.
%
% INPUTS:
%
%    a, b:  interval (a, b) of dist
%
%       m:  location parameter. If m=0, the distribution is symmetric about
%           the midpoint of the interval, (a+b)/2.
%
%       s:  scale parameter. The larger s is, the more variable the dist
%           is. s is generally taken 0-5. When s>5, the dist becomes U
%           shaped, like uniform dist. When s gets smaller, it becomes more
%           like normal dist.
%
% OUTPUTS:
%
%   fnpdf:  function handle of bounded dist probability density function
%
%      mn:  mean of this bounded dist
%
%      sd:  standard deviation of this bounded dist; this is a first-order
%           approximation. See Fenton & Griffiths (2008) P60 for details.
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
% Created:      Mar 26, 2024
% Version:      v1.0
% Revision history:
%  none
%


%%
coef1 = @(x) sqrt(pi) * (b-a) ./ sqrt(2) ./ s ./ (x-a) ./ (b-x);
coef2 = 1/2 / s^2;
coef3 = @(x) pi * log((x-a) ./ (b-x));

fnpdf = @(x) coef1(x) .* exp(-coef2 * (coef3(x) - m).^2);

mn = (a+b) / 2;
sd = 0.46 * (b-a) * s / sqrt(4*pi^2 + s^2);     % empirical adjustment of
                                                % the first-order approximation
