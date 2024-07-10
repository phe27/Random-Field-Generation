function fun = corrfun(thx, thy)
%%
% CORRFUN: Define a correlation function, Markov here.
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

fun = @(t1,t2,eta1,eta2) ...
      1 + expm1(-2*sqrt( ((t1-eta1)/thx).^2 + ((t2-eta2)/thy).^2) );

