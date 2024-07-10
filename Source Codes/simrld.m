function rld = simrld(ld, n)
%%
% SIMRLD: Simulate a random load that is lognormally distributed.
%       
%   Statistics and Machine Learning Toolbox is required (`lognrnd`).
%
% INPUTS:
%
%   ld:     load statistics: ld(1): load mean; ld(2): load standard
%           deviation.
%
%   n:      number of simulation.
%
% OUTPUTS:
%
%   rld:    vector containing n simulated loads.
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
mn = ld(1);
sd = ld(2);

sdln = sqrt(log(1 + sd^2/mn^2));                % sd of ln-values
mnln = log(mn) - 0.5*sdln^2;                    % mean of ln-values

rld = lognrnd(mnln, sdln, n, 1);

