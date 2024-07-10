function [rfmn, rfstd] = chkfld(simrf, idfle)
%%
% CHKFLD: Compute mean and sd of the simulated random field and print.
%
% INPUTS:
%
%   simrf:  simulated random field
%
%   idfle:  ID of the file where to print
%
% OUTPUTS:
%
%   rfmn:   mean of the simulated random field.
%
%   rfstd:  standard deviation of the simulated random field.
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
% Version:      v1.1
% Revision history:
%  1.1  added file ID to output (Mar 22/2024)
%


%%
rfmn  = mean(simrf(:));
rfstd = std(simrf(:));

fprintf(idfle, 'RF MEAN = %8.3f, SD = %8.3f\n', rfmn, rfstd);
