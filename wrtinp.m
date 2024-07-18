function wrtinp(outfid, rf1, rf2, thx, thy, crcoef, dx, dy, nxe, nye, kseed, nsim)
%%
% WRTINP: Write input parameters to *.out file.
% 
% INPUTS:
%
%   rf1:    vector contains 7 numbers, giving the statistics of random
%           field 1. See the descriptions above.
%
%   rf2:    vector contains 7 numbers, giving the statistics of random
%           field 2. See the descriptions above.
%
%   thx:    correlation length in x-direction.
% 
%   thy:    correlation length in y-direction.
%
%   crcoef: cross correlation coefficient between random fields #1 and #2.
%
%   dx:     physical element size in x-direction.
% 
%   dy:     physical element size in y-direction.
%
%   nxe:    number of elements in x-direction.
% 
%   nye:    number of elements in y-direction.
%
%   nsim:   total number of simulations.
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
% Created:      Jul 18, 2024
% Version:      v2.0
% Revision history:
%  2.0	adapted for rf geneation only. (Jul 18/24)
%


%%
% some conversions
if rf1(3)==0
	rf1s = 'Deterministic';
elseif rf1(3)==1
	rf1s = 'Normal';
elseif rf1(3)==2
	rf1s = 'Lognormal';
elseif rf1(3)==3
	rf1s = 'Bounded';
end

if rf2(3)==0
	rf2s = 'Deterministic';
elseif rf2(3)==1
	rf2s = 'Normal';
elseif rf2(3)==2
	rf2s = 'Lognormal';
elseif rf2(3)==3
	rf2s = 'Bounded';
end

% write inputs
fprintf(outfid, "Input Parameter:\n");
fprintf(outfid, "------------------------------------------------------" + ...
				"------------\n");
fprintf(outfid, "Distribution of random field 1 . . . . . . . . %s\n", rf1s);
fprintf(outfid, "Mean and SD of random field 1  . . . . . . . . %5.3e %5.3e\n", ...
				 rf1(1), rf1(2));
fprintf(outfid, "Distribution of random field 2 . . . . . . . . %s\n", rf2s);
fprintf(outfid, "Mean and SD of random field 2  . . . . . . . . %5.3e %5.3e\n", ...
				 rf2(1), rf2(2));
fprintf(outfid, "Correlation lengths in x- and y-dir  . . . . . %5.3e %5.3e\n", ...
				 thx, thy);
fprintf(outfid, "Cross-correlation between the two RFs  . . . . %5.3e\n", ...
				 crcoef);
fprintf(outfid, "Element sizes in x- and y-dir  . . . . . . . . %5.3e %5.3e\n", ...
				 dx, dy);
fprintf(outfid, "Number of elements in x- and y-dir . . . . . . %d %d\n", ...
				 nxe, nye);
fprintf(outfid, "Initial seed number  . . . . . . . . . . . . . %d\n", ...
				 kseed);
fprintf(outfid, "Total number of realizations   . . . . . . . . %d\n", ...
				 nsim);


