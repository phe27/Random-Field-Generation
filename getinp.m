function [rf1, rf2, thx, thy, crcoef, dx, dy, nxe, nye, kseed, nsim, ...
		  lplt, ino, irf, ldbg] = getinp(inpfnm)
%%
% GETINP: Read input parameters from an input data file.
% 
%   It reads the numbers at the end of each line from column 60. Each line
%   can contain multiple numbers. The format is:
%
% =========================================================================
%   1.  mean and standard deviation of random field #1  . . . rf1(1)... [rf1(4)... rf1(7)]
%   2.  mean and standard deviation of random field #2  . . . rf2(1)... [rf2(4)... rf2(7)]
%   3.  correlation lengths of random fields  . . . . . . . . thx  thy
%   4.  cross-correlation between the two fields  . . . . . . crcoef
%   5.  element sizes in x- and y-direction . . . . . . . . . dx  dy
%   6.  number of elements in x- and y-direction  . . . . . . nxe nye
%   7.  initial seed number . . . . . . . . . . . . . . . . . kseed
%   8.  total number of realizations  . . . . . . . . . . . . nsim
%   9.  plot a random field?  . . . . . . . . . . . . . . . . lplt ino  irf
%   10. debug on? . . . . . . . . . . . . . . . . . . . . . . ldbg
% =========================================================================
% 
%   An example is as follows:
% =========================================================================
%   1.  mean and standard deviation of random field #1  . . . 30.0 9.0 2 0 0 0 0
%   2.  mean and standard deviation of random field #2  . . . 30.0 9.0 2 0 0 0 0
%   3.  correlation lengths of random fields  . . . . . . . . 30.0 30.0
%   4.  cross-correlation between the two fields  . . . . . . 0.5
%   5.  element sizes in x- and y-direction . . . . . . . . . 0.15 0.15
%   6.  number of elements in x- and y-direction  . . . . . . 150  50
%   7.  initial seed number . . . . . . . . . . . . . . . . . 111
%   8.  total number of realizations  . . . . . . . . . . . . 10
%   9.  plot a random field?  . . . . . . . . . . . . . . . . 1 1 1
%   10. debug on? . . . . . . . . . . . . . . . . . . . . . . 1
% =========================================================================
%
%
% DESCRIPTIONS:
%
%   line 1: contains at least 7 numbers (rf1) giving the statistics of
%           random fields #1. Notably,
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
%   line 2: contains at least 7 numbers (rf2) giving the statistics of
%           random fields #2. See line 2 for more info.
%
%   line 3: contains 2 numbers (thx, thy) giving the correlation lengths in
%           x- and y-directions.
%
%   line 4: contains 1 number (crcoef), giving the cross correlation
%           between random fields #1 and #2.
%
%   line 5: contains 2 integers (dx, dy), giving the physical element sizes
%           in x- and y-direction.
%
%   line 6: contains 2 integers (nxe, nye), indicating the number of
%           elements in x- and y-direction.
%
%   line 7: contains 1 integer (kseed), giving the initial seed number for
%           the random number generator.
%
%   line 8: contains 1 integer (nsim), giving the total number of
%           realizations in the analysis.
%
%   line 9: contains 3 numbers (lplt, ino, irf), giving if a random field
%           will be plotted. the 1st one should be 1 or 0 to indicate
%			plot (1) or not (0). If yes, `ino` is a postive integer giving 
%           which realization will be plotted, and `irf` gives which random 
%           field to be plotted. `irf` should be 1 or 2. If `lplt` is 0, no 
%			plot is plotted, and `ino` and `irf` will be ignored.
%
%  line 10: contains 1 integer (idbg), giving if debug mode is needed. If
%           on, `ldbg` should be set to be 1; otherwise 0.
%
%
% INPUTS:
%
%   inpfnm: input data file name.
%
%
% OUTPUTS:
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
%   kseed:  initial seed number for the random number generator.
%
%   nsim:   total number of simulations.
%
%	lplt:	plot a random field?
%
%   ino:    which no. of simulation is showing the random field? If 0, not
%           no plot will be shown.
%
%   irf:    if ino~=0, this specifies which random field will be shown.
%
%   ldbg:   debug mode on or not?
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
%  1.1	changed Line 11 to 3 numbers with the first one indicating plotting
%		or not. (Jun 06/24)
%  2.0	adapted for rf geneation only. (Jul 18/24)


%%
rcol = 59;                                  % read from this column

if exist(inpfnm, 'file') ~= 2
    error('GETINP: Input file does not exist!\n');
else
    % extract lines from cell array
    fid = fopen(inpfnm, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    lines = lines{1};
    fclose(fid);
end


%% for each line...
for ii = 1 : numel(lines)
    % extract numbers starting from column #rcol
    line  = lines{ii};
    numbs = sscanf(line(rcol:end), '%f');
    if isempty(numbs)
        error('GETINP: No data extracted from line %d! Please check!\n', ii);
    end

    % append the numbers to variable
    switch ii
        case 1
            rf1 = numbs;
        case 2
            rf2 = numbs;
        case 3
            thx = numbs(1);
            thy = numbs(2);
        case 4
            crcoef = numbs;
        case 5
            dx = numbs(1);
            dy = numbs(2);
        case 6
            nxe = numbs(1);
            nye = numbs(2);
        case 7
            kseed = numbs;
        case 8
            nsim = numbs;
        case 9
			lplt = numbs(1);
            ino  = numbs(2);
            irf  = numbs(3);
        case 10
            ldbg = numbs;
    end

end

