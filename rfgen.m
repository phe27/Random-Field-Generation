function [simrf1, simrf2] = rfgen(inpfnm)
%%
% RFGEN: Main program to simulate two-dimensional random fields that can be
%		 cross-correlated using the Local Average Subdivision (LAS) method.
%
%   In this programme, at most two random fields (`simrf1` & `simrf2`) can
%	be considered, and can be correlated if necessary. Refer to Fenton &
%	Griffiths's book ('Risk assessment in geotechnical engineering'),
%	Chapter 3. 
%
%   The generated rf from `sim2bc` follows matlab matrix layout with the
%   size of nrfy X nrfx ( the example below is `z(1:3,1:5)`):
%
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         |z(1,1) |z(1,2) |z(1,3) |z(1,4) |z(1,5) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         |z(2,1) |z(2,2) |z(2,3) |z(2,4) |z(2,5) |
%         |-------|-------|-------|-------|-------|
%         |       |       |       |       |       |
%         |z(3,1) |z(3,2) |z(3,3) |z(3,4) |z(3,5) |
%         |-------|-------|-------|-------|-------|
%
% INPUTS:
%
%   inpfnm: input data file name
%
% OUTPUTS:
%
%   simrf1: simulated random field #1. size of nye*nxe*nsim
%
%   simrf2: simulated random field #2. size of nye*nxe*nsim. If only 1 
%           random field is specified in the input data file, simrf==[]
%
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
%  1.1	added `lplt` to indicate plot a random field or not. (Jun/06/24)
%  1.2	fixed a 'lplt' bug. (Jun/07/24)
%  2.0	adapted for rf geneation only. (Jul 18/24)
%


%%
mxk = 256;								% max of k1*k2
mxm = 10;                               % max of m
nrf = 2;								% number of random fields, 
										% MUST be 2 for this example

% create *.out and open in append mode
[fpath, fname, ~] = fileparts(inpfnm);
% fpath  = pwd;
outext = '.out';
outfnm = fullfile(fpath, [fname,outext]);
if exist(outfnm,'file')==2
	outfid = fopen(outfnm, 'w');        % if exists, clear contents
	fclose(outfid);
end
outfid = fopen(outfnm, 'a');            % append mode

% print current time
t = datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y  HH:mm:ss  Z');
fprintf(outfid, '%s\n', char(t));
fprintf(outfid, 'This is RFGEN (Version 2.0): %s\n', t);
fprintf(outfid, "------------------------------------------------------" + ...
				"------------\n");
fprintf(outfid, '\n');


%% input data
[rf1, rf2, thx, thy, crcoef, dx, dy, nxe, nye, kseed, nsim, ...
					 lplt, ino, irf, ldbg] = getinp(inpfnm);
r = [1, crcoef;
	 crcoef, 1];                    % cross correlation coef matrix
rc = chol(r);                       % Cholesky decomposition

if rf1(3)==3						% mean and sd of bounded dist
	rf1(1) = (rf1(4) + rf1(5)) / 2;
	rf1(2) = 0.46*(rf1(5)-rf1(4))*rf1(7) / sqrt(4*pi^2+(rf1(7))^2);
end
if rf2(3)==3						% mean and sd of bounded dist
	rf2(1) = (rf2(4) + rf2(5)) / 2;
	rf2(2) = 0.46*(rf2(5)-rf2(4))*rf2(7) / sqrt(4*pi^2+(rf2(7))^2);
end

if ~ismember(ino,0:nsim)
	fprintf(outfid, "RFGEN: 'ino' should be 0 to nsim=%d!\n" + ...
					"       Current 'ino'=%d.\n", nsim, ino);
	error("RFGEN: 'ino' should be 0 to nsim=%d!\n" + ...
		  "       Current 'ino'=%d.\n", nsim, ino);
end

if ino~=0
	if ~ismember(irf,[1,2])
		fprintf(outfid, "RFGEN: 'irf' should be 1 or 2!\n" + ...
						"       Current 'irf'=%d.\n", irf);
		error("RFGEN: 'irf' should be 1 or 2!\n" + ...
			  "       Current 'irf'=%d.\n", irf);
	end
else
	irf = [];
end

% write input parameters to *.out
wrtinp(outfid, rf1, rf2, thx, thy, crcoef, dx, dy, nxe, nye, kseed, nsim);


%% pre-processing
% find nrfx, nrfy with corresponding k1, k2, m
[nrfx, nrfy] = chknxe(nxe, nye, mxk);

% if debug on, output some validation results
if ldbg
	% mean and std of generated random fields for validation
	rfmn1  = zeros(nsim, 1);            % mean of rf1 for each realization
	rfsd1 = zeros(nsim, 1);				% std of rf1 for each realization
	rfmn2  = zeros(nsim, 1);			% mean of rf2 for each realization
	rfsd2 = zeros(nsim, 1);				% std of rf2 for each realization

	% correlation structure of generated random fields for validation
	nrfm = min(nrfx, nrfy);
	cor1rf1 = zeros(nrfx, 2);           % corr structure in x-dir, rf1
	cor2rf1 = zeros(nrfy, 2);           % corr structure in y-dir, rf1
	cor0rf1 = zeros(nrfm, 2);           % corr structure in diagonal, rf1
	cor1rf2 = zeros(nrfx, 2);			% corr structure in x-dir, rf2
	cor2rf2 = zeros(nrfy, 2);			% corr structure in y-dir, rf12
	cor0rf2 = zeros(nrfm, 2);			% corr structure in diagonal, rf2

	% correlation coefficient between random fields for validaton
	crcf = zeros(nsim, 1);
end

simrf1 = zeros(nye, nxe, nsim);
simrf2 = zeros(nye, nxe, nsim);
hatqu  = zeros(nsim, 1);


%% now, for each realization...
if ldbg
	fprintf(outfid, "\n");
	fprintf(outfid, "\n");
	fprintf(outfid, "Realizations:\n");
	fprintf(outfid, "-----------------------------------------------" + ...
					"-------------------\n");
end

tic;
for isim = 1 : nsim
	rng(kseed);                         % for reproduction
	lfst = false;
	if isim==1
		lfst = true;
	end

	lpltrf = false;
	if ino==isim
		lpltrf = true;
	end

	% simulate random fields; rf size of nrfy*nrfx, NOT nrfx*nrfy!
	[simrf1i, simrf2i] = sim2bc(rf1, rf2, rc, thx, thy, dx, dy, nrfx, nrfy, ...
								nxe, nye, mxm, mxk, lfst, lpltrf, irf, nrf);
	simrf1(:,:,isim) = simrf1i(1:nye, 1:nxe);	% extract nye*nxe
	simrf2(:,:,isim) = simrf2i(1:nye, 1:nxe);

	kseed = kseed + 1;


	%% if debug on ...
	% check correlation structure
	if ldbg
		lplti = false;                   % suppress plot in `corplt`
		[~, ~, cor1rf1i, cor2rf1i, cor0rf1i] = ...
						 dcor2d(simrf1i, dx, dy, thx, thy, lplti);
		cor1rf1 = cor1rf1 + cor1rf1i;
		cor2rf1 = cor2rf1 + cor2rf1i;
		cor0rf1 = cor0rf1 + cor0rf1i;
		[~, ~, cor1rf2i, cor2rf2i, cor0rf2i] = ...
						 dcor2d(simrf2i, dx, dy, thx, thy, lplti);
		cor1rf2 = cor1rf2 + cor1rf2i;
		cor2rf2 = cor2rf2 + cor2rf2i;
		cor0rf2 = cor0rf2 + cor0rf2i;

		% check correlation coefficient between two random fields
		fprintf(outfid, 'Realization %d: \n', isim);
		if rf1(3)~=0 && rf2(3)~=0
			crcfi = corrcoef(simrf1i(:), simrf2i(:));
			crcf(isim) = crcfi(1, 2);
			fprintf(outfid, '    Corr coef between RFs is: %3.2f\n', ...
					crcf(isim));
		else		% `corrcoef` not work for zero SD (det RF)
			fprintf(outfid, '    Corr coef between RFs is: N/A (one rf is det)\n');
		end
	end

	% check random field mean and std
	if ldbg
		fprintf(outfid, '    Check random field #1: ');
		[rfmn1(isim), rfsd1(isim)] = chkfld(simrf1i, outfid);
		fprintf(outfid, '    Check random field #2: ');
		[rfmn2(isim), rfsd2(isim)] = chkfld(simrf2i, outfid);
		fprintf(outfid, '\n');
	end
end                     % end of isim
tsim = toc;             % total time for all nsim

% statistics
if ldbg
	rfmn1a = mean(rfmn1);				% mean of all `nsim` RF1 MEANs
	rfsd1a = mean(rfsd1);				% mean of all `nsim` RF1 SDs
	rfmn2a = mean(rfmn2);
	rfsd2a = mean(rfsd2);
	crcfa  = mean(crcf);				% mean of all `nsim` crcf
end
qumn = mean(hatqu);
qusd = std(hatqu);

fprintf(outfid, '\n');
fprintf(outfid, '\n');
fprintf(outfid, 'Statistics:\n');
fprintf(outfid, "------------------------------------------------------" + ...
				"------------\n");
fprintf(outfid, 'Total Realizations:	%d\n', nsim);
if ldbg
	fprintf(outfid, '\n');
	fprintf(outfid, 'MEAN of all RF1 MEANs:	%7.5e (actual %7.5e)\n', ...
			rfmn1a, rf1(1));
	fprintf(outfid, 'MEAN of all RF1 SDs:	%7.5e (actual %7.5e)\n', ...
			rfsd1a, rf1(2));
	fprintf(outfid, 'MEAN of all RF2 MEANs:	%7.5e (actual %7.5e)\n', ...
			rfmn2a, rf2(1));
	fprintf(outfid, 'MEAN of all RF2 SDs:	%7.5e (actual %7.5e)\n', ...
			rfsd2a, rf2(2));
	fprintf(outfid, 'MEAN of all crs-corr:	%7.5e (actual %7.5e)\n', ...
			crcfa, crcoef);
end
fprintf(outfid, '\n');
fprintf(outfid, 'MEAN of capacity:	%7.5e\n', qumn);
fprintf(outfid, 'SD   of capacity:	%7.5e\n', qusd);

fprintf(outfid, '\n');
fprintf(outfid, 'Total elapsed time:	%7.5e sec.\n', tsim);
fprintf(1, 'Total elapsed time: %7.5e sec.\n', tsim);

% average correlation structures over all realizations
if ldbg
	cor1rf1 = cor1rf1 / nsim;
	cor2rf1 = cor2rf1 / nsim;
	cor0rf1 = cor0rf1 / nsim;
	cor1rf2 = cor1rf2 / nsim;
	cor2rf2 = cor2rf2 / nsim;
	cor0rf2 = cor0rf2 / nsim;
end


%% if plot on, some plots here
if lplt
	% histograms of simulated RFs
	histplt(simrf1(:), rf1, 'Random field 1');
	histplt(simrf2(:), rf2, 'Random field 2');

	% correlation structures averaged over all realizations
	if ldbg
		corplt(cor1rf1, cor2rf1, cor0rf1, cor1rf2, cor2rf2, cor0rf2, ...
			   dx, dy, thx, thy, nrfx, nrfy, nrf);
	end
end

fclose(outfid);

