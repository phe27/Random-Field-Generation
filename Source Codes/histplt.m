function histplt(hdata, stat, xname)
%%
% HISTPLT: Plot a histogram for sample data. Theoretical curve only takes
%          deterministic (0), normal (1), lognormal (2), and bounded (3)
%          distributions so far.
%
% INPUTS:
%
%   hdata:  data for histogram
%
%    stat:  real vector of length at least 7 containing the parameters of
%           theoretical curves. Notably,
%              stat(1) = mean,
%              stat(2) = standard deviation,
%              stat(3) = distribution type;
%                        = 0.0 if deterministic (at mean value)
%                        = 1.0 if normally distributed
%                        = 2.0 if lognormally distributed
%                        = 3.0 if bounded distributed
%              stat(4) = lower bound (if bounded)
%              stat(5) = upper bound (if bounded)
%              stat(6) = m parameter (if bounded)
%              stat(7) = s parameter (if bounded)
%           If this is bounded, stat(1) and stat(2) are ignored and the
%           values stat(4) through stat(7) completely describe the dist.
%
%   xname:  x axis name
%
% OUTPUTS:
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
% Created:      Mar 22, 2024
% Version:      v1.1
% Revision history:
%  1.1  added 'bounded' distribution. (Mar 16/2024)
%


%%
dist = stat(3);
% sample histogram
nbin = ceil( 1 + log2(length(hdata)) ) + 3;
SetFigProp('s', xname, 'PDF');
histogram(hdata, nbin, Normalization="pdf");

spmn = mean(hdata);                         % sample mean
spsd = std(hdata);                          % sample std

% theoretical dist curve
mn  = stat(1);                              % theoretical mean
sd  = stat(2);                              % theoretical std
cov = sd / mn;                              % theoretical COV
if dist==0                                          % deterministic dist
    sd = 0;
elseif dist==1                                      % normal dist
    pd = makedist("Normal", mu=mn, sigma=sd);
elseif dist==2                                      % lognormal dist
    sdln = sqrt(log(1 + cov^2));            % sd of log values
    mnln = log(mn) - 0.5*sdln^2;            % mean of log values
    pd = makedist("Lognormal", mu=mnln, sigma=sdln);
elseif dist==3                                      % bounded dist
    a = stat(4);
    b = stat(5);
    m = stat(6);
    s = stat(7);
    [fnpdf, mn, sd] = bddist(a, b, m, s);
else
    error('HISTPLT: Unknown distribution type!')
end

if dist==0
    warning('HISTPLT: Deterministic dist; NO theoretical curve!');
elseif dist==1 || dist==2
    xmin = max([mn-4*sd, 0]);
    xmax = mn + 5*sd;
    x = linspace(xmin, xmax, 200);
    y = pdf(pd, x);
    hold on;
    plot(x, y, 'r-', LineWidth=1);
    hold off;
elseif dist==3
    xmin = a;
    xmax = b;
    x = linspace(xmin, xmax, 200);
    y = fnpdf(x);
    hold on;
    plot(x, y, 'r-', LineWidth=1);
    hold off;
end

if dist~=0
    txtstr = sprintf(['Sample mean: %.2f\nSample std: %.2f\nTheoretical mean: ' ...
                      '%.2f\nTheoretical std: %.2f'], spmn, spsd, mn, sd);
    annotation('textbox', [0.571 0.735 0.325 0.180], ...
               'String', txtstr, ...
               'Color',[1 0 0],...
               'Interpreter','latex',...
               'FontSize',11,...
               'FitBoxToText','off',...
               'EdgeColor',[1 0 0]);
end


