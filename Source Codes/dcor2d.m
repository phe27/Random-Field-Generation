function [mn, sd, cor1, cor2, cor] = dcor2d(rf, dx, dy, thx, thy, lplt)
%%
% DCOR2D: Estiamte the correlation function of a 2-dimension random field.
%
%   Directional (x-, y-, and diagonal direction) correlation unctions are
%   evaluated. Sample mean and standard deviation are also calculated. See
%   Fenton and Griffiths (2008) Chapter 5.3.6 for more details.
%
% INPUTS:
%
%   rf:     random field for evaluation
%
%   dx:     element size in x-direction
%
%   dy:     element size in y-direction
%
%   thx:    correlation length in x-direction
%
%   thy:    correlation length in y-direction
%
%   lplt:   plot figure?
%
% OUTPUTS:
%
%   mn:     sample mean
%
%   sd:     sample standard deviation
%
%   cor1:   sample correlation function in x-direction
%
%   cor2:   sample correlation function in y-direction
%
%   cor:    sample isotropic correlation function
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
%  1.1  changed diagonal correlation structure (Mar 25/2024)
%


%% mean, variance, standard deviation
mn = mean(rf(:));
vr = var(rf(:));
sd = std(rf(:));


%% correlation functions
[n2, n1] = size(rf);
% xl  = n1 * dx;
% yl  = n2 * dy;
% tau = sqrt(xl^2+yl^2) / max(n1,n2);

% x-direction
tau1 = dx;
cor1 = zeros(n1, 2);

for jj = 0 : n1-1
    % isum = 0;
    for kk = 1 : n2
        for ii = 1 : n1-jj
            cor1tmp = (rf(kk,ii)-mn) * (rf(kk,ii+jj)-mn);
            cor1(jj+1,2) = cor1(jj+1,2) + cor1tmp;
            % isum = isum + 1;
        end
    end
    % fprintf(1, 'DCOR2D (x-dir): isum=%d; isum2=%d\n', isum, (n2*(n1-jj)-1));
    coef = 1 / vr / (n2*(n1-jj)-1);
    cor1(jj+1,2) = cor1(jj+1,2) * coef;
    cor1(jj+1,1) = jj * tau1;
end

% y-direction
tau2 = dy;
cor2 = zeros(n2, 2);

for jj = 0 : n2-1
    % isum = 0;
    for kk = 1 : n1
        for ii = 1 : n2-jj
            cor2tmp = (rf(ii,kk)-mn) * (rf(ii+jj,kk)-mn);
            cor2(jj+1,2) = cor2(jj+1,2) + cor2tmp;
            % isum = isum + 1;
        end
    end
    % fprintf(1, 'DCOR2D (y-dir): isum=%d; isum2=%d\n', isum, (n1*(n2-jj)-1));
    coef = 1 / vr / (n1*(n2-jj)-1);
    cor2(jj+1,2) = cor2(jj+1,2) * coef;
    cor2(jj+1,1) = jj * tau2;
end

% diagonal
tau = sqrt(dx^2 + dy^2);
nm  = min(n1, n2);
cor = zeros(nm, 2);

for jj = 0 : nm-1
    isum = 0;                               % counter of total data pairs
    for kk = 1 : n2-jj
        for ii = 1 : n1-jj
            cortmp = (rf(kk,ii)-mn) * (rf(kk+jj,ii+jj)-mn);
            cor(jj+1,2) = cor(jj+1,2) + cortmp;
            isum = isum + 1;
        end
    end
    % fprintf(1, 'DCOR2D (diag-dir): isum=%d\n', isum);
    coef = 1 / vr / isum;
    cor(jj+1,2) = cor(jj+1,2) * coef;
    cor(jj+1,1) = jj * tau;
end


%% plot and compare with theoretical curve
if lplt
    corplt(cor1, cor2, cor, [], [], [], dx, dy, thx, thy, n1, n2, 1);
end
