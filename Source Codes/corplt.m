function corplt(cor1rf1, cor2rf1, cor0rf1, cor1rf2, cor2rf2, cor0rf2, ...
                dx, dy, thx, thy, nrfx, nrfy, nrf)
%%
% CORPLT: Plot estimated correlation structures of random fields.
%
%   plot estimated correlation structures of simulated random fields in x-,
%   y-, and diagonal directions, and compare with theoretical curves
%   (Markov correlation function).
%
% INPUTS:
%
%  cor1rf1: sample correlation function of random field #1 in x-dir
%
%  cor2rf1: sample correlation function of random field #1 in y-dir
%
%  cor0rf1: sample correlation function of random field #1 in diagonal dir
%
%  cor1rf2: sample correlation function of random field #2 in x-dir.
%           if nrf==1, it is [].
%
%  cor2rf2: sample correlation function of random field #2 in y-dir.
%           if nrf==1, it is [].
%
%  cor0rf2: sample correlation function of random field #2 in diagonal dir.
%           if nrf==1, it is [].
%
%   dx:     element size in x-dir
%
%   dy:     element size in y-dir
%
%   thx:    correlation length in x-dir
%
%   thy:    correlation length in y-dir
%
%   nrf:    how many random fields? 1 or 2
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
%  1.1  changed line styles (Mar 25/2024)
%


% correlation structures averaged over nsim realizations
tau1 = dx;
tau2 = dy;
tau  = sqrt(dx^2 + dy^2);
n1   = nrfx;
n2   = nrfy;

if thx==thy                             % isotropic case
    % random field #1
    SetFigProp('s', 'Lag, $\tau$ [m]', '$\rho(\tau)$');
    fun = @(x) exp(-2 * x / thx);
    h4 = fplot(fun, [0, max(n1*tau1,n2*tau2)], 'k--', LineWidth=1);

    h1 = plot(cor1rf1(:,1), cor1rf1(:,2), 'b-', LineWidth=1);
    h2 = plot(cor2rf1(:,1), cor2rf1(:,2), 'r-', LineWidth=1);
    h3 = plot(cor0rf1(:,1), cor0rf1(:,2), 'k-', LineWidth=1);
    legend([h1,h2,h3,h4], ...
           {'Simulation: $x$ direction', ...
            'Simulation: $y$ direction', ...
            'Simulation: diagonal', ...
            'Theoretical'}, ...
           'Location', 'northeast', ...
           'Interpreter', 'latex');
    title('Random field #1');
    
    % random field #2
    if nrf==2
        SetFigProp('s', 'Lag, $\tau$ [m]', '$\rho(\tau)$');
        fun = @(x) exp(-2 * x / thx);
        h4 = fplot(fun, [0, max(n1*tau1,n2*tau2)], ...
                   'k--', LineWidth=1);

        h1 = plot(cor1rf2(:,1), cor1rf2(:,2), 'b-', LineWidth=1);
        h2 = plot(cor2rf2(:,1), cor2rf2(:,2), 'r-', LineWidth=1);
        h3 = plot(cor0rf2(:,1), cor0rf2(:,2), 'k-', LineWidth=1);
        legend([h1,h2,h3,h4], ...
               {'Simulation: $x$ direction', ...
                'Simulation: $y$ direction', ...
                'Simulation: diagonal', ...
                'Theoretical'}, ...
               'Location', 'northeast', ...
               'Interpreter', 'latex');
        title('Random field #2');
    end
    
else                                    % cross-anisotropic case
    % random field #1
    SetFigProp('s', 'Lag, $\tau$ [m]', '$\rho(\tau)$');
    funx = @(x) exp(-2 * x / thx);
    funy = @(y) exp(-2 * y / thy);
    a = 1 / sqrt(1 + (dy/dx)^2);
    b = a * (dy/dx);
    fun  = @(z) exp( -2 * sqrt((a*z/thx).^2 + (b*z/thy).^2) );
    h4 = fplot(funx, [0, n1*tau1], 'b--', LineWidth=1);
    h5 = fplot(funy, [0, n2*tau2], 'r--', LineWidth=1);
    nmin = min(n1, n2);
    h6 = fplot(fun, [0, nmin*tau], 'k--', LineWidth=1);

    h1 = plot(cor1rf1(:,1), cor1rf1(:,2), 'b-', LineWidth=1);
    h2 = plot(cor2rf1(:,1), cor2rf1(:,2), 'r-', LineWidth=1);
    h3 = plot(cor0rf1(:,1), cor0rf1(:,2), 'k-', LineWidth=1);
    legend([h1,h2,h3,h4,h5,h6], ...
           {'Simulation: $x$ direction', ...
            'Simulation: $y$ direction', ...
            'Simulation: diagonal', ...
            'Theoretical: $x$ direction', ...
            'Theoretical: $y$ direction', ...
            'Theoretical: diagonal'}, ...
           'Location', 'northeast', ...
           'Interpreter', 'latex');
    title('Random field #1');
    
    % random field #2
    if nrf==2
        SetFigProp('s', 'Lag, $\tau$ [m]', '$\rho(\tau)$');
        funx = @(x) exp(-2 * x / thx);
        funy = @(y) exp(-2 * y / thy);
        a = 1 / sqrt(1 + (dy/dx)^2);
        b = a * (dy/dx);
        fun  = @(z) exp( -2 * sqrt((a*z/thx).^2 + (b*z/thy).^2) );
        h4 = fplot(funx, [0, n1*tau1], 'b--', LineWidth=1);
        h5 = fplot(funy, [0, n2*tau2], 'r--', LineWidth=1);
        nmin = min(n1, n2);
        h6 = fplot(fun, [0, nmin*tau], 'k--', LineWidth=1);

        h1 = plot(cor1rf2(:,1), cor1rf2(:,2), 'b-', LineWidth=1);
        h2 = plot(cor2rf2(:,1), cor2rf2(:,2), 'r-', LineWidth=1);
        h3 = plot(cor0rf2(:,1), cor0rf2(:,2), 'k-', LineWidth=1);
        legend([h1,h2,h3,h4,h5,h6], ...
               {'Simulation: $x$ direction', ...
                'Simulation: $y$ direction', ...
                'Simulation: diagonal', ...
                'Theoretical: $x$ direction', ...
                'Theoretical: $y$ direction', ...
                'Theoretical: diagonal'}, ...
               'Location', 'northeast', ...
               'Interpreter', 'latex');
        title('Random field #2');
    end

end
