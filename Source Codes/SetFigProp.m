function axes_h = SetFigProp(normal_large, str_xlabel, str_ylabel)
%%
%	SETFIGHPROP: only used to setup some figure properties. It can be used
%				 for all generated figures in a consistent way. This is not
%				 necessary for random generation.
%


% Create figure
h = figure;

% Create axes
axes_h = axes('Parent', h);
hold(axes_h, 'on');

% Create xlabel and ylabel
xlabel(str_xlabel, 'FontName', 'Times New Roman', 'interpreter', 'latex');
ylabel(str_ylabel, 'FontName', 'Times New Roman', 'interpreter', 'latex');

% grid and box
box(axes_h, 'on');
grid(axes_h, 'on');

% fullscreen or normal size, font size
if (normal_large == 'n')
    set(h, 'WindowStyle', 'normal', 'Position', get(0,'defaultfigureposition'));
    set(axes_h, 'FontSize', 12);
    set(axes_h, 'BoxStyle', 'back', 'LineWidth', 0.5);
elseif (normal_large == 's')
    set(h, 'WindowStyle', 'normal', 'Position', get(0,'defaultfigureposition'));
    set(axes_h, 'FontSize', 14);
    set(axes_h, 'BoxStyle', 'back', 'LineWidth', 0.5);
elseif (normal_large == 'l')
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    set(axes_h, 'FontSize', 22);
    set(axes_h, 'BoxStyle', 'back', 'LineWidth', 0.5);
else
    error('Figure size error!');
end

% Set the remaining axes properties
set(axes_h, ...
    'FontName', 'Times New Roman', ...
    'GridAlpha', 0.15, ...
    'GridLineStyle', '--');

