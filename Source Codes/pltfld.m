function im = pltfld(simrf, nxe, nye, dx, dy, tl)
%%
% PLTFLD: Plot a random field graphically.
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
%   In matlab, the built-in function 'imagsc' plot it following the order
%   of the input matrix, e.g., (1,1) of the matrix will be placed at the
%   top left corner, (1,2) is located at one element to right, (2,1), is
%   located at one element to downward. So, the original simrf(nxe,nye)
%   need to manipulated to ensure the physical elements being properly
%   placed. That is, in the above example, (1,4) in the original matrix
%   simrf should be transformed to be (1,1) of the input matrix of 'imagsc'.
%
% INPUTS:
%
%   simrf:  simulated random field of size nye*nxe
%
%   nxe:    number of element in x-dir
%
%   nye:    number of element in y-dir
%
%   dx:     element size in x-dir
%
%   dy:     element size in y-dir
%
%   tl:		figure title
%
% OUTPUTS:
%
%   im:     image object
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
% Version:      v1.2
% Revision history:
%  1.1  Remove matrix manipulation. Now it is done in `sim2bc`.
%       (Mar 22/2024)
%  1.2	Added title and changed font. (Jun 20/2024)
%


% for mesh: x and y locations of the center of the first element (top left)
% and the center of the last element (bottom right)
x = [dx/2, nxe*dx-dx/2];
y = [dy/2, nye*dy-dy/2];

% manipulate the original matrix
% simrf = flipud(simrf');

% plot grayscale image with specified element sizes
figure;
im = imagesc(x, y, simrf);

axis equal tight;                       % aspect ratio preserved
colormap gray;                          % grayscale
% colormap turbo;
colorbar;

title(tl,'FontSize', 12,'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex', ...
		'XAxisLocation', 'top', ...
		'XMinorTick', 'on', ...
		'YMinorTick','on');
