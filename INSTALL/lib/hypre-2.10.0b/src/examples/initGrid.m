%======================
%  initGrid
%======================
%  initialize grid
%
% INPUT:
%       N      = # gpd 
%       ax, bx = lower and upper bound in x-dir.
%       ay, by = lower and upper bound in y-dir.
%       CCgrid = 1 for cell centred grid, 0 for Vertex-Centred
%
% OUTPUT
%   X,Y  = meshgrid
%     h  = spacing


function [X,Y,h] = initGrid(N,ax,bx,ay,by,CCgrid);


switch CCgrid
    case {0}  % Vertex Centred
        nx = N;
        ny = nx;
        h = (bx-ax)/(nx-1);
        ii = 1:nx; jj = 1:ny;
        x = ax + (ii-1)*h;
        y = ay + (jj-1)*h;
        [X,Y] = ndgrid(x,y);	% like meshgrid but x index is first, y is second;
        
    case {1}  % Cell-Centred
        l = floor(log(N)/log(2));  % to use max. levels in multigrid method
        nx = 2^l + 2;
        ny = nx;
        h = (bx-ax)/(2^l);
        ii = 0:nx-1;
        jj = 0:nx-1;
        x=ax + (ii - 0.5)*h;
        y=ay + (jj - 0.5)*h;
        [X,Y] = ndgrid(x,y);
        
    otherwise
        'grid specification is wrong, set CCGrid to 1 or 0'
end;