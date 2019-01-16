function [  ] = quiverC2D(x,y,u,v,maxNumArrows)
%quiverC2D creates a 2D quiver plot and adds a color coding. The color coding is
%given by the absolut values of the component vectors. Large values result in colors 
%from the upper end of the used colormap. Plotting parameters have to be changed within 
%the function in this version. Values equal to NaN or inf are set to zero.
%In case of complex arrays the absolute value is used.
% 
%   INPUT:
%       x - 2D matrix, x components of initial points
%       y - 2D matrix, y components of initial points
%       u - 2D matrix, x components of arrows
%       v - 2D matrix, y components of arrows
%       maxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If maxNumArrows is not given its value is set to 1000.
% 
%   OUTPUT:
%       none
% 
%   WARNING!: Using large datasets in combination with choosing maxNumArrows 
%       very large might result in this script running forever.
% 
% --------------------------------------------------------------------------------------
% 
%   EXAMPLE: 
%       [x,y] = meshgrid(linspace(0,10,100),linspace(0,10,100));
%       u = exp(-0.2*(x-5).^2 - 0.2*(y-5).^2);
%       v = -u;
%       quiverC2D(x,y,u,v,1000);
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER, LINESPEC, COLORMAP.
% 
% 

%% prearrangements

narginchk(4,5);
n_in = nargin;

if ~isequal(size(x),size(y),size(u),size(v))
    error('X,Y,U,V have to be matrices of the same size.');
end


%% assignments

% maximum number of arrows if necessary
if n_in == 4
    maxNumArrows = 1000;
end
% Line width
lw = 2;
% Maximum of arrow head size
hs = 2;
% Colormap
colormap parula;

%% initialization
if numel(u) > maxNumArrows
    N = ceil(sqrt(numel(u)/maxNumArrows));
    
    x = x(1:N:end,1:N:end);
    y = y(1:N:end,1:N:end);
    u = u(1:N:end,1:N:end);
    v = v(1:N:end,1:N:end);
end

s = size(u);

%% taking care of possible issues

x = issues(x);
y = issues(y);
u = issues(u);
v = issues(v);

if ~isequal(u,zeros(s)) 
    u = u./norm(u);
end
if ~isequal(v,zeros(s))
    v = v./norm(v);
end



%% colormap definition
I = sqrt(u.^2 + v.^2);
Ic = round( I/max(max(I))*64);
Ic( Ic == 0) = 1;
C = colormap;


%% plotting
hold on;

for n = 1:s(1)
    for m = 1:s(2)
        if u(n,m) > 0 || v(n,m) > 0
            quiver(x(n,m),y(n,m),u(n,m),v(n,m),'Color',C(Ic(n,m),:),'LineWidth',lw,'maxheadsize',hs);
        end
    end
end

hold off;

end


