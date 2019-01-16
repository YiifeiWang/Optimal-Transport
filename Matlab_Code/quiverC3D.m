function [  ] = quiverC3D(x,y,z,u,v,w,maxNumArrows)
%quiverC3D creates a 3D quiver plot and adds a color coding. The color coding is
%given by the absolut values of the component vectors. Large values result in colors 
%from the upper end of the used colormap. Plotting parameters have to be changed within 
%the function in this version. Values equal to NaN or inf are set to zero.
%In case of complex arrays the absolute value is used.
% 
%   INPUT:
%       x - 2D matrix, x components of initial points
%       y - 2D matrix, y components of initial points
%       z - 2D matrix, z components of initial points
%       u - 2D matrix, x components of arrows
%       v - 2D matrix, y components of arrows
%       w - 2D matrix, z components of arrows
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
%       z = sin(x.^2 + y.^2);
%       u = exp(-0.2*(x-5).^2 - 0.2*(y-5).^2);
%       v = -u;
%       w = u.*v;
%       quiverC3D(x,y,z,u,v,w,500);
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER3, LINESPEC, COLORMAP.
% 

%% prearrangements

narginchk(6,7);
n_in = nargin;

if ~isequal(size(x),size(y),size(z),size(u),size(v),size(w))
    error('X,Y,Z,U,V,W have to be matrices of the same size.');
end


%% assignments

% maximum number of arrows if necessary
if n_in == 6
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
    z = z(1:N:end,1:N:end);
    u = u(1:N:end,1:N:end);
    v = v(1:N:end,1:N:end);
    w = w(1:N:end,1:N:end);
end

s = size(u);

%% taking care of possible issues

x = issues(x);
y = issues(y);
z = issues(z);
u = issues(u);
v = issues(v);
w = issues(w);

if ~isequal(u,zeros(s)) 
    u = u./norm(u);
end
if ~isequal(v,zeros(s))
    v = v./norm(v);
end
if ~isequal(w,zeros(s))
    w = w./norm(w);
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
            quiver3(x(n,m),y(n,m),z(n,m),u(n,m),v(n,m),w(n,m),'Color',C(Ic(n,m),:),'LineWidth',lw,'maxheadsize',hs);
        end
    end
end

hold off;

end


