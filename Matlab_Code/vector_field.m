function [] = vector_field(opts)
C = colormap(gray);
hold on;
source = opts.source;
target = opts.target;
s = opts.m;
solx = opts.x;
solx = reshape(round(solx/max(max(solx))*24), s, s);
for n = 1:s
    for m = 1:s
        if solx(n, m) > 0
            quiver(source(m,1),source(m,2),target(n,1)-source(m,1),target(n,2)-source(m,2),1,'Color',C(64-solx(n,m),:),'LineWidth',0.75,'maxheadsize',0.1);
        end
    end
end
scatter(source(:, 1), source(:, 2), 'r', '.');
scatter(target(:, 1), target(:, 2), 'b', '.');
alpha(0.5);
axis equal;
end

