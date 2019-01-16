function [out] = LP_augmented_lagrangian(model)

m = model.m;
n = model.n;
obj = model.obj;
t = model.t;
cst = model.cst;
mu = cst(1:m);
vu = cst(m+1:m+n);
omg = model.omg;
if ~isfield(model, 'X')
    X = zeros(m, n); 
else
    X = model.X;
end
iter = 0;
X_pre = X;
while iter < model.iter
    %while norm(diff) > model.eps
    %diff_pre = diff;
    %diff = c - (pi + lmd') + t * (sum(X, 2) + sum(X, 1) - (mu + vu')) + t * max(-sign(X), 0);%min(0, X+omg/t);
    X = X + (iter-1)/(iter+2) * (X - X_pre);
    X_pre = X;
    R = -obj + t * (mu + vu') + omg * X;
    Z = (R - t*(sum(R, 2)/(n*t+omg) + sum(R, 1)/(m*t+omg)) + t^2*(1/(n*t+omg)+1/(m*t+omg))/(m*t+n*t+omg) * sum(sum(R)))/omg;
    Y = max(2*Z-X, 0);% + min(2*Z-X+t*model.l1_pen, 0);
    X = X + Y - Z;
        %X = X - model.step * diff;
        %d = diff - diff_pre;
        %d = sum(sum(diff .* diff))/sum(sum(diff .* d));
        %model.step = max(min(model.step*d, model.max), model.min);
        %fprintf("%f\n", norm(diff));
    %end
    %pi = pi + t * (mu-sum(X, 2));
    %lmd = lmd + t * (vu-sum(X, 1)');
    iter = iter + 1;
    gpi = norm([sum(X, 1)-vu', sum(X, 2)'-mu']);
    fprintf("%f %f\n", sum(sum(obj.*X)), gpi);
end
out.m = m;
out.n = n;
out.t = t;
out.obj = obj;
out.cst = cst;
out.omg = omg;
out.X = X;
out.iter = model.iter;
out.gpi = gpi;
out.objval = sum(sum(obj.*X));

end

