function [out] = LP_admm_split_(model)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

m = model.m;
n = model.n;
obj = reshape(model.obj, m, n);
if ~isfield(model, 't')
    model.t = 100 * (m+n) * mean(mean(obj)); 
end
if ~isfield(model, 'tol_sol')
    model.tol_sol = 1e-6;
end
t = model.t;
cst = model.cst;
tol = model.tol;
tol_sol = model.tol_sol;
mu = cst(1:m);
vu = cst(m+1:m+n);
gamma = zeros(m, 1);
lmd = zeros(n, 1);
omg = zeros(m, n);
X = zeros(m, n);
XX = zeros(m, n);
iter = 1;
if ~isfield(model, 'iter')
    model.iter = 20000;
end
if ~isfield(model, 'alpha')
    model.alpha = 1;
end
tic;
while (iter <= model.iter)
    RX = XX + mu - (0.5 * obj + gamma + omg)/t;
    X = max(RX - sum(RX, 2)/(n+1), 0);
    RXX = X + vu' - (0.5 * obj + lmd' - omg)/t;
    XX = max(RXX - sum(RXX, 1)/(m+1), 0);
    gamma = gamma + t * (sum(X, 2) - mu);
    lmd = lmd + t * (sum(X, 1)' - vu);
    omg = omg + model.alpha * t * (X-XX);
    diff = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
    norm_diff = norm(X-XX, 1);
    if mod(iter, 10) == 0
        fprintf("%d: %.9f %.9f %.9f\n", iter, sum(sum(obj.*X)), diff, norm_diff);
    end
    %if mod(iter, 250) == 0
    %    t = max(t / 2, (m+n) * mean(mean(obj)));
    %end
    iter = iter + 1;
    if (diff < tol && iter >= 2000)
        break; 
    end
end
time = toc;
out.m = m;
out.n = n;
out.t = t;
out.obj = obj;
out.cst = cst;
out.iter = iter - 1;
out.X = (X+XX)/2;
out.objval = sum(sum(obj.*out.X));
out.vltcst = norm([sum(out.X, 1)-vu', sum(out.X, 2)'-mu'], 1);
out.time = time;

end

