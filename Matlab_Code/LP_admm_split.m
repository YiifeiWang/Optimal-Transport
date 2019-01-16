function [out] = LP_admm_split(model)

m = model.m;
n = model.n;
obj = reshape(model.obj, m, n);
if ~isfield(model, 't')
    model.t = 4 * (m+n) * mean(mean(obj)); 
end
t = model.t;
cst = model.cst;
tol = model.tol;
mu = cst(1:m);
vu = cst(m+1:m+n);
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
    X = XX - (0.5 * obj + omg)/t;
    X = LP_l1_proj(X', mu')';
    XX = X + (omg - 0.5 * obj)/t;
    XX = LP_l1_proj(XX, vu');
    omg = omg + model.alpha * t * (X-XX);
    diff = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
    norm_diff = norm(X-XX, 1);
    if mod(iter, 1000) == 0
        fprintf("Split - Iter: %d objval: %.9f vltcst: %.9f D-P gap: %.9f\n", iter, sum(sum(obj.*(X+XX)/2)), diff, norm_diff);
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

