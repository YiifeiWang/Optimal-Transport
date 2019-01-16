function [out] = LPER_admm(model)

m = model.m;
n = model.n;
obj = reshape(model.obj, m, n);
cst = model.cst;
eps = model.eps;
t = model.t;
mu = cst(1:m);
vu = cst(m+1:m+n);
if ~isfield(model, 'gamma')
    gamma = zeros(m, 1);
else
    gamma = model.gamma;
end
if ~isfield(model, 'lmd')
    lmd = zeros(n, 1);
else
    lmd = model.lmd;
end
if ~isfield(model, 'omg')
    omg = zeros(m, n);
else
    omg = model.omg;
end
if ~isfield(model, 'X')
    X = zeros(m, n);
else
    X = model.X;
end
if ~isfield(model, 'XX')
    XX = zeros(m, n);
else
    XX = model.XX;
end
iter = 1;
if ~isfield(model, 'alpha')
    model.alpha = 1;
end
tic;
while iter <= model.iter
    R = (omg + (gamma + lmd') - obj)/t + (mu + vu') + XX;
    X = (R - (sum(R, 1)/(m+1)+sum(R, 2)/(n+1)) + (1/(n+1)+1/(m+1))/(m+n+1) * sum(sum(R)));
    XX = XX+1e-16;
    tmp = XX - (eps*log(XX)+omg+t.*XX-t*X)./(eps./XX+t); 
    %tmp = X-omg/t;
    XX = max(tmp, 0);
    gamma = gamma + model.alpha * t * (mu-sum(X, 2));
    lmd = lmd + model.alpha * t * (vu-sum(X, 1)');
    omg = omg + model.alpha * t * (XX-X);
    diff = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
    if mod(iter, 1000) == 0
        norm_diff = norm(X-XX, 1);
        fprintf("Admm - Iter: %d: objval: %.9f vltcst %.9f D-P diff: %.9f\n", iter, sum(sum(obj.*X)), diff, norm_diff);
    end
    iter = iter + 1;
    if (diff < model.tol)
        break;
    end
end
time = toc;
out.m = m;
out.n = n;
out.obj = obj;
out.t = t;
out.cst = cst;
out.eps = eps;
out.gamma = gamma;
out.lmd = lmd;
out.omg = omg;
out.X = X;
out.XX = XX;
out.iter = iter - 1;
out.vltcst = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
out.objval = sum(sum(obj.*X));
out.entval = out.objval + eps * sum(sum(X.*(log(X)-1)));
out.time = time;

end
