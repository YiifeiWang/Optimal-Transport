function [out] = LP_admm_primal(model)
%-------------------------------------------------------------------------
% This program implements ADMM for LP primal problem
% 
% Input:
%     model --- the LP model structure with fields:
%				--- general ---
%               m, n   dimension of rows and cols
%				obj    matrix C
%               cst    constraints
%               --- admm ---
%               t      the multiplier of the second order term
%               alpha  shrinkage parameters when updating multipliers
%               
%
% Output:
%       out --- the output structure with fields:
%				--- general ---
%				X     optimal solution
%				objval objective value
%               vltcst violation of constraints
%				time   time elapsed
%				--- admm ---
%               m, n   dimensions of rows and cols
%               t      the multiplier of the second order term
%               obj    matrix C
%               cst    constraints
%				iter   number of iteration
%               
%
% Author: Yifei Wang & Feng Zhu
% Version 1.1 .... 2018/12
%%-------------------------------------------------------------------------

m = model.m;
n = model.n;
obj = reshape(model.obj, m, n);
tol = model.tol;
t = model.t;
cst = model.cst;
mu = cst(1:m);
vu = cst(m+1:m+n);
pi = zeros(m, 1);
lmd = zeros(n, 1);
omg = zeros(m, n);
X = zeros(m, n);
XX = zeros(m, n);
iter = 1;
if ~isfield(model, 'alpha')
    model.alpha = 1;
end
tic;
while iter <= model.iter
    R = (omg + (pi + lmd') - obj)/t + (mu + vu') + XX;
    X = (R - (sum(R, 1)/(m+1)+sum(R, 2)/(n+1)) + (1/(n+1)+1/(m+1))/(m+n+1) * sum(sum(R)));
    XX = max(X-omg/t, 0);
    pi = pi + model.alpha * t * (mu-sum(X, 2));
    lmd = lmd + model.alpha * t * (vu-sum(X, 1)');
    omg = omg + model.alpha * t * (XX-X);
    diff = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
    if mod(iter, 1000) == 0
        norm_diff = norm(X-XX, 1);
        fprintf("Primal - Iter: %d objval: %.9f vltcst: %.9f D-P gap: %.9f\n", iter, sum(sum(obj.*X)), diff, norm_diff);
    end
    iter = iter + 1;
    if (diff < tol && iter > 1000)
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
out.X = X;
out.objval = sum(sum(obj.*out.X));
out.vltcst = norm([sum(out.X, 1)-vu', sum(out.X, 2)'-mu'], 1);
out.time = time;

end

