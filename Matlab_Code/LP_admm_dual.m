function [out] = LP_admm_dual(model)
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
t = model.t;
obj = reshape(model.obj, m, n);
cst = model.cst;
tol = model.tol;
mu = cst(1:m);
vu = cst(m+1:m+n);
if ~isfield(model, 'alpha')
    model.alpha = 1;
end
z = zeros(m, n);
omg = zeros(m, n);
lmd = zeros(n, 1);
iter = 1;
tic;
while iter <= model.iter
    a = sum(z, 2) - mu/t - sum(omg, 2)/t - sum(obj, 2);
    b = sum(z, 1)' - vu/t - sum(omg, 1)'/t - sum(obj, 1)';
    b = b(1:end-1);
    pi_sum = sum(a) - sum(b);
    lmd_sum = (n*sum(b) - (n-1)* sum(a))/m;
    pi = (a-lmd_sum)/n;
    lmd(1:end-1) = (b-pi_sum)/m;
    z = max(obj + pi + lmd' + omg/t, 0);
    omg = omg + model.alpha * t * (obj + pi + lmd'-z);
    iter = iter + 1;
	diff = norm([sum(omg, 1)+vu', sum(omg, 2)'+mu'], 1);
    if mod(iter, 1000) == 0
        fprintf("Dual - Iter: %d objval: %.9f vltcst: %.9f\n", iter, -pi'*mu - lmd'*vu, diff);
    end
    if (iter > 1000 && diff < tol)
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
out.X = -omg;
out.objval = sum(sum(obj.*out.X));
out.vltcst =  norm([sum(out.X, 1)-vu', sum(out.X, 2)'-mu'], 1);
out.time = time;

end