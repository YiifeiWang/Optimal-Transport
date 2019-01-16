function [out] = LP_smpf(model)
%-------------------------------------------------------------------------
% This program implements splitting method with penalty functions
% for LP primal problem
% 
% Input:
%     model --- the LP model structure with fields:
%				--- general ---
%               m, n   dimension of rows and cols
%				obj    matrix C
%               cst    constraints
%               --- admm ---
%               t      the multiplier of the second order term
%               omg    proximal projection term
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
%               t      the coefficient of penalty function
%               omg    the coefficient of proximal projection term
%               obj    matrix C
%               cst    constraints
%				iter   maximum number of iteration
%               
%
% Author: Yifei Wang & Feng Zhu
% Version 1.1 .... 2018/12
%%-------------------------------------------------------------------------

m = model.m;
n = model.n;
obj = model.obj;
cst = model.cst;
t = model.t;
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
tic;
while iter < model.iter
    X = X + (iter-1)/(iter+2) * (X - X_pre);
    X_pre = X;
    R = -obj + t * (mu + vu') + omg * X;
    Z = (R - t*(sum(R, 2)/(n*t+omg) + sum(R, 1)/(m*t+omg)) + t^2*(1/(n*t+omg)+1/(m*t+omg))/(m*t+n*t+omg) * sum(sum(R)))/omg;
    Y = max(2*Z-X, 0);
    X = X + Y - Z;
    iter = iter + 1;
    fprintf("%f %f\n", sum(sum(obj.*X)), norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1));
end
time = toc;

out.m = m;
out.n = n;
out.t = t;
out.omg = omg;
out.obj = obj;
out.cst = cst;
out.iter = model.iter;
out.X = X;
out.vltcst = norm([sum(X, 1)-vu', sum(X, 2)'-mu'], 1);
out.objval = sum(sum(obj.*X));
out.time = time;

end

