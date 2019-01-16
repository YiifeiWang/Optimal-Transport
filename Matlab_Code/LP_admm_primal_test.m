function [out, out1] = LP_admm_primal_test(model, opts)
%-------------------------------------------------------------------------
% This program tests ADMM for LP primal problem
% 
% Input:
%     model --- the LP model structure with fields:
%				--- general ---
%               m, n   dimension of rows and cols
%				obj    matrix C
%               cst    constraints
%      opts --- the option structures with fields
%               t1     (admm) the multiplier of the second order term
%               iter1  (admm) the number of iterations
%               omg    (smpf) the coefficient of proximal projection term
%               t2     (smpf) the coefficient of penalty function
%               iter2  (smpf) the number of iterations
%               
%
% Output:
%       out --- the output structure with fields:
%				--- general ---
%				pi     optimal solution
%				objval objective value
%               vltcst violation of constraints
%				time   time elapsed
%				--- admm ---
%               pi_     optimal solution (only admm)
%				objval_ objective value (only admm)
%               vltcst_ violation of constraints (only admm)
%				time_   time elapsed (only admm)
%				iter    number of total iteration
%
%
% Author: Yifei Wang & Feng Zhu
% Version 1.1 .... 2018/12
%%-------------------------------------------------------------------------

if ~isfield(opts, 't1')  
    opts.t1 = 5 * (model.m+model.n)*mean(mean(model.obj)); end
if ~isfield(opts, 'iter1')
    opts.iter1 = 20000;  end
% if ~isfield(opts, 'omg')  
%     opts.omg = 1e20; end
% if ~isfield(opts, 't2')  
%     opts.t2 = 1e20; end
% if ~isfield(opts, 'iter2')
%     opts.iter2 = 1000;  end

model.t = opts.t1;
model.iter = opts.iter1;
model.tol = opts.tol;
out1 = LP_admm_primal(model);
% out1.t = opts.t2;
% out1.omg = opts.omg;
% out2 = LP_smpf(out1);

out.pi = out1.X;
out.objval = out1.objval;
out.vltcst = out1.vltcst;
out.time = out1.time;
out.iter = out1.iter;
% out.pi = out2.X;
% out.objval = out2.objval;
% out.vltcst = out2.vltcst;
% out.time = out2.time + out1.time;
% out.iter = opts.iter1 + opts.iter2;

end

