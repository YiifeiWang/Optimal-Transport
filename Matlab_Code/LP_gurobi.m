function [out] = LP_gurobi(model,opts)
	%-------------------------------------------------------------------------
	% This program calls gorubi to solve the standard LP probelm
	% from 'model'
	% 
	% Input:
	%   opts --- option structure with fields:
	% 			 method     -1: automatic
    % 						 0: primal simplex
	% 						 1: dual simplex
	% 						 2: barrier
	% 						 3: concurrent
	% 						 4: deterministic concurrent
	% 						 5: deterministic concurrent simplex
	%
	% Output:
	%	 out --- output structure with fields:
	%			 --- general ---
	%			 pi			 optimal solution (size:m*n)
	%			 objval		 object value
	%			 vltcst		 violation of constraints
	%			 time		 time elapsed
	%            --- gurobi ---
	%			 result      result from gurobi
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------
	tic;
	if nargin <2; opts = []; end
	% options for gorubi
	if ~isfield(opts, 'method');			opts.method = 0; end

	% copy parameter
	method = opts.method;
	m = model.m;
	n = model.n;
	c = reshape(model.obj,m*n,1);
	a = OT_matrix(m,n);

	out = [];

	g_model.A = a;
	g_model.obj = c;
	g_model.rhs = model.cst;
	g_model.sense = '=';
	g_model.vtype = 'C';
	g_model.modelsense = 'min';

	params.Method = method;
	params.outputflag = 1;

	result = gurobi(g_model, params);
	%disp(result);
	x0 = result.x;
	time = toc;

	out.result = result;
	out.objval = c'*x0;
	out.vltcst = norm(a*x0-model.cst,1);
	out.pi = reshape(x0,m,n);
	out.time = time;
end