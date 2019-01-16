function [out] = LP_mosek(model, opts)
	%-------------------------------------------------------------------------
	% This program calls mosek to solve the standard LP probelm
	% from 'model'
	% 
	% Input:
	%   opts --- option structure with fields:
	% 			 method      0: primal simplex
	% 						 1: dual simplex
	% 						 2: iterior point
	% Output:
	%	 out --- output structure with fields:
	%			 --- general ---
	%			 pi			 optimal solution (size:m*n)
	%			 objval		 object value
	%			 vltcst		 violation of constraints
	%			 time		 time elapsed
	%            --- moesek ---
	%			 result      result from mosek
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------

	if nargin <2; opts = []; end
	if ~isfield(opts, 'method');			opts.method = 0; end

	% copy parameter
	method = opts.method;

	out = [];

	m = model.m;
	n = model.n;
	a = OT_matrix(m,n);
    c = reshape(model.obj,m*n,1);
    cst = model.cst;
    blc = model.cst;
    buc = model.cst;
    blx = zeros([model.m*model.n, 1]);
    bux = [];

	param = [];

	switch method
		case 0
			param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_PRIMAL_SIMPLEX';
		case 1
			param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_DUAL_SIMPLEX';
		case 2
			param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
		otherwise
			error('wrong input for opts.method')
	end

	tic;
	result = msklpopt(c,a,blc,buc,blx,bux,param);
	time = toc;

	switch method
		case {0,1}
			x0 = result.sol.bas.xx;
		case 2
			x0 = result.sol.itr.xx;
	end

	out.result = result;

	out.objval = c'*x0;
	out.vltcst = norm(a*x0-cst,1);
	out.time = time;
	out.pi = reshape(x0,m,n);
end