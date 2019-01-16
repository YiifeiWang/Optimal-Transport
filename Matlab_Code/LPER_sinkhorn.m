function [out] = LPER_sinkhorn(model,opts)
	%-------------------------------------------------------------------------
	% This program implements sinkhorn algorithm for LP problem with entropic 
	% regularization.
	% 
	% Input:
	%     model --- the LP model structure with fields:
	%				--- general ---
	%               m, n   		dimension of rows and cols
	%				obj    		matrix C
	%               cst    		constraints
	%
	%	   opts --- the option structure with fields:
	%				epsilon   	parameter of entropic regularization
	%				itermax   	maximum number of iteration
	%				add_coup  	whether to do additional coupling in the end of the algorithm
	%				tqdm		whether to use tqdm
	%               
	%
	% Output:
	%       out --- the output structure with fields:
	%				--- general ---
	%				pi     		optimal solution
	%				objval 		objective value
	%               vltcst    	violation of constraints
	%				time   		time elapsed
	%				--- sinkhorn ---
	%				maxiter maximum number of iteration
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------

	tic;
	if nargin<2; opts=[]; end
	if ~isfield(opts, 'epsilon');  opts.epsilon = 1;   	  end
	if ~isfield(opts, 'maxiter');  opts.maxiter = 1e3;    end
	if ~isfield(opts, 'add_coup'); opts.add_coup = 0;     end
	if ~isfield(opts, 'tqdm');	   opts.tqdm = 0; 		  end

	% copy parameter
	epsilon = opts.epsilon;
	maxiter = opts.maxiter;
	add_coup = opts.add_coup;
	tqdm = opts.tqdm;

	m = model.m;
	n = model.n;
	C = reshape(model.obj,m,n);
	C_nrm = norm(C,inf);
	cst = model.cst;
	mu = cst(1:m);
	nu = cst(m+1:m+n);
	logmu = log(mu);
	lognu = log(nu);

	K = exp(-C/epsilon);
	if any(isnan(K))
	    error('NaN encountered.')
	end

	b = ones(n,1);

	iter = 1;
	while iter <= maxiter
		a = mu./(K*b);
		b = nu./(K'*a);
	    iter = iter+1;
	    if tqdm && mod(iter,maxiter/10)==0
	    	fprintf(['(',mat2str(iter),'/',mat2str(maxiter),')']);toc;
	    end
	end
	time = toc;

	if any(isnan(a)) || any(isnan(b))
	    error('NaN encountered.')
	end

	if add_coup
		a = a.*min(mu./(a.*(K*b)),1);
		b = b.*min(nu./(b.*(K'*a)),1);
		delta1 = mu-a.*(K*b);
		delta2 = nu-b.*(K'*a);
		Pi = a.*K.*b'+delta1*delta2'/norm(delta1,1);
	else
		Pi = (a.*K).*b';
	end

	out.pi = Pi;
	out.vltcst = norm([sum(Pi, 1)-nu', sum(Pi, 2)'-mu'],1);
	out.objval = sum(sum(C.*Pi));
	out.time = time;


end
