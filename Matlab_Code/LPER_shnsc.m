function [out] = LPER_shnsc(model,opts)
	%-------------------------------------------------------------------------
	% This program implements our proposed sinkhorn algorithm with numerical stabilty 
	% and continuation strategy for LP problem with entropic regularization.
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
	%				itermax   	maximum number of sub-iteration
	%				epsPow		decay of epsilon
	%				epsilon0	start value for epsilon
	%				add_coup  	whether to do additional coupling in the end of the algorithm
	%               
	%
	% Output:
	%       out --- the output structure with fields:
	%				--- general ---
	%				pi    		optimal solution
	%				objval 		objective value
	%               vltcst    	violation of constraints
	%				time   		time elapsed
	%				iter        number of iteration
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2019/1
	%%-------------------------------------------------------------------------

	tic;
	if nargin<2; opts=[]; end
	if ~isfield(opts, 'epsilon');  opts.epsilon = 1e-10;   end
	if ~isfield(opts, 'epsilon0'); opts.epsilon0 = 1e-1;  end
	if ~isfield(opts, 'epsPow');   opts.epsPow = 0.1;     end
	if ~isfield(opts, 'maxiter');  opts.maxiter = 1e3;    end
	if ~isfield(opts, 'add_coup'); opts.add_coup = 1;     end
	if ~isfield(opts, 'tqdm');	   opts.tqdm = 0; 		  end

	% copy parameter
	epsilon = opts.epsilon;
	epsilon0 = opts.epsilon0;
	epsPow = opts.epsPow;
	maxiter = opts.maxiter;
	add_coup = opts.add_coup;
	tqdm = opts.tqdm;

	m = model.m;
	n = model.n;
	C = reshape(model.obj,m,n);
	cst = model.cst;
	mu = cst(1:m);
	nu = cst(m+1:m+n);
	logmu = log(mu);
	lognu = log(nu);

	g = zeros(n,1);

	epsilonk = epsilon0;
	k = 0;
	total_iter = 0;
	while epsilonk>epsilon
		for iter = 1:maxiter
			fhat = max(-C+g',[],2);
			f = epsilonk*(logmu-log(sum(exp((-C+g'-fhat)/epsilonk),2)))-fhat;
			ghat = max(-C'+f',[],2);
			g = epsilonk*(lognu-log(sum(exp((-C'+f'-ghat)/epsilonk),2)))-ghat;
		end
		if tqdm 
	    	fprintf('(%.2e/%.2e)',epsilonk,epsilon);toc;
	    end
		epsilonk = max(epsilonk*epsPow,epsilon);
		k = k+1;
	end
	for iter = 1:maxiter
		fhat = max(-C+g',[],2);
		f = epsilon*(logmu-log(sum(exp((-C+g'-fhat)/epsilon),2)))-fhat;
		ghat = max(-C'+f',[],2);
		g = epsilon*(lognu-log(sum(exp((-C'+f'-ghat)/epsilon),2)))-ghat;
	end

	if add_coup
		fhat = max(-C+g',[],2);
		f = f+min(epsilon*(logmu-log(sum(exp((-C+g'-fhat)/epsilon),2)))-fhat-f,0);
		ghat = max(-C'+f',[],2);
		g = g+min(epsilon*(lognu-log(sum(exp((-C'+f'-ghat)/epsilon),2)))-ghat-g,0);
		fhat = max(-C+g',[],2);

		Pip = exp((f-C+g')/epsilon);
		delta1 = mu-sum(Pip,2);
		delta1 = delta1/norm(delta1,1);
		delta2 = nu-sum(Pip',2);
		Pi = Pip+delta1*delta2';
	else
		Pi = exp((f-C+g')/epsilon);
	end

	time = toc;

	if any(isnan(Pi))
	    error('NaN encountered.')
	end

	out.pi = Pi;
	out.vltcst = norm([sum(Pi, 1)-nu', sum(Pi, 2)'-mu'],1);
	out.objval = sum(sum(C.*Pi));
	out.entval = out.objval+epsilon*sum(arrayfun(@discrete_entropy,Pi(:)));
	out.time = time;
	out.iter = k*maxiter;

	function Y = discrete_entropy(p)
		if p < eps
			Y = 0;
		else
			Y = p*(log(p)-1);
		end
	end

end
