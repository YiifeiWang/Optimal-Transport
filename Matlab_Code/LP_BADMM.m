function [out] = LP_BADMM(model,opts)
	%-------------------------------------------------------------------------
	% This program implements Bregman ADMM for LP problem.
	% 
	% Input:
	%     model --- the LP model structure with fields:
	%				--- general ---
	%               m, n   		dimension of rows and cols
	%				obj    		matrix C
	%               cst    		constraints
	%
	%	   opts --- the option structure with fields:
	%				rho   		parameter of ADMM
	%				maxiter   	maximum number of iteration
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
	%				--- BADMM ---
	%				maxiter maximum number of iteration
	%				iter 		number of iteration
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------
	tic;
	if nargin<2; opts=[]; end
	if ~isfield(opts, 'maxiter');  	opts.maxiter = 2e4; end
	if ~isfield(opts, 'tqdm');	   	opts.tqdm = 0; 		end
	if ~isfield(opts, 'rho');  		opts.rho = 10; 		end
	if ~isfield(opts, 'tol');		opts.tol = 1e-6; 	end

	% copy parameter
	maxiter = opts.maxiter;
	tqdm = opts.tqdm;
	rho = opts.rho;
	tol = opts.tol;

	m = model.m;
	n = model.n;
	C = reshape(model.obj,m,n);
	cst = model.cst;
	mu = max(cst(1:m),eps);
	nu = max(cst(m+1:m+n),eps);

	Cmin = min(C(:));
	C = C-Cmin;
	Cmean = mean(C(:));

	Chat = exp(-C/rho/Cmean);
	C = C+Cmin;

	Pi2 = mu*nu';
	Sigma = zeros(m,n);

	iter = 0;
	time = 0;
	while iter < maxiter
		eSigma = exp(Sigma);
		Pi1 = Pi2.*Chat./eSigma;
		Pi1 = bsxfun(@times,Pi1,mu./sum(Pi1, 2));
		Pi2 = Pi1 .* eSigma;
		Pi2 = bsxfun(@times,Pi2,nu'./sum(Pi2, 1));
		Sigma = Sigma+Pi1-Pi2;
	    iter = iter+1;
	    time = time+toc;
	    if tqdm && mod(iter,maxiter/20)==0
	    	fprintf(['(',mat2str(iter),'/',mat2str(maxiter),'), time elapsed: %.2e\n'],time);
	    	Pi = 0.5*(Pi1+Pi2);
	    	fprintf('%.2e %.2e %.2e\n',norm([sum(Pi, 1)-nu', sum(Pi, 2)'-mu'],1),sum(sum(C.*Pi)),norm(Pi1-Pi2,1));
	    end
	    Pi = 0.5*(Pi1+Pi2);
	    if iter>500 && norm(Pi1-Pi2,1)<tol
	    	tic;break;
	    end
	    tic;
	end
	Pi = 0.5*(Pi1+Pi2);
	time = time+toc;

	out.pi = Pi;
	out.vltcst = norm([sum(Pi, 1)-nu', sum(Pi, 2)'-mu'],1);
	out.objval = sum(sum(C.*Pi));
	out.time = time;
	out.iter = iter;
end
