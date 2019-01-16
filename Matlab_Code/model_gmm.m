function [model] = model_ellipse(N,opts)
	%-------------------------------------------------------------------------
	% This program returns model for ellipse 
	% 
	% Input:
	%      N  ---  the size of the transport matrix
	%
	% Output:
	%   model ---  the OT model
	%              m, n   dimension of rows and cols
	%			   obj    matrix C
	%              cst    constraints
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------

	mu1 = opts.mu1;
	mu2 = opts.mu2;
	sigma1 = opts.sigma1;
	sigma2 = opts.sigma2;
	rho1 = opts.rho1;
	rho2 = opts.rho2;
	model = [];
	p1 = gmdistribution(mu1, sigma1, rho1);
	p2 = gmdistribution(mu2, sigma2, rho2);
	index = [0:1/(N-1):1]';
	mu = p1.pdf(index) + 0.1; 
	mu = mu/sum(mu);
	nu = p2.pdf(index) + 0.1; 
	nu = nu/sum(nu);
	cst = [mu;nu];
	C=pdist2(index, index).^2;

	model.obj = reshape(C,N^2,1);
	model.cst = cst;
	model.m = N;
    model.n = N;


end

