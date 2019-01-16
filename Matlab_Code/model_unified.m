function [model] = model_unified(opts)
	%-------------------------------------------------------------------------
	% This is a unified program to generate model
	% based on different dataset
	% 
	% Input:
	%   opts --- option structure with fields:
	%			 --- dataset ---
	% 			 data 	      0: RGD
	% 						  1: DOTmark
	% 						  2: caff
	%						  3: ellipse
	%						  4: gaussian mixture model
	%
    %			 seed		  random seed
	% 			 --- options for RGD ---
	%            m, n         the number of rows and cols of the transport matrix
	%            
	%            --- options for DOTmark ---
	%            nd 		  resolution of image is nd * nd
	%			 num  		  class of DOTmark
	%					   	  1: CauchyDensity
	%						  2: GRFmoderate
	%						  3: WhiteNoise
	%		     image1       number for image1 (range: 1 to 10)
	%			 image2		  number for image2 (range: 1 to 10)
    %
    %            --- options for ellipse/caff ---
    %            size        the number of source/target points
    %
    %			 --- options for gaussian mixture model ---
    %			 N   		 the transport matrix is N*N
    %			 gmmopt      options structure with fields:
    %						 mu1,mu2       means
    %						 sigma1,sigma2 covariances
    %						 rho1,rho2     component proportion
    %
	% Output:
	%    model --- the LP model structure with fields:
	%				--- general ---
	%               m, n   dimension of rows and cols
	%				obj    matrix C
	%               cst    constraints
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12

	if nargin<1; opts=[]; end;
	if ~isfield(opts, 'data');    			opts.data     = 0; end
	if ~isfield(opts, 'seed');    			opts.seed     = 0; end
	% options for RGD
	if ~isfield(opts, 'm');					opts.m        = 32; end
	if ~isfield(opts, 'n');					opts.n        = 32; end
	% options for DOTMARK
	if ~isfield(opts, 'nd');				opts.nd       = 32; end
	if ~isfield(opts, 'num');				opts.num      = 1; end
	if ~isfield(opts, 'image1');			opts.image1   = 1; end
	if ~isfield(opts, 'image2');			opts.image2   = 2; end
    % options for ellipse/caff
    if ~isfield(opts, 'size');              opts.size     = 20; end
    % options for gaussian mixture model
    if ~isfield(opts, 'N');              	opts.N     = 128; end
    if ~isfield(opts, 'gmmopt');			opts.gmmopt   = []; end
    gmmopt = opts.gmmopt;
    if ~isfield(gmmopt, 'mu1');				gmmopt.mu1    = [0.3;0.5]; end
    if ~isfield(gmmopt, 'mu2');				gmmopt.mu2    = [0.6;0.7]; end
    if ~isfield(gmmopt, 'sigma1');			gmmopt.sigma1 = reshape([0.05, 0.03].^2, 1, 1, 2); end
    if ~isfield(gmmopt, 'sigma2');			gmmopt.sigma2 = reshape([0.03, 0.05].^2, 1, 1, 2); end
    if ~isfield(gmmopt, 'rho1');			gmmopt.rho1   = [0.5,0.5]; end
    if ~isfield(gmmopt, 'rho2');			gmmopt.rho2   = [0.6,0.4]; end


    if ~isfield(opts, 'alg');				opts.alg 	  = 0; end
    if ~isfield(opts, 'algopt');			opts.algopt   = []; end


    % copy parameter
	data = opts.data;
	alg = opts.alg;
	algopt = opts.algopt;
	size_ = opts.size;
	N = opts.N;

	out = [];

	% set random seed
	rng(opts.seed);
	switch data
		case 0
			model = model_RGD(opts.m,opts.n);
		case 1
			model = model_DOTmark(opts);
        case 2
            model = model_ellipse(size_);
        case 3
            model = model_caff(size_);
        case 4
        	model = model_gmm(N,gmmopt);
		otherwise
			error('wrong input for opts.data')
	end

end