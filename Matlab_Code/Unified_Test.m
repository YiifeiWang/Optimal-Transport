function [out] = Unified_Test(opts)
	%-------------------------------------------------------------------------
	% This is a unified test program to test various algorihtms on LP problem
	% based on different dataset
	% 
	% Input:
	%   opts --- option structure with fields:
	%			 --- dataset ---
	% 			 data 	      0: RGD
	% 						  1: DOTmark
	% 						  2: caff
	%						  3: ellipse
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
    %            size       the number of source/target points
    %
    %			 --- options for gaussian mixture model ---
    %			 N   		 the transport matrix is N*N
    %			 gmmopt      options structure with fields:
    %						 mu1,mu2       means
    %						 sigma1,sigma2 covariances
    %						 rho1,rho2     component proportion
    %
    %			 --- algorithm ---
	%			 alg          0: mosek
	%						  1: gurobi
	%						  2: ADMM-primal
	%						  3: ADMM-dual
	%						  4: sinkhorn
	%						  5: sinkhorn with numerical stability
	%						  6: Bregman Admm
	% 			 algopt       option structure with different algorihms
	%						 
	% Output:
	%    out --- output structure with fields:
	%			 --- general ---
	%			 x 			 optimal solution
	%			 objval		 object value
	%			 gpi		 violation of constraints
	%			 time		 time elapsed
	%			 --- caff and ellipse ---
	%            size
	%			 source
	%			 target
	%            --- mosek ---
	%			 result      result from mosek
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12

	if nargin<1; opts=[]; end
	

    if ~isfield(opts, 'alg');				opts.alg 	  = 0; end
    if ~isfield(opts, 'algopt');			opts.algopt   = []; end


    % copy parameter
	alg = opts.alg;
	algopt = opts.algopt;

	% the default setting for opts is given in model_unified.m
	model = model_unified(opts);

	switch alg
		case 0
			[out] = LP_mosek(model,algopt);
		case 1
			[out] = LP_gurobi(model,algopt);
		case 2
			[out] = LP_admm_primal_test(model, algopt);
		case 3
			[out] = LP_admm_dual_test(model, algopt);
        case 4
            [out] = LP_admm_split_test(model, algopt);
		case 5
			[out] = LPER_sinkhorn(model, algopt);
		case 6
			[out] = LPER_shnsc(model, algopt);
		case 7
			[out] = LP_BADMM(model,algopt);
        case 8
            [out] = LPER_admm_test(model, algopt);
		otherwise
			error('wrong input for opts.alg')


end