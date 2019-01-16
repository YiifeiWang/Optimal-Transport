function [model] = model_RGD(m, n)
	%-------------------------------------------------------------------------
	% This program returns model for randomly generated data (RGD) 
	% 
	% Input:
	%     m, n --- the number of rows and cols of the transport matrix
	%
	% Output:
	%    model --- the OT model structure with fields:
	%              m, n   dimension of rows and cols
	%			   obj    matrix C
	%              cst    constraints
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------

	model = [];
	model.m = m;
	model.n = n;
	C = normrnd(0, 1, m*n, 1);
	C = C-min(C);
	model.obj = reshape(C,m,n);
	model.cst = OT_cst(m, n);

end