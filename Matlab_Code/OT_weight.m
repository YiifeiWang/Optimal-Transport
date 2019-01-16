function [w]=OT_weight(n,opts)
	%-------------------------------------------------------------------------
	% This program construct the weight matrix for optimal transport
	% 
	% Input:
	%	      n --- dimension of resolution
	%      opts --- the option structure with fields:
	%			    p   Wasserstein-p
	%				C   mutiply constant
	%
	% Output:
	%         w --- the weight matrix
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------
	if nargin < 2; opts = [];    end
	if ~isfield (opts,'p');    		opts.p           = 2;   end
	if ~isfield (opts,'C');    		opts.C           = 1;   end
	p = opts.p;
	C = opts.C;

	piece1 = reshape([1:n],1,n);
	aux1   = zeros(n,1);
	block1 = piece1+aux1;
	row1 = reshape(block1,1,n^2);
	col1 = row1';
	mat1 = abs(row1-col1).^p;
	piece2 = reshape([1:n],n,1);
	aux2 = zeros(1,n);
	block2 = piece2+aux2;
	row2 = reshape(block2,1,n^2);
	col2 = row2';
	mat2 = abs(row2-col2).^p;
	w = C*(mat1+mat2);
end