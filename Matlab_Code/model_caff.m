function [model] = model_caff(size_)
	%-------------------------------------------------------------------------
	% This program returns model for caff
	% 
	% Input:
	%   size_ ---  size of the transportation matrix
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

	x = 2 * rand(size_, 2) - 1;
	x = x(x(:, 1).^2 + x(:, 2).^2<=1, :);
	x1 = x(x(:, 1)<=0, :) - [2, 0];
	x2 = x(x(:, 1)>0, :) + [2, 0];
	size_ = length(x);
	model.source = x;
	model.target = [x1; x2];
	distance_ = zeros(size_, size_);
	for i = 1:size_
	    for j = 1:size_
	        distance_(i, j) = norm(model.source(j, :)-model.target(i, :), 2)^2;
	    end
	end
	distance_ = reshape(distance_, size_^2, 1);
	model.obj = distance_;
	model.cst = OT_cst(size_, size_);%ones(2*size, 1);
	model.m = size_;
	model.n = size_;

end

