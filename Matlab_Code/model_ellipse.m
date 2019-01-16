function [model] = model_ellipse(size_)
	%-------------------------------------------------------------------------
	% This program returns model for ellipse 
	% 
	% Input:
	%   m, n  ---  the number of rows and cols of the transport matrix
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

	model = [];
	x = 2 * pi * rand(size_, 1);
	data = [cos(x), sin(x)] + 0.1 * randn(size_, 2);
	model.source = [2, 0.5].*data;
	x = 2 * pi * rand(size_, 1);
	data = [cos(x), sin(x)] + 0.1 * randn(size_, 2);
	model.target = [0.5, 2].*data;
	distance_ = zeros(size_, size_);
	for i = 1:size_
	    for j = 1:size_
	        distance_(i, j) = norm(model.source(j, :)-model.target(i, :), 2)^2;
	    end
	end
	model.obj = distance_;
	model.cst = OT_cst(size_, size_);%ones(2*size, 1);
	model.m = size_;
    model.n = size_;

end

