function [A] = OT_matrix(m, n)
	%-------------------------------------------------------------------------
	% This program generates the coefficient matrix of constraints
	% 
	% Input:
	%      m, n --- dimension of rows and cols of the transport matrix
	%
	% Output:
	%         A --- a (m+n)*(m*n) sparse matrix as the coefficient of constraints
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------
	x = kron(sparse(diag(ones(1, n))), sparse(ones(1, m)));
	y = kron(sparse(ones(1, n)), sparse(diag(ones(1, m))));
	A = cat(1, y, x);
end

