function [cst]  = OT_cst(m, n)
	%-------------------------------------------------------------------------
	% This program generates the constraints of LP problem
	% 
	% Input:
	%      m, n --- dimension of rows and cols of the transport matrix
	%
	% Output:
	%       cst --- a (m+n)*1 vector as constraints
	%               
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------
	
	T = abs(normrnd(0, 1, m, n));
	cst = zeros([m+n,1]);
	Tsum = sum(T(:));
	cst(1:m) = sum(T,2)/Tsum;
	cst(m+1:m+n) = sum(T,1)'/Tsum;
end