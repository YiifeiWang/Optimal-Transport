function [model] = model_DOTmark(opts)
	%-------------------------------------------------------------------------
	% This program returns model for DOTmark
	% 
	% Input:
	%    opts --- option structure with fields:
	%			  nd 		resolution of image is nd * nd
	%			  num  		class of DOTmark
	%						1: CauchyDensity
	%						2: ClassicImages
	%						3: GRFmoderate
	%						4: GRFrough
	%						5: GRFsmooth
	%						6: LogGRF
	%						7: LogitGRF
	%						8: MicroscopyImages
	%						9: Shapes
	%					   10: WhiteNoise
	%
	%
	%			  image1    number for image1 (range: 1 to 10)
	%			  image2	number for image2 (range: 1 to 10)
	%			  
	%
	% Output:
	%   model ---  the OT model structure with fields:
	%              m, n   dimension of rows and cols
	%			   obj    matrix C
	%              cst    constraints
	%
	% Author: Yifei Wang & Feng Zhu
	% Version 1.1 .... 2018/12
	%%-------------------------------------------------------------------------

	classes = [{'CauchyDensity'}, {'ClassicImages'}, {'GRFmoderate'}, {'GRFrough'}, {'GRFsmooth'}, {'LogGRF'}, {'LogitGRF'}, {'MicroscopyImages'}, {'Shapes'}, {'WhiteNoise'}];
	nd = opts.nd;
	num = opts.num;
	image1 = opts.image1;
	image2 = opts.image2;
	switch image1
		case {1,2,3,4,5,6,7,8,9}
			file1 = ['DOTmark_1/Data/', classes{num}, '/data',mat2str(nd), '_100',mat2str(image1), '.csv'];
		case 10
			file1 = ['DOTmark_1/Data/', classes{num}, '/data',mat2str(nd), '_1010.csv'];
		otherwise
			error('wrong input for image1')
	end
	switch image2
		case {1,2,3,4,5,6,7,8,9}
			file2 = ['DOTmark_1/Data/', classes{num}, '/data',mat2str(nd), '_100',mat2str(image2), '.csv'];
		case 10
			file2 = ['DOTmark_1/Data/', classes{num}, '/data',mat2str(nd), '_1010.csv'];
		otherwise
			error('wrong input for image2')
	end

	model = [];
	W = OT_weight(nd);
	model.obj = reshape(W,nd^4,1);
	model.mat = OT_matrix(nd^2,nd^2);
	f1 = csvread(file1);
	f2 = csvread(file2);
	f1_sum = sum(f1(:));
	cst = zeros(2*nd^2,1);
	cst(1:nd^2) = reshape(f1,nd^2,1);
	cst(nd^2+1:2*nd^2) = reshape(f2,nd^2,1);
	cst = max(cst./f1_sum, eps);
	model.cst = cst;
	model.m = nd^2;
	model.n = nd^2;

end