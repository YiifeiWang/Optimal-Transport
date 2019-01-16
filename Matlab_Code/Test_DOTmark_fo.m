clear all;

test_times = 1;
num_data = 10;
num_method = 5;
class_num = [1,2,3,4,5,6,7,8,9,10];

BADMM_tol = [4e-6];
ADMM_primal_tol = [5e-7];
ADMM_dual_tol = [2e-4];
ADMM_split_tol = [4e-7];

time = zeros(num_data, num_method);
iter = zeros(num_data, num_method);
objval = zeros(num_data, num_method);
vltcst = zeros(num_data, num_method);

% DOTmark
for i=1:num_data
	opts = [];
	opts.data = 1;
	opts.num = class_num(i);
	opts.nd = 32;

    % call mosek
	opts.alg = 0;
	out = Unified_Test(opts);
	time(i,1) = out.time;
    iter(i,1) = 0;
	objval(i,1) = out.objval;
	vltcst(i,1) = out.vltcst;
    
	 % call admm primal
	opts.alg = 2;
    opts.algopt.tol = ADMM_primal_tol(1);
	out = Unified_Test(opts);
	time(i,2) = out.time;
    iter(i,2) = out.iter;
	objval(i,2) = out.objval;
	vltcst(i,2) = out.vltcst;
    
    % call admm dual
	opts.alg = 3;
    opts.algopt.tol = ADMM_dual_tol(1);
	out = Unified_Test(opts);
	time(i,3) = out.time;
    iter(i,3) = out.iter;
	objval(i,3) = out.objval;
	vltcst(i,3) = out.vltcst;
    
    % call admm split
	opts.alg = 4;
    opts.algopt.tol = ADMM_split_tol(1);
	out = Unified_Test(opts);
	time(i,4) = out.time;
    iter(i,4) = out.iter;
	objval(i,4) = out.objval;
	vltcst(i,4) = out.vltcst;

	% Bregman ADMM
	opts.alg = 7;
	algopt = [];
	algopt.maxiter = 2e4;
	algopt.tol = BADMM_tol(1);
	opts.algopt = algopt;
	out = Unified_Test(opts);
	time(i,5) = out.time;
	iter(i,5) = out.iter;
	objval(i,5) = out.objval;
	vltcst(i,5) = out.vltcst;
	
end

for i=1:num_data
	fprintf('&mosek: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e\n', time(i,1), objval(i,1), vltcst(i,1));
	fprintf('&admm-primal: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,2),(objval(i,2)-objval(i,1))/abs(objval(i,1)), vltcst(i,2), iter(i,2));
	fprintf(' &admm-dual: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,3),(objval(i,3)-objval(i,1))/abs(objval(i,1)), vltcst(i,3), iter(i,3));
	fprintf('&admm-split: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,4),(objval(i,4)-objval(i,1))/abs(objval(i,1)), vltcst(i,4), iter(i,4));
	fprintf('&BADMM: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,5),(objval(i,5)-objval(i,1))/abs(objval(i,1)), vltcst(i,5), iter(i,5));
end

fprintf('latex');
for i=1:num_data
	fprintf('DOTmark-%d\n',i);
	fprintf('&time(s) ');
	for j=1:num_method
		fprintf('& %3.2f ',time(i,j));
	end
	fprintf('\\\\\\cline{2-7}\n');
	fprintf('&iter ');
	for j=1:num_method
		fprintf('& %d ',iter(i,j));
	end
	fprintf('\\\\\\cline{2-7}\n');
	fprintf('&objval ');
	for j=1:num_method
		if j==1
			fprintf('& %3.2e ',objval(i,j));
		else
			fprintf('& %3.2e ',(objval(i,j)-objval(i,1))/objval(i,1));
		end
	end
	fprintf('\\\\\\cline{2-7}\n');
	fprintf('&vltcst ');
	for j=1:num_method
		fprintf('& %3.2e ',vltcst(i,j));
	end
	fprintf('\\\\\\hline\n');
end
%filename = strcat('./data/results_DOTmark_tt',num2str(test_times));
%save(filename,'time','iter','objval','vltcst');