clear all;

test_times = 1;
num_data = 4;
num_method = 5;

ms = [128, 256, 512, 1024];
BADMM_tol = [1e-6];
ADMM_primal_tol = [5e-7];
ADMM_dual_tol = [2e-4];
ADMM_split_tol = [1e-7];

time = zeros(num_data, num_method);
iter = zeros(num_data, num_method);
objval = zeros(num_data, num_method);
vltcst = zeros(num_data, num_method);

% RGD
for i=1:num_data
	opts = [];
	opts.data = 0;
	opts.m = ms(i);
	opts.n = ms(i);

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
	fprintf('RGD-%d\n',ms(i));
	fprintf('      mosek: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e\n', time(i,1), objval(i,1), vltcst(i,1));
    fprintf('admm-primal: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,2),(objval(i,2)-objval(i,1))/abs(objval(i,1)), vltcst(i,2), iter(i,2));
	fprintf('  admm-dual: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,3),(objval(i,3)-objval(i,1))/abs(objval(i,1)), vltcst(i,3), iter(i,3));
    fprintf(' admm-split: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,4),(objval(i,4)-objval(i,1))/abs(objval(i,1)), vltcst(i,4), iter(i,4));
    fprintf('      BADMM: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,5),(objval(i,5)-objval(i,1))/abs(objval(i,1)), vltcst(i,5), iter(i,5));
end

fprintf('latex');
for i=1:num_data
	fprintf('RGD-%d\n',ms(i));
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
%filename = strcat('./data/results_RGD_tt',num2str(test_times));
%save(filename,'time','iter','objval','vltcst');