clear all;

test_times = 1;
num_data = 10;
num_method = 2;

eps = [1e-2, 1e-4, 1e-6];
[~, m] = size(eps);
class_num = [1,2,3,4,5,6,7,8,9,10];

shnsc_maxiter = [4000];
shnsc_tol = [5e-7];
admm_tol = [2e-7];

time = zeros(num_data, 1+m*num_method);
iter = zeros(num_data, 1+m*num_method);
objval = zeros(num_data, 1+m*num_method);
entval = zeros(num_data, 1+m*num_method);
vltcst = zeros(num_data, 1+m*num_method);

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
    
    for j = 1:m
    % call sinkhorn
	opts.alg = 6;
    opts.algopt.epsilon = eps(j);
    opts.algopt.epsilon0 = eps(j)*10^4;
    opts.algopt.maxiter = shnsc_maxiter(1);
	out = Unified_Test(opts);
	time(i, 2*j) = out.time;
    iter(i, 2*j) = out.iter;
	objval(i, 2*j) = out.objval;
    entval(i, 2*j) = out.entval;
	vltcst(i, 2*j) = out.vltcst;
    
    % call admm
	opts.alg = 8;
    opts.algopt.tol = admm_tol(1);
    opts.algopt.eps = eps(j);
	out = Unified_Test(opts);
	time(i, 2*j+1) = out.time;
    iter(i, 2*j+1) = out.iter;
	objval(i, 2*j+1) = out.objval;
    entval(i, 2*j+1) = out.entval;
	vltcst(i, 2*j+1) = out.vltcst;
    
    end
    
end

% for i=1:num_data
% 	fprintf('RGD-%d\n',ms(i));
% 	fprintf('      mosek: cpu: %5.2f, objval: %3.2e, vltcst: %3.2e\n', time(i,1), objval(i,1), vltcst(i,1));
%     fprintf('admm-primal: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,2),(objval(i,2)-objval(i,1))/abs(objval(i,1)), vltcst(i,2), iter(i,2));
% 	fprintf('  admm-dual: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,3),(objval(i,3)-objval(i,1))/abs(objval(i,1)), vltcst(i,3), iter(i,3));
%     fprintf(' admm-split: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,4),(objval(i,4)-objval(i,1))/abs(objval(i,1)), vltcst(i,4), iter(i,4));
%     fprintf('      BADMM: cpu: %5.2f, objval-to-mosek: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,5),(objval(i,5)-objval(i,1))/abs(objval(i,1)), vltcst(i,5), iter(i,5));
% end

fprintf('latex\n');
for i=1:num_data
	fprintf('DOTmark-%d\n',i);
	fprintf('& time(s) ');
	for j=1:m*num_method+1
		fprintf('& %3.2f ',time(i,j));
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('& iter ');
	for j=1:m*num_method+1
		fprintf('& %d ',iter(i,j));
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('& objval ');
	for j=1:m*num_method+1
		if j==1
			fprintf('& %3.2e ',objval(i,j));
		else
			fprintf('& %3.2e ',(objval(i,j)-objval(i,1))/objval(i,1));
		end
	end
	fprintf('\\\\\\cline{2-8}\n');
    fprintf('& entval ');
	for j=1:m*num_method+1
		fprintf('& %3.2e ',entval(i,j));
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('& vltcst ');
	for j=1:m*num_method+1
		fprintf('& %3.2e ',vltcst(i,j));
	end
	fprintf('\\\\\\hline\n');
end
filename = strcat('./data/results_DOTmark_er',num2str(test_times));
save(filename, 'time', 'iter', 'objval', 'entval', 'vltcst');