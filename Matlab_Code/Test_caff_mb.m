clear all;

test_times = 1;
num_method = 6;


ms = [128,256,512,1024];%
num_data = 4;

time = zeros(num_data, num_method);
iter = zeros(num_data, num_method);
objval = zeros(num_data, num_method);
vltcst = zeros(num_data, num_method);

% RGD
for i=1:num_data
	opts = [];
	opts.data = 3;
	opts.size = ms(i);

	% call mosek prim
	opts.alg = 0;
	algopt=[];
	algopt.method = 0;
	opts.algopt=algopt;
	out = Unified_Test(opts);
	time(i,1) = out.time;
	objval(i,1) = out.objval;
	vltcst(i,1) = out.vltcst;

	% call mosek dual
	opts.alg = 0;
	algopt.method = 1;
	opts.algopt=algopt;
	out = Unified_Test(opts);
	time(i,2) = out.time;
	objval(i,2) = out.objval;
	vltcst(i,2) = out.vltcst;

	% call mosek int
	algopt.method = 2;
	opts.algopt=algopt;	
	out = Unified_Test(opts);
	time(i,3) = out.time;
	objval(i,3) = out.objval;
	vltcst(i,3) = out.vltcst;

	% call gurobi prim
	opts.alg = 1;
	algopt.method = 0;
	opts.algopt=algopt;
	out = Unified_Test(opts);
	time(i,4) = out.time;
	objval(i,4) = out.objval;
	vltcst(i,4) = out.vltcst;
	iter(i,4) = out.result.itercount;

	% call gurobi dual
	algopt.method = 1;
	opts.algopt=algopt;
	out = Unified_Test(opts);
	time(i,5) = out.time;
	objval(i,5) = out.objval;
	vltcst(i,5) = out.vltcst;
	iter(i,5) = out.result.itercount;

	% call gurobi int
	algopt.method = 2;
	opts.algopt=algopt;
	out = Unified_Test(opts);
	time(i,6) = out.time;
	objval(i,6) = out.objval;
	vltcst(i,6) = out.vltcst;
	iter(i,6) = out.result.baritercount;
	
end

for i=1:num_data
	fprintf('Caff-%d\n',ms(i));
	fprintf('prim(M): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e\n', time(i,1), objval(i,1), vltcst(i,1));
	fprintf('dual(M): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,2),(objval(i,2)-objval(i,1))/abs(objval(i,1)), vltcst(i,2), iter(i,2));
	fprintf(' int(M): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,3),(objval(i,3)-objval(i,1))/abs(objval(i,1)), vltcst(i,3), iter(i,3));
	fprintf('prim(G): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,4),(objval(i,4)-objval(i,1))/abs(objval(i,1)), vltcst(i,4), iter(i,4));
	fprintf('dual(G): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,5),(objval(i,5)-objval(i,1))/abs(objval(i,1)), vltcst(i,5), iter(i,5));
	fprintf(' int(G): cpu: %5.2f, objval: %3.2e, vltcst: %3.2e, iter: %d\n', time(i,6),(objval(i,6)-objval(i,1))/abs(objval(i,1)), vltcst(i,6), iter(i,6));
end

fprintf('latex');
for i=1:num_data
	fprintf('Caff-%d\n',ms(i));
	fprintf('&time(s)');
	for j=1:num_method
		fprintf('& %3.2f ',time(i,j));
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('&iter');
	for j=1:num_method
		fprintf('& %d ',iter(i,j));
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('&objval');
	for j=1:num_method
		if j==1
			fprintf('& %3.2e ',objval(i,j));
		else
			fprintf('& %3.2e ',(objval(i,j)-objval(i,1))/objval(i,1));
		end
	end
	fprintf('\\\\\\cline{2-8}\n');
	fprintf('&vltcst');
	for j=1:num_method
		fprintf('& %3.2e ',vltcst(i,j));
	end
	fprintf('\\\\\\hline\n');
end
filename = strcat('./data/results_caff_mg_tt',num2str(test_times));
save(filename,'time','iter','objval','vltcst');